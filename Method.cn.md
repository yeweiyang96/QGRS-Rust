## 1. BFS 候选扩展算法

BFS（广度优先搜索）候选扩展是 QGRS-Rust 核心搜索逻辑，用于枚举序列中所有合法的 G-quadruplex 结构。算法分为四个阶段：种子生成、广度优先循环、Loop 发现、评分筛选。

### 1.1 种子生成（seed_queue）

**目标**：扫描序列识别所有潜在的 G-run（连续 G 碱基），并在「允许的 tetrad 数量 × 允许的偏移」笛卡尔积上生成初始候选。所有后续 BFS 扩展都建立在这些种子之上。

**实现细节**：
```rust
fn seed_queue(
    cands: &mut VecDeque<G4Candidate>,
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    limits: ScanLimits,
) {
    // 限制 tetrad 数，避免生成超过 max_g_run 或 max_g4_length/4 的种子
    let mut max_tetrads_allowed = limits.max_g_run;
    if limits.max_g4_length >= 4 {
        max_tetrads_allowed = max_tetrads_allowed.min(limits.max_g4_length / 4);
    }
    if max_tetrads_allowed < min_tetrads {
        return;
    }

    for (run_start, run_len) in GRunScanner::new(&seq.normalized, min_tetrads) {
        let max_tetrads_for_run = run_len.min(max_tetrads_allowed);
        let mut tetrads = min_tetrads;
        while tetrads <= max_tetrads_for_run {
            if tetrads * 4 > limits.max_g4_length {
                break; // tetrad 自身已超出整体长度限制
            }
            let max_offset = run_len.saturating_sub(tetrads);
            for offset in 0..=max_offset {
                let start = run_start + offset;
                cands.push_back(G4Candidate::new(seq.clone(), tetrads, start, limits));
            }
            tetrads += 1;
        }
    }
}
```

**关键优化**：
- `GRunScanner::new(&seq.normalized, min_tetrads)` 仍基于 `memchr2(b'g', b'G')`，但现在仅返回长度≥`min_tetrads` 的 run，减少了上层过滤成本；SIMD 扫描相较逐字节检查快约 10×。
- `max_tetrads_allowed = min(max_g_run, max_g4_length/4)` 把 G-run 长度与整体长度约束结合起来，防止生成无法通过后续 `max_length` 检查的冗余种子。
- 当 run 长度超过 `max_g_run` 时，`run_len.min(max_tetrads_allowed)` 会截断可用 tetrad 数，但随后的偏移枚举依旧覆盖整段 G-run，使得“大于最大连续 G 限制”的长 run 只会以受限窗口形式进入 BFS，而不会完全丢失。
- 对每个 run，在有效 tetrad 数范围内再枚举所有偏移，确保诸如 `GGGGG` 这类长 run 会覆盖所有子区间，而不用多次扫描序列。

**示例**（`min_tetrads=2`, `max_g_run=5`, `max_g4_length=40`）：
- 序列片段：`GGGGG...`（run_len=5）
- 有效 tetrad 数：2、3、4、5（受 `max_g4_length` 约束，4×5=20 ≤ 40）
- 对于 tetrad=3，允许的偏移为 `0..=2`，因此会生成 3 个起点：`start`, `start+1`, `start+2`
- chunk 模式额外在窗口层面限制 `start < primary_end`，以避免窗口重叠区重复发射 raw hit，但 `seed_queue` 本身保持与 stream/inline 路径一致

### 1.2 广度优先循环（find_raw_with_sequence）

**核心逻辑**：从队列依次取出候选，检查是否完整（三个 loop 都已分配），若完整则评分收录，否则扩展填充下一个 loop。

**伪代码流程**：
```rust
let mut cands: VecDeque<G4Candidate> = seed_queue;
let mut raw_g4s: Vec<G4> = Vec::new();

while let Some(cand) = cands.pop_front() {
    if cand.complete() {           // 检查 y1、y2、y3 是否都已赋值
        if cand.viable(min_score) { // 通过最小分数 & 最大长度检查
            raw_g4s.push(G4::from_candidate(&cand));
        }
    } else {
        // 扩展候选：为下一个未分配的 loop 枚举所有可能长度
        for expanded in cand.expand(limits) {
            if expanded.partial_length() <= max_length {
                cands.push_back(expanded);
            }
        }
    }
}
```

**实际代码**（`mod.rs` 行 850-863）：
```rust
while let Some(cand) = cands.pop_front() {
    if cand.y3.is_some() {
        if cand.viable(&limits) {
            let g4 = G4::from_candidate(&cand, data.clone());
            raw_g4s.push(g4);
        }
    } else {
        for expanded in cand.expand(data, &limits) {
            cands.push_back(expanded);
        }
    }
}
```

**完整性判断**：
- `y1.is_none()` → 需要填充第一个 loop
- `y2.is_none()` → 需要填充第二个 loop
- `y3.is_some()` → 所有 loop 已填充，候选完整

**可行性检查**（`viable`）：
```rust
impl G4Candidate {
    fn viable(&self, limits: &ScanLimits) -> bool {
        let total_len = self.total_length();
        total_len <= limits.max_length && self.score(limits) >= limits.min_score
    }
}
```

**BFS 保证**：
- **完备性**：所有合法组合都会被枚举（队列确保每个候选都被扩展）
- **正确性**：通过 `partial_length()` 剪枝避免无效扩展（超出 max_length 的候选不入队）
- **顺序无关**：同一坐标的候选可能以不同顺序入队（取决于 seed 顺序），但最终去重确保唯一性

### 1.3 Loop 发现机制（expand）

**目标**：从候选当前状态（已分配部分 loop）向前扫描下一个 G-run，为下一个 loop 枚举所有合法长度。

**实现细节**：
```rust
impl G4Candidate {
    fn expand(&self, data: &SequenceData, limits: &ScanLimits) -> Vec<G4Candidate> {
        let cursor = self.next_cursor_position();  // 当前 tetrad 后的位置
        let loop_lengths = find_loop_lengths_from(cursor, data, limits, self);
        
        loop_lengths.into_iter()
            .map(|y| {
                let mut new_cand = self.clone();
                new_cand.assign_next_loop(y);  // 填充 y1/y2/y3
                new_cand
            })
            .collect()
    }
}
```

**cursor 计算规则**：
- y1 未分配 → cursor = tetrad1 结束位置（`start + tetrad_len`）
- y2 未分配 → cursor = tetrad1 + y1 + tetrad2 结束位置
- y3 未分配 → cursor = tetrad1 + y1 + tetrad2 + y2 + tetrad3 结束位置

**find_loop_lengths_from 逻辑**：
```rust
fn find_loop_lengths_from(cursor: usize, data: &SequenceData, limits: &ScanLimits, cand: &G4Candidate) -> Vec<usize> {
    let mut lengths = Vec::new();
    let min_len = cand.min_acceptable_loop_length();  // 0 或 1（若已有 loop=0）
    let max_len = limits.max_length - cand.partial_length() - cand.tetrad_len;
    
    // 从 cursor 向前扫描，寻找长度为 tetrad_len 的 G-run
    for y in min_len..=max_len {
        let tetrad_start = cursor + y;
        if is_valid_tetrad_at(tetrad_start, cand.tetrad_len, data) {
            lengths.push(y);
        }
    }
    lengths
}
```

**约束条件**：
1. **最小长度**：`min_acceptable_loop_length()` 返回 0（默认）或 1（已有 loop=0 时，避免连续两个 0-length loop）
2. **最大长度**：确保 `partial_length() + y + tetrad_len ≤ max_length`（避免候选超出长度限制）
3. **tetrad 匹配**：`is_valid_tetrad_at()` 检查 `cursor+y` 处是否有连续 `tetrad_len` 个 G

**示例扩展**：
- 候选状态：`start=212, tetrad_len=3, y1=5` (已分配第一个 loop)
- cursor = 212 + 3 + 5 = 220（第二个 tetrad 应从此开始）
- 假设在 220 和 222 处各找到长度为 3 的 G-run
- 生成两个新候选：
  - `{start=212, tetrad=3, y1=5, y2=0}` (loop2 长度为 0)
  - `{start=212, tetrad=3, y1=5, y2=2}` (loop2 长度为 2)

### 1.4 评分公式（score）

**Legacy 公式**（继承自 C++ 版本，见 `qgrs.h`）：
```rust
impl G4Candidate {
    fn score(&self, limits: &ScanLimits) -> i32 {
        let gmax = (limits.max_length - (self.tetrad_len * 4 + 1)) as f64;
        let gavg = self.loop_variance_mean();
        let bonus = gmax * (self.tetrad_len as f64 - 2.0);
        (gmax - gavg + bonus).floor() as i32
    }
    
    fn loop_variance_mean(&self) -> f64 {
        let y1 = self.y1.unwrap();
        let y2 = self.y2.unwrap();
        let y3 = self.y3.unwrap();
        let d12 = (y1 as i32 - y2 as i32).abs();
        let d23 = (y2 as i32 - y3 as i32).abs();
        let d13 = (y1 as i32 - y3 as i32).abs();
        (d12 + d23 + d13) as f64 / 3.0
    }
}
```

**参数含义**：
- `gmax`：剩余空间（max_length 减去 4 个 tetrad + 1）
- `gavg`：三个 loop 长度差的平均值（衡量均匀性）
- `bonus`：tetrad 数量加成（2-tetrad 无加成，3-tetrad 加 1×gmax）

**评分特性**：
1. **长度惩罚**：总长度越接近 max_length，gmax 越小，分数越低
2. **均匀性奖励**：loop 长度越接近（gavg 越小），分数越高
3. **tetrad 加成**：更多 tetrad（3+）获得额外加成

**典型分数对比**（max_length=45, min_tetrads=2）：
- `{tetrad=3, y1=5, y2=5, y3=5}`：gmax=32, gavg=0, bonus=32 → score=64
- `{tetrad=3, y1=5, y2=2, y3=5}`：gmax=32, gavg=2, bonus=32 → score=62
- `{tetrad=3, y1=5, y2=1, y3=6}`：gmax=32, gavg=2.67, bonus=32 → score=61

### 1.5 与 C++ 版本的对比

**C++ 实现**（`qgrs.cpp`）：
```cpp
void findCandidates(const char* seq, int len) {
    queue<Candidate> q;
    vector<G4> results;
    // 种子生成：简单字符串扫描（无 SIMD）
    for (int i = 0; i < len; ++i) {
        if (seq[i] == 'G' || seq[i] == 'g') {
            int run_len = 1;
            while (i+run_len < len && (seq[i+run_len]=='G' || seq[i+run_len]=='g'))
                ++run_len;
            for (int off = 0; off <= run_len - min_tetrads; ++off)
                q.push(Candidate(i+off, min_tetrads+off));
        }
    }
    // BFS 循环：与 Rust 版本逻辑完全一致
    while (!q.empty()) {
        Candidate c = q.front(); q.pop();
        if (c.complete()) {
            if (c.viable()) results.push_back(c);
        } else {
            for (auto& exp : c.expand()) q.push(exp);
        }
    }
}
```

**关键差异**：
| 维度 | Rust | C++ |
|------|------|-----|
| G-run 扫描 | `memchr2` SIMD (~10x faster) | 逐字节 while 循环 |
| 队列类型 | `VecDeque<G4Candidate>` | `std::queue<Candidate>` |
| 内存管理 | `Arc<Vec<u8>>` 零拷贝 | `std::string` 每次复制 |
| 并行化 | Rayon 自动分块 | 单线程顺序执行 |
| 评分公式 | **完全相同**（逐行翻译） | 原始 legacy 公式 |

**性能测试结果**（`big.txt` 10MB 序列）：
- C++：~45 秒（单线程）
- Rust (stream)：~8 秒（8 核心）
- Rust (mmap)：~6 秒（8 核心 + 零拷贝）

### 1.6 典型示例解析

**场景**：在位置 212 检测到 5 个连续 G
- 序列：`...GGGGGACGTGGGACGTGGG...`
- 参数：`min_tetrads=2, max_g_run=5, max_length=45`

**BFS 扩展过程**：
1. **种子**：生成 4 个候选（offset 0-3 对应 tetrad=2,3,4,5）
2. **第一轮扩展**（填充 y1）：
   - 候选 `{212, tetrad=3}` 在位置 220 找到 G-run
   - 发现可行 y1 值：[5, 6, 7]（对应位置 220, 221, 222 的 tetrad）
3. **第二轮扩展**（填充 y2）：
   - 候选 `{212, tetrad=3, y1=5}` 扩展为：
     - `{..., y1=5, y2=2}` (gscore=19)
     - `{..., y1=5, y2=5}` (gscore=17)
4. **第三轮扩展**（填充 y3）：
   - 候选 `{..., y1=5, y2=2}` 最终扩展为：
     - `{..., y1=5, y2=2, y3=5}` (gscore=19, length=27)
5. **收录**：通过 `viable()` 检查（length≤45, score≥0），加入 `raw_g4s`

**去重结果**：
- 原始产出包含多个 `start=212` 候选（不同 loop 配置）
- `consolidate_g4s` 第一阶段去重保留 gscore=19 版本（HashMap Entry API）
- 最终输出：`start=212, length=27, tetrads=3, y1=5, y2=2, y3=5, gscore=19`

---
## 2. 家族合并（consolidate_g4s）详解
### 阶段1 - 去重:
用 HashMap<DedupKey, G4> 按 (start,end,sequence) 去重

**DedupKey 定义**：
```rust
#[derive(Eq, PartialEq, Hash)]
struct DedupKey {
    start: usize,           // G4 起始位置（1-based）
    end: usize,             // G4 结束位置（exclusive）
    slice: SequenceSlice,   // 包含实际序列字节的切片引用
}
```
`SequenceSlice` 实现了基于字节内容的哈希和相等比较，确保相同坐标+相同序列的候选被识别为重复

**去重逻辑**：
```rust
let mut best_by_key: HashMap<DedupKey, G4> = HashMap::new();
for g in raw_g4s.into_iter() {
    let key = DedupKey::new(&g);
    match best_by_key.entry(key) {
        Entry::Vacant(slot) => slot.insert(g),
        Entry::Occupied(mut slot) => {
            if g.gscore > slot.get().gscore {
                slot.insert(g);  // 替换为更高分版本
            }
        }
    }
}
```
对同 key 的多个候选保留最高 gscore 版本

**为何需要**：BFS 可能产生坐标相同但 loop 配置不同的候选（如 start=212的 `y1=5,y2=2,y3=5 gscore=19` vs `y1=5,y2=1,y3=6 gscore=17`），HashMap 确保只保留最优配置. HashMap 查找是 O(1)，而 C++ 的 set 是 O(log n)


### 阶段2 - 分组：
线性扫描现有家族，用 belongs_in 判断新候选是否与家族任一成员重叠

**重叠判定（overlapped）**：
```rust
fn overlapped(a: &G4, b: &G4) -> bool {
    let a_start = a.start as isize;
    let a_end = (a.start + a.length) as isize;
    let b_start = b.start as isize;
    let b_end = (b.start + b.length) as isize;
    // 四种区间相交情况任一为真即重叠
    (a_start >= b_start && a_start <= b_end) ||  // a起点在b内
    (a_end >= b_start && a_end <= b_end) ||      // a终点在b内
    (b_start >= a_start && b_start <= a_end) ||  // b起点在a内
    (b_end >= a_start && b_end <= a_end)         // b终点在a内
}
```
判定两个 G4 区间 `[start, start+length)` 是否共享任何碱基位置

**贪心归属策略**：
```rust
for g4 in deduped.into_iter() {
    let mut inserted = false;
    for family in &mut families {
        if belongs_in(&g4, family) {  // 与家族任一成员重叠
            family.push(g4.clone());
            inserted = true;
            break;  // 第一个匹配就停止，不再检查后续家族
        }
    }
    if !inserted { 
        families.push(vec![g4]);  // 创建新家族
    }
}
```
- 候选按 `(start, end)` 排序后依次检查现有家族
- 遇到第一个匹配的家族立即加入并停止（单次归属）
- 若与所有家族都不重叠，自成一族
- **传递性连接**：若 A 与 B 重叠、B 与 C 重叠，即使 A 与 C 不直接重叠，三者仍会归入同一家族

**顺序依赖性**：去重后必须按 `(start, end)` 排序，否则同一批候选在 chunk/stream 模式下因迭代顺序不同可能被分到不同家族，导致输出不一致

### 阶段3 - 择优：
每个家族选最高 gscore 成员输出
生物学原理：重叠结构可能代表同一特征，选最佳配置

### 2.1 元数据导出与 `--overlap`

`consolidate_g4s` 现在同时返回 `(Vec<G4>, Vec<(usize, usize)>)`：第一项是家族优胜者，第二项记录每个家族的 `[start,end]` 范围。CLI 的 `--overlap` 标志会在写入主 CSV/Parquet 的同时追加两类调试文件：

1. `{base}.overlap.csv`：使用 `render_csv_results` 直接 dump 排序后的 raw hits；字段与主 CSV 完全一致，方便做 diff。
2. `{base}.family.csv`：通过 `render_family_ranges_csv` 输出 `family_index,start,end`，可迅速定位哪些家族被合并。

流式管线通过 `StreamChromosomeResults` 把 `hits`、`family_ranges` 以及可选的 `raw_hits` 传给回调：

```rust
pub struct StreamChromosomeResults {
    pub hits: Vec<G4>,
    pub family_ranges: Vec<(usize, usize)>,
    pub raw_hits: Option<Vec<G4>>, // 仅在 --overlap 时填充
}
```

`--overlap` 默认为关闭，以避免无谓的 `Vec` clone。当用户启用该选项：

- inline 模式必须提供 `--output`，否则无法推导出调试文件的基名；
- FASTA/mmap 模式为每条染色体写三份文件（主输出 + `.overlap.csv` + `.family.csv`），所有文件都沿用 `sanitize_name`/去重后缀；
- stream 模式在每个染色体完成时即时 flush 附加 CSV，确保内存占用不会随输入规模增长。

这一机制让你在排查 chunk vs stream 或 Rust vs C++ 差异时，可以直接 diff raw hits 或家族范围，而无需修改核心算法。