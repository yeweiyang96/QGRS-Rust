## 1. BFS 候选扩展算法
* 种子生成：使用 SIMD 加速的 GRunScanner 扫描 G-run，为每个符合条件的 run 在所有有效偏移处生成初始候选
* 广度优先循环：从队列取候选 → 若完整且符合条件则收录 → 否则调用 expand() 填充下一个未分配的 loop
* Loop 发现机制：详细解释 find_loop_lengths_from 如何从当前 cursor 位置向前扫描匹配的 tetrad run，以及长度/位置约束
* 评分公式：完整展示 legacy 公式及其组成部分（gmax、gavg、bonus）
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

**为何需要**：BFS 可能产生坐标相同但 loop 配置不同的候选（如 start=212的 `y1=5,y2=2,y3=5 gscore=19` vs `y1=5,y2=1,y3=6 gscore=17`），HashMap 确保只保留最优配置

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