## 1. BFS 候補展開アルゴリズム

BFS（幅優先探索）による候補展開は、QGRS-Rust の中核となる探索ロジックであり、配列中のすべての妥当な G-quadruplex 構造を列挙します。アルゴリズムは 4 つの段階に分かれます。シード生成、幅優先ループ、ループ検出、スコアによるフィルタリングです。

### 1.1 シード生成（`seed_queue`）

**目的**: 配列を走査してすべての潜在的な G-run（連続した G 塩基列）を見つけ、「許可された tetrad 数 x 許可されたオフセット」の直積上で初期候補を生成します。以降の BFS 展開はすべてこれらのシードを起点にします。

**実装詳細**:
```rust
fn seed_queue(
    cands: &mut VecDeque<G4Candidate>,
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    limits: ScanLimits,
) {
    // max_g_run または max_g4_length / 4 を超えるシードを作らないように
    // tetrad 数を制限する
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
                break; // tetrad 自体の長さだけで全体長制限を超える
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

**主な最適化**:
- `GRunScanner::new(&seq.normalized, min_tetrads)` は引き続き `memchr2(b'g', b'G')` を使いますが、長さが `min_tetrads` 以上の run だけを返すようになっており、上位層のフィルタコストを減らしています。SIMD スキャンは 1 バイトずつ調べる方法より約 10 倍高速です。
- `max_tetrads_allowed = min(max_g_run, max_g4_length/4)` により、G-run 長の制約と全体長の制約を結合し、後段の `max_length` 判定を絶対に通らない冗長なシードの生成を防ぎます。
- run 長が `max_g_run` を超える場合でも、`run_len.min(max_tetrads_allowed)` で利用可能な tetrad 数だけを切り詰め、後続のオフセット列挙で run 全体をカバーします。つまり、「最大連続 G 制限」を超える長い run でも、制限付きウィンドウとして BFS に入るため、完全には失われません。
- 各 run について、有効な tetrad 数の範囲ごとにすべてのオフセットを列挙するため、`GGGGG` のような長い run でも配列を何度も再走査することなく、すべての部分区間を網羅できます。

**例**（`min_tetrads=2`, `max_g_run=5`, `max_g4_length=40`）:
- 配列断片: `GGGGG...`（`run_len=5`）
- 有効な tetrad 数: 2、3、4、5（`max_g4_length` により制約される。`4x5=20 <= 40`）
- `tetrad=3` の場合、許可されるオフセットは `0..=2` なので、開始位置は `start`, `start+1`, `start+2` の 3 つ生成されます
- chunk モードでは、重なったウィンドウから raw hit が重複出力されないようにウィンドウ層で追加の `start < primary_end` 制約を掛けますが、`seed_queue` 自体は stream / inline 経路と同じです

### 1.2 幅優先ループ（`find_raw_with_sequence`）

**中核ロジック**: キューから候補を順に取り出します。候補が完全であれば（3 つの loop がすべて割り当て済み）、スコア計算を行って結果に格納します。未完成であれば、次の loop を埋める形で展開します。

**擬似コード**:
```rust
let mut cands: VecDeque<G4Candidate> = seed_queue;
let mut raw_g4s: Vec<G4> = Vec::new();

while let Some(cand) = cands.pop_front() {
    if cand.complete() {            // y1, y2, y3 がすべて割り当て済みか確認
        if cand.viable(min_score) { // 最小スコアと最大長の条件を満たす
            raw_g4s.push(G4::from_candidate(&cand));
        }
    } else {
        // 候補を展開: 次に未割り当ての loop について取り得る長さをすべて列挙
        for expanded in cand.expand(limits) {
            if expanded.partial_length() <= max_length {
                cands.push_back(expanded);
            }
        }
    }
}
```

**実際のコード**（`mod.rs` 850-863 行）:
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

**完全性判定**:
- `y1.is_none()` -> 1 つ目の loop をまだ埋める必要がある
- `y2.is_none()` -> 2 つ目の loop をまだ埋める必要がある
- `y3.is_some()` -> すべての loop が埋まっており、候補は完成している

**妥当性判定**（`viable`）:
```rust
impl G4Candidate {
    fn viable(&self, limits: &ScanLimits) -> bool {
        let total_len = self.total_length();
        total_len <= limits.max_length && self.score(limits) >= limits.min_score
    }
}
```

**BFS が保証すること**:
- **完全性**: キューによってすべての候補が展開されるため、妥当な組み合わせは漏れなく列挙される
- **正しさ**: `partial_length()` による枝刈りで無効な展開を防ぎ、`max_length` を超える候補はキューに戻されない
- **順序非依存性**: 同一座標の候補がシード順により異なる順番でキューに入ることはあるが、最終的な重複排除で一意性が保証される

### 1.3 ループ検出メカニズム（`expand`）

**目的**: 現在の候補状態（すでに一部の loop が割り当て済み）から前方に走査し、次の G-run を見つけて、次の loop の妥当な長さをすべて列挙します。

**実装詳細**:
```rust
impl G4Candidate {
    fn expand(&self, data: &SequenceData, limits: &ScanLimits) -> Vec<G4Candidate> {
        let cursor = self.next_cursor_position();  // 現在の tetrad の直後の位置
        let loop_lengths = find_loop_lengths_from(cursor, data, limits, self);

        loop_lengths.into_iter()
            .map(|y| {
                let mut new_cand = self.clone();
                new_cand.assign_next_loop(y);  // y1 / y2 / y3 を埋める
                new_cand
            })
            .collect()
    }
}
```

**cursor の計算規則**:
- `y1` 未割り当て -> `cursor =` 1 つ目の tetrad の終端位置（`start + tetrad_len`）
- `y2` 未割り当て -> `cursor =` 1 つ目の tetrad + `y1` + 2 つ目の tetrad の終端
- `y3` 未割り当て -> `cursor =` 1 つ目の tetrad + `y1` + 2 つ目の tetrad + `y2` + 3 つ目の tetrad の終端

**`find_loop_lengths_from` のロジック**:
```rust
fn find_loop_lengths_from(cursor: usize, data: &SequenceData, limits: &ScanLimits, cand: &G4Candidate) -> Vec<usize> {
    let mut lengths = Vec::new();
    let min_len = cand.min_acceptable_loop_length();  // 0 または 1（すでに loop=0 がある場合）
    let max_len = limits.max_length - cand.partial_length() - cand.tetrad_len;

    // cursor から前方に走査し、長さ tetrad_len の G-run を探す
    for y in min_len..=max_len {
        let tetrad_start = cursor + y;
        if is_valid_tetrad_at(tetrad_start, cand.tetrad_len, data) {
            lengths.push(y);
        }
    }
    lengths
}
```

**制約条件**:
1. **最小長**: `min_acceptable_loop_length()` は通常 0 を返し、すでに 0 長 loop が存在する場合は 1 を返して、0 長 loop が連続しないようにする
2. **最大長**: `partial_length() + y + tetrad_len <= max_length` を保証し、候補が全体長制限を超えないようにする
3. **tetrad 一致**: `is_valid_tetrad_at()` は `cursor + y` の位置に `tetrad_len` 個連続する G があるかを確認する

**展開例**:
- 候補状態: `start=212, tetrad_len=3, y1=5`（1 つ目の loop は割り当て済み）
- `cursor = 212 + 3 + 5 = 220`（2 つ目の tetrad はここから始まる）
- 220 と 222 の位置にそれぞれ長さ 3 の G-run が見つかったとする
- 次の 2 つの新しい候補が生成される:
  - `{start=212, tetrad=3, y1=5, y2=0}`（`loop2` の長さが 0）
  - `{start=212, tetrad=3, y1=5, y2=2}`（`loop2` の長さが 2）

### 1.4 スコア式（`score`）

**Legacy 公式**（C++ 版から継承。`qgrs.h` を参照）:
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

**パラメータの意味**:
- `gmax`: 残りの長さ余裕（`max_length` から 4 つの tetrad と 1 を引いたもの）
- `gavg`: 3 つの loop 長の差の平均で、均一性を表す
- `bonus`: tetrad 数に応じた加点（2-tetrad では加点なし、3-tetrad では `1 x gmax` の加点）

**スコアの性質**:
1. **長さペナルティ**: 全長が `max_length` に近づくほど `gmax` は小さくなり、スコアも下がる
2. **均一性報酬**: loop 長が近いほど（`gavg` が小さいほど）スコアは高くなる
3. **tetrad 加点**: tetrad が多い候補（3 以上）は追加ボーナスを得る

#### 1.4.1 `--max-g4-length` が tetrads ごとに与える実際の影響

QGRS-Rust において `--max-g4-length` は単なる出力フィルタではない。候補生成、loop 展開、可否判定、スコア計算のすべてに影響する。各候補内部で実際に使われる `max_length` は次の通りである:

```text
max_length(n) = min(legacy_cap(n), --max-g4-length)
legacy_cap(n) = 30, n < 3 のとき
legacy_cap(n) = 45, n >= 3 のとき
```

ここで `n` は候補の `tetrads` 値を表す。したがって、その影響はまず次のように要約できる:

| tetrads | 候補が実際に使う `max_length` | `--max-g4-length` の効果 |
| --- | --- | --- |
| `n < 3` | `min(30, L)` | `L < 30` のときにのみ、2-tetrad 候補の許容全長、スコア上限、loop 展開空間が実際に縮小する |
| `n >= 3` | `min(45, L)` | `L < 45` のときにのみ、3-tetrad 以上の候補の許容全長、スコア上限、loop 展開空間が実際に縮小する |

ここで `L = --max-g4-length` である。既定値 `L=45` に注目すると、旧来の QGRS の段階的ルールは次の `tetrads -> 許容全長 -> スコアへの影響` の表に書き換えられる:

| tetrads | 許容全長上限 | スコアへの影響 |
| --- | --- | --- |
| 2 | `<= 30` | `gscore = floor(21 - gavg)` |
| 3 | `<= 45` | `gscore = floor(64 - gavg)` |
| 4 | `<= 45` | `gscore = floor(84 - gavg)` |
| 5 | `<= 45` | `gscore = floor(96 - gavg)` |
| 6 | `<= 45` | `gscore = floor(100 - gavg)` |
| 7 | `<= 45` | `gscore = floor(96 - gavg)` |
| 8 | `<= 45` | `gscore = floor(84 - gavg)` |
| 9 | `<= 45` | `gscore = floor(64 - gavg)` |
| 10 | `<= 45` | `gscore = floor(36 - gavg)` |

補足:
- ここでの `gavg = (|y1-y2| + |y2-y3| + |y1-y3|) / 3` は、3 つの loop 長の不均衡に対するペナルティ項である。
- 2-tetrad 候補は既定で `30 bp` の長さ予算しか持たないため、`--max-g4-length` を `45` より大きくしてもスコアや長さ判定はそれ以上緩和されない。
- 3-tetrad 以上の候補は既定で `45 bp` の長さ予算を共有する。`--max-g4-length < 45` の場合は、スコア上限と loop 展開空間が同時に小さくなる。
- `--max-g4-length` は seed 段階で列挙できる最大 tetrads も制限する: `max_tetrads_allowed = min(max_g_run, floor(L/4))`。そのため `L` が小さいと、高 tetrad 候補の一部は BFS 展開前に生成されなくなる。

**典型的なスコア比較**（`max_length=45`, `min_tetrads=2`）:
- `{tetrad=3, y1=5, y2=5, y3=5}`: `gmax=32`, `gavg=0`, `bonus=32` -> `score=64`
- `{tetrad=3, y1=5, y2=2, y3=5}`: `gmax=32`, `gavg=2`, `bonus=32` -> `score=62`
- `{tetrad=3, y1=5, y2=1, y3=6}`: `gmax=32`, `gavg=2.67`, `bonus=32` -> `score=61`

### 1.5 C++ 版との比較

**C++ 実装**（`qgrs.cpp`）:
```cpp
void findCandidates(const char* seq, int len) {
    queue<Candidate> q;
    vector<G4> results;
    // シード生成: 単純な文字列走査（SIMD なし）
    for (int i = 0; i < len; ++i) {
        if (seq[i] == 'G' || seq[i] == 'g') {
            int run_len = 1;
            while (i+run_len < len && (seq[i+run_len]=='G' || seq[i+run_len]=='g'))
                ++run_len;
            for (int off = 0; off <= run_len - min_tetrads; ++off)
                q.push(Candidate(i+off, min_tetrads+off));
        }
    }
    // BFS ループ: ロジック自体は Rust 版と完全に同じ
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

**主な差分**:
| 観点 | Rust | C++ |
|------|------|-----|
| G-run 走査 | `memchr2` SIMD（約 10 倍高速） | 1 バイトずつの `while` ループ |
| キュー型 | `VecDeque<G4Candidate>` | `std::queue<Candidate>` |
| メモリ管理 | `Arc<Vec<u8>>` によるゼロコピー | 毎回 `std::string` をコピー |
| 並列化 | Rayon による自動分割 | 単一スレッド逐次実行 |
| スコア式 | **完全に同一**（逐語的に移植） | 元の legacy 公式 |

**性能測定結果**（`big.txt`、10MB 配列）:
- C++: 約 45 秒（単一スレッド）
- Rust (`stream`): 約 8 秒（8 コア）
- Rust (`mmap`): 約 6 秒（8 コア + ゼロコピー）

### 1.6 典型例の解析

**シナリオ**: 位置 212 で 5 個連続の G を検出
- 配列: `...GGGGGACGTGGGACGTGGG...`
- パラメータ: `min_tetrads=2, max_g_run=5, max_length=45`

**BFS 展開の流れ**:
1. **シード**: 4 個の候補を生成（`offset 0-3` が `tetrad=2,3,4,5` に対応）
2. **第 1 展開**（`y1` を埋める）:
   - 候補 `{212, tetrad=3}` が位置 220 に G-run を見つける
   - 実行可能な `y1` の値: `[5, 6, 7]`（220, 221, 222 の位置に対応する tetrad）
3. **第 2 展開**（`y2` を埋める）:
   - 候補 `{212, tetrad=3, y1=5}` は次のように展開される:
     - `{..., y1=5, y2=2}`（`gscore=19`）
     - `{..., y1=5, y2=5}`（`gscore=17`）
4. **第 3 展開**（`y3` を埋める）:
   - 候補 `{..., y1=5, y2=2}` は最終的に:
     - `{..., y1=5, y2=2, y3=5}`（`gscore=19`, `length=27`）
5. **収集**: `viable()` の条件（`length <= 45`, `score >= 0`）を満たすため、`raw_g4s` に追加される

**重複排除の結果**:
- 生の出力には、異なる loop 構成を持つ `start=212` の候補が複数含まれる
- `consolidate_g4s` の第 1 段階では、HashMap Entry API により `gscore=19` の版が残る
- 最終出力: `start=212, length=27, tetrads=3, y1=5, y2=2, y3=5, gscore=19`

---
## 2. ファミリー統合（`consolidate_g4s`）詳細
### 第 1 段階 - 重複排除:
`HashMap<DedupKey, G4>` を使い、`(start, end, sequence)` で重複を取り除きます。

**`DedupKey` の定義**:
```rust
#[derive(Eq, PartialEq, Hash)]
struct DedupKey {
    start: usize,           // G4 の開始位置（1-based）
    end: usize,             // G4 の終了位置（exclusive）
    slice: SequenceSlice,   // 実際の配列バイトを含むスライス参照
}
```
`SequenceSlice` はバイト内容に基づくハッシュと等値比較を実装しているため、同じ座標かつ同じ配列を持つ候補は重複として認識されます。

**重複排除ロジック**:
```rust
let mut best_by_key: HashMap<DedupKey, G4> = HashMap::new();
for g in raw_g4s.into_iter() {
    let key = DedupKey::new(&g);
    match best_by_key.entry(key) {
        Entry::Vacant(slot) => slot.insert(g),
        Entry::Occupied(mut slot) => {
            if g.gscore > slot.get().gscore {
                slot.insert(g);  // より高スコアの版に置き換える
            }
        }
    }
}
```
同じ key を持つ候補が複数ある場合は、最も高い `gscore` を持つものだけを残します。

**これが必要な理由**: BFS は、`start=212, y1=5, y2=2, y3=5, gscore=19` と `y1=5, y2=1, y3=6, gscore=17` のように、座標は同じでも loop 構成が異なる候補を生成し得ます。HashMap により最良構成だけが残ります。HashMap の検索は `O(1)` で、C++ の `set` は `O(log n)` です。


### 第 2 段階 - グループ化:
既存のファミリーを線形走査し、`belongs_in` で新しい候補がファミリー内のいずれかのメンバーと重なるかを判定します。

**重なり判定（`overlapped`）**:
```rust
fn overlapped(a: &G4, b: &G4) -> bool {
    let a_start = a.start as isize;
    let a_end = (a.start + a.length) as isize;
    let b_start = b.start as isize;
    let b_end = (b.start + b.length) as isize;
    // 4 つの区間交差パターンのいずれかが真なら重なる
    (a_start >= b_start && a_start <= b_end) ||  // a の開始点が b の中にある
    (a_end >= b_start && a_end <= b_end) ||      // a の終了点が b の中にある
    (b_start >= a_start && b_start <= a_end) ||  // b の開始点が a の中にある
    (b_end >= a_start && b_end <= a_end)         // b の終了点が a の中にある
}
```
これは 2 つの G4 区間 `[start, start+length)` が少なくとも 1 塩基位置を共有するかどうかを判定します。

**貪欲な帰属戦略**:
```rust
for g4 in deduped.into_iter() {
    let mut inserted = false;
    for family in &mut families {
        if belongs_in(&g4, family) {  // ファミリー内のいずれかと重なる
            family.push(g4.clone());
            inserted = true;
            break;  // 最初に一致したファミリーで打ち切る
        }
    }
    if !inserted {
        families.push(vec![g4]);  // 新しいファミリーを作る
    }
}
```
- 候補は `(start, end)` でソートされたあと、既存ファミリーに順に照合される
- 最初に一致したファミリーが見つかった時点でそこに追加し、後続のファミリーは見ない
- どのファミリーとも重ならなければ、新しいファミリーを作る
- **推移的な連結**: A と B が重なり、B と C が重なるなら、A と C が直接重ならなくても 3 つとも同じファミリーに入る

**順序依存性**: 重複排除後の候補は `(start, end)` でソートしておく必要があります。そうしないと、同じ候補集合でも chunk / stream モードで反復順が異なり、異なるファミリー分けになって出力が一致しなくなります。

### 第 3 段階 - 最良候補の選択:
各ファミリーから `gscore` が最も高いメンバーを出力します。
生物学的な意味としては、重なった構造は同じ特徴を表している可能性があるため、最良の構成を残します。

### 2.1 メタデータ出力と `--overlap`

`consolidate_g4s` は現在 `(Vec<G4>, Vec<(usize, usize)>)` を返します。前者は各ファミリーの勝者、後者は各ファミリーの `[start, end]` 範囲です。CLI の `--overlap` フラグを有効にすると、主 CSV / Parquet 出力と一緒に 2 種類のデバッグ用ファイルが追加で書き出されます。

1. `{base}.overlap.csv`: `render_csv_results` を使ってソート済み raw hits をそのまま出力します。列は主 CSV と完全に同じなので diff が容易です。
2. `{base}.family.csv`: `render_family_ranges_csv` により `family_index,start,end` を出力し、どのファミリーが統合されたかを素早く確認できます。

ストリーミングパイプラインでは、`StreamChromosomeResults` を通じて `hits`、`family_ranges`、および任意の `raw_hits` をコールバックに渡します。

```rust
pub struct StreamChromosomeResults {
    pub hits: Vec<G4>,
    pub family_ranges: Vec<(usize, usize)>,
    pub raw_hits: Option<Vec<G4>>, // --overlap が有効なときのみ埋まる
}
```

`--overlap` のデフォルトは無効で、不要な `Vec` clone を避けます。ユーザーが有効化した場合:

- inline モードでは `--output` も必須で、そうでないとデバッグファイルのベースパスを推定できない
- FASTA / `mmap` モードでは各染色体ごとに 3 つのファイル（主出力 + `.overlap.csv` + `.family.csv`）を書き出し、すべて `sanitize_name` と重複サフィックス規則を引き継ぐ
- stream モードでは各染色体の処理が終わるたびに追加 CSV を即時 flush するため、入力サイズに応じてメモリ使用量が増えない

この仕組みにより、chunk vs. stream や Rust vs. C++ の差異を調べる際に、コアアルゴリズムを変更せず raw hits や family ranges を直接 diff できます。
