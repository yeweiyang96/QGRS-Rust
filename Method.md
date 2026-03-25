## 1. BFS Candidate Expansion Algorithm

BFS (breadth-first search) candidate expansion is the core search logic of QGRS-Rust. It enumerates every valid G-quadruplex structure in a sequence. The algorithm has four stages: seed generation, the breadth-first loop, loop discovery, and score filtering.

### 1.1 Seed Generation (`seed_queue`)

**Goal**: scan the sequence to identify every potential G-run (a contiguous stretch of G bases), then generate initial candidates over the Cartesian product of "allowed tetrad counts x allowed offsets". Every later BFS expansion starts from these seeds.

**Implementation details**:
```rust
fn seed_queue(
    cands: &mut VecDeque<G4Candidate>,
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    limits: ScanLimits,
) {
    // Limit the tetrad count so we do not generate seeds beyond max_g_run
    // or max_g4_length / 4
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
                break; // the tetrad span alone already exceeds the overall length limit
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

**Key optimizations**:
- `GRunScanner::new(&seq.normalized, min_tetrads)` still relies on `memchr2(b'g', b'G')`, but now returns only runs with length >= `min_tetrads`, which reduces filtering work in upper layers. The SIMD scan is about 10x faster than checking one byte at a time.
- `max_tetrads_allowed = min(max_g_run, max_g4_length/4)` combines the G-run length bound with the global length bound, preventing redundant seeds that could never pass the later `max_length` check.
- When a run is longer than `max_g_run`, `run_len.min(max_tetrads_allowed)` truncates the available tetrad counts, but the later offset enumeration still covers the whole run. That means a long run that exceeds the "maximum contiguous G" limit still enters BFS through bounded windows instead of being dropped entirely.
- For each run, all offsets are enumerated again for every valid tetrad count. This ensures that long runs such as `GGGGG` cover every sub-interval without rescanning the sequence multiple times.

**Example** (`min_tetrads=2`, `max_g_run=5`, `max_g4_length=40`):
- Sequence fragment: `GGGGG...` (`run_len=5`)
- Valid tetrad counts: 2, 3, 4, 5 (constrained by `max_g4_length`, because `4x5=20 <= 40`)
- For `tetrad=3`, the allowed offsets are `0..=2`, so three start positions are generated: `start`, `start+1`, `start+2`
- In chunk mode, the window layer additionally enforces `start < primary_end` to avoid emitting duplicate raw hits from overlapping windows, but `seed_queue` itself stays identical to the stream and inline paths

### 1.2 Breadth-First Loop (`find_raw_with_sequence`)

**Core logic**: pop candidates from the queue in order. If a candidate is complete (all three loops are assigned), score it and keep it. Otherwise, expand it by filling in the next loop.

**Pseudo-code**:
```rust
let mut cands: VecDeque<G4Candidate> = seed_queue;
let mut raw_g4s: Vec<G4> = Vec::new();

while let Some(cand) = cands.pop_front() {
    if cand.complete() {            // check whether y1, y2, and y3 are all assigned
        if cand.viable(min_score) { // pass minimum score and maximum length checks
            raw_g4s.push(G4::from_candidate(&cand));
        }
    } else {
        // Expand the candidate: enumerate every possible length for the next unassigned loop
        for expanded in cand.expand(limits) {
            if expanded.partial_length() <= max_length {
                cands.push_back(expanded);
            }
        }
    }
}
```

**Actual code** (`mod.rs` lines 850-863):
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

**Completeness check**:
- `y1.is_none()` -> the first loop still needs to be assigned
- `y2.is_none()` -> the second loop still needs to be assigned
- `y3.is_some()` -> all loops are assigned, so the candidate is complete

**Viability check** (`viable`):
```rust
impl G4Candidate {
    fn viable(&self, limits: &ScanLimits) -> bool {
        let total_len = self.total_length();
        total_len <= limits.max_length && self.score(limits) >= limits.min_score
    }
}
```

**What BFS guarantees**:
- **Completeness**: every valid combination is enumerated, because the queue ensures that every candidate is expanded
- **Correctness**: pruning by `partial_length()` avoids invalid expansions; candidates that exceed `max_length` are never pushed back into the queue
- **Order independence**: candidates at the same coordinates may enter the queue in different orders (depending on seed order), but the final deduplication stage still guarantees uniqueness

### 1.3 Loop Discovery Mechanism (`expand`)

**Goal**: from the current candidate state (with some loops already assigned), scan forward to the next G-run and enumerate all valid lengths for the next loop.

**Implementation details**:
```rust
impl G4Candidate {
    fn expand(&self, data: &SequenceData, limits: &ScanLimits) -> Vec<G4Candidate> {
        let cursor = self.next_cursor_position();  // position after the current tetrad
        let loop_lengths = find_loop_lengths_from(cursor, data, limits, self);

        loop_lengths.into_iter()
            .map(|y| {
                let mut new_cand = self.clone();
                new_cand.assign_next_loop(y);  // fill y1 / y2 / y3
                new_cand
            })
            .collect()
    }
}
```

**Cursor calculation rules**:
- `y1` unassigned -> `cursor =` end position of tetrad 1 (`start + tetrad_len`)
- `y2` unassigned -> `cursor =` end of tetrad 1 + `y1` + tetrad 2
- `y3` unassigned -> `cursor =` end of tetrad 1 + `y1` + tetrad 2 + `y2` + tetrad 3

**Logic in `find_loop_lengths_from`**:
```rust
fn find_loop_lengths_from(cursor: usize, data: &SequenceData, limits: &ScanLimits, cand: &G4Candidate) -> Vec<usize> {
    let mut lengths = Vec::new();
    let min_len = cand.min_acceptable_loop_length();  // 0 or 1 (if a previous loop is already 0)
    let max_len = limits.max_length - cand.partial_length() - cand.tetrad_len;

    // Scan forward from cursor and look for a G-run of length tetrad_len
    for y in min_len..=max_len {
        let tetrad_start = cursor + y;
        if is_valid_tetrad_at(tetrad_start, cand.tetrad_len, data) {
            lengths.push(y);
        }
    }
    lengths
}
```

**Constraints**:
1. **Minimum length**: `min_acceptable_loop_length()` returns 0 by default, or 1 if an earlier loop is already 0, to prevent two consecutive zero-length loops
2. **Maximum length**: ensures `partial_length() + y + tetrad_len <= max_length`, so the candidate never exceeds the global length limit
3. **Tetrad match**: `is_valid_tetrad_at()` checks whether there are `tetrad_len` consecutive G's at `cursor + y`

**Expansion example**:
- Candidate state: `start=212, tetrad_len=3, y1=5` (the first loop is already assigned)
- `cursor = 212 + 3 + 5 = 220` (the second tetrad should start here)
- Assume length-3 G-runs are found at positions 220 and 222
- Two new candidates are generated:
  - `{start=212, tetrad=3, y1=5, y2=0}` (`loop2` length is 0)
  - `{start=212, tetrad=3, y1=5, y2=2}` (`loop2` length is 2)

### 1.4 Scoring Formula (`score`)

**Legacy formula** (inherited from the C++ version; see `qgrs.h`):
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

**Parameter meanings**:
- `gmax`: remaining span (`max_length` minus four tetrads plus 1)
- `gavg`: average pairwise loop-length difference, which measures uniformity
- `bonus`: tetrad-count bonus (no bonus for 2 tetrads, `1 x gmax` bonus for 3 tetrads)

**Scoring characteristics**:
1. **Length penalty**: the closer the total length gets to `max_length`, the smaller `gmax` becomes, and the lower the score gets
2. **Uniformity reward**: the more similar the loop lengths are (smaller `gavg`), the higher the score
3. **Tetrad bonus**: candidates with more tetrads (3+) receive an extra bonus

#### 1.4.1 Practical effect of `--max-g4-length` across tetrad counts

In QGRS-Rust, `--max-g4-length` is not just an output filter. It simultaneously affects candidate generation, loop expansion, viability checks, and scoring. The `max_length` actually used inside each candidate is:

```text
max_length(n) = min(legacy_cap(n), --max-g4-length)
legacy_cap(n) = 30, when n < 3
legacy_cap(n) = 45, when n >= 3
```

Here `n` is the candidate's `tetrads` value. Its effect can therefore be summarized first as:

| tetrads | effective candidate `max_length` | effect of `--max-g4-length` |
| --- | --- | --- |
| `n < 3` | `min(30, L)` | Only when `L < 30` does it actually reduce the allowed total length, score ceiling, and loop expansion space for 2-tetrad candidates |
| `n >= 3` | `min(45, L)` | Only when `L < 45` does it actually reduce the allowed total length, score ceiling, and loop expansion space for 3-tetrad and larger candidates |

Here `L = --max-g4-length`. If we focus on the default setting `L=45`, the legacy QGRS tiered rule can be rewritten as the following `tetrads -> allowed total length -> score impact` table:

| tetrads | allowed total length | score impact |
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

Notes:
- Here `gavg = (|y1-y2| + |y2-y3| + |y1-y3|) / 3`, the penalty term for loop-length imbalance.
- A 2-tetrad candidate gets only a `30 bp` length budget by default, so raising `--max-g4-length` above `45` does not further relax its scoring or length check.
- Candidates with 3 tetrads or more share a default `45 bp` length budget; when `--max-g4-length < 45`, both their score ceiling and loop expansion space shrink together.
- `--max-g4-length` also limits the maximum tetrad count that can be seeded: `max_tetrads_allowed = min(max_g_run, floor(L/4))`. When `L` is small, some high-tetrad candidates disappear before BFS expansion even starts.

**Typical score comparison** (`max_length=45`, `min_tetrads=2`):
- `{tetrad=3, y1=5, y2=5, y3=5}`: `gmax=32`, `gavg=0`, `bonus=32` -> `score=64`
- `{tetrad=3, y1=5, y2=2, y3=5}`: `gmax=32`, `gavg=2`, `bonus=32` -> `score=62`
- `{tetrad=3, y1=5, y2=1, y3=6}`: `gmax=32`, `gavg=2.67`, `bonus=32` -> `score=61`

### 1.5 Comparison with the C++ Version

**C++ implementation** (`qgrs.cpp`):
```cpp
void findCandidates(const char* seq, int len) {
    queue<Candidate> q;
    vector<G4> results;
    // Seed generation: simple string scan (no SIMD)
    for (int i = 0; i < len; ++i) {
        if (seq[i] == 'G' || seq[i] == 'g') {
            int run_len = 1;
            while (i+run_len < len && (seq[i+run_len]=='G' || seq[i+run_len]=='g'))
                ++run_len;
            for (int off = 0; off <= run_len - min_tetrads; ++off)
                q.push(Candidate(i+off, min_tetrads+off));
        }
    }
    // BFS loop: logically identical to the Rust version
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

**Key differences**:
| Dimension | Rust | C++ |
|------|------|-----|
| G-run scanning | `memchr2` SIMD (~10x faster) | byte-by-byte `while` loop |
| Queue type | `VecDeque<G4Candidate>` | `std::queue<Candidate>` |
| Memory management | `Arc<Vec<u8>>` zero-copy | `std::string` copy each time |
| Parallelism | Rayon automatic chunking | single-threaded sequential execution |
| Scoring formula | **exactly the same** (line-by-line translation) | original legacy formula |

**Performance results** (`big.txt`, 10 MB sequence):
- C++: about 45 seconds (single-threaded)
- Rust (`stream`): about 8 seconds (8 cores)
- Rust (`mmap`): about 6 seconds (8 cores + zero-copy)

### 1.6 Typical Example Walkthrough

**Scenario**: five consecutive G's are detected at position 212
- Sequence: `...GGGGGACGTGGGACGTGGG...`
- Parameters: `min_tetrads=2, max_g_run=5, max_length=45`

**BFS expansion process**:
1. **Seeds**: generate 4 candidates (`offset 0-3` corresponding to `tetrad=2,3,4,5`)
2. **First expansion round** (fill `y1`):
   - Candidate `{212, tetrad=3}` finds a G-run at position 220
   - Viable `y1` values discovered: `[5, 6, 7]` (corresponding to tetrads at positions 220, 221, 222)
3. **Second expansion round** (fill `y2`):
   - Candidate `{212, tetrad=3, y1=5}` expands into:
     - `{..., y1=5, y2=2}` (`gscore=19`)
     - `{..., y1=5, y2=5}` (`gscore=17`)
4. **Third expansion round** (fill `y3`):
   - Candidate `{..., y1=5, y2=2}` finally expands into:
     - `{..., y1=5, y2=2, y3=5}` (`gscore=19`, `length=27`)
5. **Collection**: it passes `viable()` (`length <= 45`, `score >= 0`) and is pushed into `raw_g4s`

**Deduplication result**:
- The raw output contains multiple candidates with `start=212` but different loop layouts
- In phase 1 of `consolidate_g4s`, deduplication keeps the `gscore=19` variant by using the HashMap Entry API
- Final output: `start=212, length=27, tetrads=3, y1=5, y2=2, y3=5, gscore=19`

---
## 2. Family Consolidation (`consolidate_g4s`) in Detail
### Phase 1 - Deduplication:
Use `HashMap<DedupKey, G4>` to deduplicate by `(start, end, sequence)`.

**Definition of `DedupKey`**:
```rust
#[derive(Eq, PartialEq, Hash)]
struct DedupKey {
    start: usize,           // G4 start position (1-based)
    end: usize,             // G4 end position (exclusive)
    slice: SequenceSlice,   // slice reference containing the actual sequence bytes
}
```
`SequenceSlice` implements hashing and equality by byte content, so candidates with the same coordinates and the same sequence are recognized as duplicates.

**Deduplication logic**:
```rust
let mut best_by_key: HashMap<DedupKey, G4> = HashMap::new();
for g in raw_g4s.into_iter() {
    let key = DedupKey::new(&g);
    match best_by_key.entry(key) {
        Entry::Vacant(slot) => slot.insert(g),
        Entry::Occupied(mut slot) => {
            if g.gscore > slot.get().gscore {
                slot.insert(g);  // replace with the higher-score version
            }
        }
    }
}
```
For multiple candidates with the same key, the one with the highest `gscore` is kept.

**Why this is needed**: BFS may generate candidates with identical coordinates but different loop layouts, such as `start=212, y1=5, y2=2, y3=5, gscore=19` versus `y1=5, y2=1, y3=6, gscore=17`. The HashMap guarantees that only the best configuration survives. HashMap lookup is `O(1)`, while the C++ `set` is `O(log n)`.


### Phase 2 - Grouping:
Linearly scan the existing families and use `belongs_in` to decide whether the new candidate overlaps with any member of a family.

**Overlap test (`overlapped`)**:
```rust
fn overlapped(a: &G4, b: &G4) -> bool {
    let a_start = a.start as isize;
    let a_end = (a.start + a.length) as isize;
    let b_start = b.start as isize;
    let b_end = (b.start + b.length) as isize;
    // Any of the four interval-intersection cases means overlap
    (a_start >= b_start && a_start <= b_end) ||  // a starts inside b
    (a_end >= b_start && a_end <= b_end) ||      // a ends inside b
    (b_start >= a_start && b_start <= a_end) ||  // b starts inside a
    (b_end >= a_start && b_end <= a_end)         // b ends inside a
}
```
This checks whether two G4 intervals, `[start, start+length)`, share any base position.

**Greedy assignment strategy**:
```rust
for g4 in deduped.into_iter() {
    let mut inserted = false;
    for family in &mut families {
        if belongs_in(&g4, family) {  // overlaps with any member in the family
            family.push(g4.clone());
            inserted = true;
            break;  // stop at the first matching family
        }
    }
    if !inserted {
        families.push(vec![g4]);  // create a new family
    }
}
```
- Candidates are sorted by `(start, end)` and then checked against existing families in order
- Once the first matching family is found, the candidate is inserted there and no later families are checked
- If the candidate overlaps with no family at all, it starts a new family
- **Transitive connectivity**: if A overlaps B and B overlaps C, all three end up in the same family even when A and C do not overlap directly

**Order dependence**: after deduplication, candidates must be sorted by `(start, end)`. Otherwise the same candidate set can be assigned to different families in chunk and stream modes because of different iteration orders, which would make the outputs diverge.

### Phase 3 - Pick the Winner:
For each family, output the member with the highest `gscore`.
Biological intuition: overlapping structures may represent the same feature, so the best-scoring configuration is kept.

### 2.1 Metadata Export and `--overlap`

`consolidate_g4s` now returns `(Vec<G4>, Vec<(usize, usize)>)`: the first item is the winner of each family, and the second item records the `[start, end]` range of every family. When the CLI `--overlap` flag is enabled, two debugging files are written alongside the main CSV or Parquet output:

1. `{base}.overlap.csv`: uses `render_csv_results` to dump sorted raw hits directly. The columns are exactly the same as the main CSV, which makes diffing easier.
2. `{base}.family.csv`: uses `render_family_ranges_csv` to output `family_index,start,end`, making it easy to locate which families were merged.

The streaming pipeline passes `hits`, `family_ranges`, and optional `raw_hits` to the callback through `StreamChromosomeResults`:

```rust
pub struct StreamChromosomeResults {
    pub hits: Vec<G4>,
    pub family_ranges: Vec<(usize, usize)>,
    pub raw_hits: Option<Vec<G4>>, // only populated when --overlap is enabled
}
```

`--overlap` is disabled by default to avoid unnecessary `Vec` clones. When the user enables it:

- inline mode must also provide `--output`, otherwise the base path for the debugging files cannot be inferred
- FASTA / `mmap` mode writes three files per chromosome (main output + `.overlap.csv` + `.family.csv`), and every file still uses `sanitize_name` and duplicate suffix handling
- stream mode flushes the additional CSV files as soon as each chromosome is finished, so memory usage does not grow with input size

This mechanism lets you debug chunk vs. stream or Rust vs. C++ differences by diffing raw hits or family ranges directly, without modifying the core algorithm.
