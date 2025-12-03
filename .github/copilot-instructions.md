# QGRS-Rust AI Coding Agent Instructions

## Project Overview
G-quadruplex (G4) detection tool rewritten in Rust. Scans DNA/RNA sequences using legacy scoring heuristics with modern zero-copy architecture. Dual pipeline: memory-mapped (`mmap`) for random access, streaming for multi-GB FASTA files.

## Architecture Fundamentals

### Core Search Pipeline (`src/qgrs/mod.rs`)
**Data Flow**: Raw bytes → `SequenceData` (Arc-wrapped lowercase) → BFS candidate expansion → raw `G4` hits → deduplication → family consolidation → final results

- **Zero-copy design**: `ChromSequence` stores `Arc<Vec<u8>>` (lowercase). Output sequences lazily uppercased via `G4::sequence()` using `OnceLock<String>`.
- **Chunking strategy**: Long sequences split at `chunk_size_for_limits()` with `compute_chunk_overlap()` margin. Each chunk processed via `find_raw_on_window_bytes()`, results merged globally.
- **Deduplication**: `consolidate_g4s()` uses `HashMap<DedupKey, G4>` keeping highest `gscore` per `(start, end, sequence)` tuple before family grouping. Critical: Must sort by `(start, end)` after dedup to ensure deterministic family formation across chunk/stream modes.

### BFS Candidate Expansion Algorithm
The search uses breadth-first expansion to enumerate all valid G4 structures:

1. **Seeding** (`seed_queue`): Scan sequence with `GRunScanner` (SIMD-accelerated `memchr2`) to find G-runs meeting `min_tetrads` ≤ length ≤ `max_g_run`. For each run, generate candidates at every valid offset (e.g., run of 5 Gs with `min_tetrads=2` yields candidates at positions 0-3 within the run for 2-tetrad, 3-tetrad structures).

2. **BFS Loop** (`find_raw_with_sequence`):
   ```rust
   while let Some(cand) = cands.pop_front() {
       if cand.complete() {           // All three loops (y1,y2,y3) assigned
           if cand.viable(min_score) { // Passes score & length checks
               raw_g4s.push(G4::from_candidate(&cand));
           }
       } else {
           for expanded in cand.expand() {  // Fill next unassigned loop
               cands.push_back(expanded);
           }
       }
   }
   ```

3. **Loop Discovery** (`G4Candidate::expand`): 
   - Determine next cursor position: after tetrad1 (if y1 unset), tetrad2 (if y2 unset), or tetrad3 (if y3 unset).
   - Call `find_loop_lengths_from(cursor)` to scan forward for matching tetrad runs, collecting all valid loop lengths `y` where:
     - `y >= min_acceptable_loop_length()` (0 unless another loop is already 0, then 1)
     - `cursor + y + tetrad_len - 1 < max_length` (candidate stays within length limit)
   - Clone candidate for each valid `y`, assign to next loop slot, push to queue if `partial_length() <= max_length`.

4. **Scoring** (`G4Candidate::score`): Legacy formula `floor(gmax - gavg + gmax*(tetrads-2))` where:
   - `gmax = max_length - (tetrads*4 + 1)`
   - `gavg = mean(|y1-y2|, |y2-y3|, |y1-y3|)`
   - Higher score = more uniform loop distribution within tighter overall structure.

### Family Consolidation (`consolidate_g4s`)
Merges overlapping G4 candidates into representative families:

1. **Deduplication Phase**:
   ```rust
   // Input: raw_g4s sorted by start position
   HashMap<DedupKey, G4>  // Key = (start, end, sequence_bytes)
   ```
   - For each candidate, check if key exists: insert if new, replace if `gscore` higher.
   - Collect HashMap values, sort by `(start, end)` to stabilize iteration order.
   - **Why**: BFS can generate multiple candidates with same coordinates but different loop configs (e.g., `(5,2,5)` vs `(5,1,6)` at position 212). We keep the best scorer to avoid family merge seeing duplicates.

2. **Family Grouping**:
   ```rust
   for g4 in deduped {
       for family in &mut families {
           if belongs_in(&g4, family) {  // Overlaps any member?
               family.push(g4.clone());
               inserted = true;
               break;
           }
       }
       if !inserted { families.push(vec![g4]); }
   }
   ```
   - `overlapped(a, b)`: Returns true if intervals `[a.start, a.start+a.length]` and `[b.start, b.start+b.length]` share any base.
   - Linear scan through families; first match wins (order-dependent, hence dedup sort requirement).
   - Each family = set of G4s sharing overlapping genomic coordinates.

3. **Winner Selection**:
   ```rust
   for family in families {
       let best = family.iter().max_by_key(|g| g.gscore)?;
       results.push(best.clone());
   }
   ```
   - Emit single highest-scoring member per family.
   - **Rationale**: Overlapping structures likely represent same biological feature; choose best-scoring configuration.

### Key Functions
```rust
// Public entry points
find_owned_bytes(Arc<Vec<u8>>, ...) -> Vec<G4>              // Default limits helper
find_owned_bytes_with_limits(Arc<Vec<u8>>, ...) -> Vec<G4>  // Zero-copy, parallel chunks

// Raw (unconsolidated) helpers
find_raw_bytes_no_chunking(Vec<u8>, ...) -> Vec<G4>         // Stream workers, window slices
find_raw_on_window_bytes(Arc<Vec<u8>>, ...) -> Vec<G4>      // Chunked mmap windows

// Internal building blocks
consolidate_g4s(Vec<G4>) -> Vec<G4>              // Dedup + family merge
maximum_length(tetrads, limits) -> usize         // 2 tetrads→30bp, 3+→45bp
```

### Streaming Mode (`src/qgrs/stream.rs`)
`StreamChromosome` maintains sliding window state, dispatches chunks to Rayon workers, merges results on chromosome boundary. No consolidation per-chunk—raw hits accumulated, then single `consolidate_g4s()` call on finish.

### CLI Orchestration (`src/bin/qgrs.rs`)
- **Inline sequence**: `--sequence` → normalize to lowercase bytes → `find_owned_bytes_with_limits()`
- **FASTA file**: 
  - `mmap` mode: `load_sequences_from_path()` → `Vec<ChromSequence>` → parallel `into_par_iter()` → per-chromosome `find_owned_bytes_with_limits()` → write outputs
  - `stream` mode: `process_fasta_stream_with_limits()` → callback per chromosome → write outputs
- **Output formats**: CSV (`render_csv_results`), Parquet (`write_parquet_results`). Both consume `Vec<G4>`.

## Critical Patterns

### Testing Parity
Test `big_sequence_internal_equals_chunked` ensures chunk/non-chunk paths produce identical results. When modifying `consolidate_g4s()` or overlap logic, this MUST pass.

### Coordinate System
- **Internal (0-based)**: All search logic uses 0-based half-open `[start, end)`.
- **Output (1-based)**: `G4::from_candidate()` adds 1 to start; `end` remains exclusive. CSV/Parquet show 1-based inclusive for genome browsers.

### Performance Hotspots
- `GRunScanner`: Uses `memchr2(b'g', b'G')` for SIMD-friendly seed detection.
- Rayon parallelism: Chunks processed via `into_par_iter().flat_map_iter()`. Thread count locked to `num_cpus::get()` in `main()`.
- Avoid per-candidate `String` allocation: `G4.sequence_view` holds slice reference, uppercase conversion deferred to `sequence()` call.

## Common Workflows

### Build & Test
```bash
cargo build --release --bin qgrs        # Optimized binary
cargo test --lib                        # Unit tests (includes chunk parity)
cargo test --lib -- --nocapture         # See println! output
```

### Adding New Output Format
1. Extend `OutputFormat` enum in `src/bin/qgrs.rs`
2. Add renderer in `src/qgrs/mod.rs` (e.g., `render_wig_density(&[G4]) -> String`)
3. Update both `process_inline_sequence()` and FASTA branches (mmap + stream callback)
4. Add `extension()` case and CLI `--format` parser

### Debugging Consolidation Issues
Use `find_raw_bytes_no_chunking()` (for ad-hoc slices) or `find_raw_on_window_bytes()` to collect unconsolidated hits, then print before/after `consolidate_g4s()`. Check `DedupKey` hash collisions and family overlap logic (`overlapped()`, `belongs_in()`).

## Project-Specific Conventions

- **Chinese comments acceptable**: Legacy collaboration used mixed English/Chinese; maintain readability for both.
- **Lowercase storage, uppercase output**: Never store uppercase sequences in `SequenceData` or `SequenceSlice`. Always normalize to lowercase on ingestion.
- **Legacy scoring preserved**: `G4Candidate::score()` math must match C++ original (`qgrs.h` reference). Formula: `floor(gmax - gavg + gmax*(tetrads-2))`.
- **No `pub` leakage**: Internal helpers (`SequenceData`, `GRunScanner`) remain crate-private. Only `find_*`, `render_*`, `write_*`, `ChromSequence`, `ScanLimits` exposed.

## File Organization
```
src/
├── lib.rs              # Re-exports qgrs module
├── qgrs/
│   ├── mod.rs          # Core search, consolidation, exporters
│   └── stream.rs       # Streaming FASTA parser + chunk scheduler
└── bin/
    ├── qgrs.rs         # Main CLI
    ├── compare_modes.rs       # Dev tool: diff mmap vs stream
    └── compare_csv_outputs.rs # Dev tool: validate output parity
```

## Avoiding Common Pitfalls

1. **Dedup order matters**: After HashMap-based dedup, MUST `sort_by((start, end))` before family grouping. Otherwise chunk/stream modes diverge due to insertion order variance.
2. **Don't re-chunk in stream workers**: `find_raw_bytes_no_chunking()` consumes scheduler windows directly to prevent nested chunking.
3. **Test with `big.txt`**: Reference dataset includes edge cases (overlapping families, 2-tetrad max_length=30). Always validate against legacy `big.csv`.
4. **Arc cloning is cheap**: `Arc<Vec<u8>>` clones only increment refcount. Safe to pass into Rayon closures without copying sequences.

## Integration Points

- **FASTA loaders**: `load_sequences_mmap()` / `load_sequences_stream()` in `mod.rs`, called by CLI.
- **Parquet schema**: 9 columns (start, end, length, tetrads, y1-y3, gscore, sequence). Arrow types in `write_parquet_from_results()`.
- **Stream callback**: `FnMut(String, Vec<G4>) -> io::Result<()>` receives chromosome name + consolidated results.

## When Modifying Core Logic

- Update `maximum_length()` only if validating against updated C++ reference.
- Changes to `consolidate_g4s()` require re-running full test suite + manual `big.txt` diff.
- New CLI flags: add to `usage()`, `run_env()` parser, and corresponding struct (`OutputFormat`, etc.).
- Rayon thread count intentionally fixed; avoid `RAYON_NUM_THREADS` env var to ensure reproducibility.
