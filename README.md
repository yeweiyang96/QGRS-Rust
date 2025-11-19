# find-repeat-G

High-throughput scanner for locating every stretch of consecutive `G` characters that meets a configurable minimum length ("tetrads"). Designed around Apple Silicon (ARM64) memory bandwidth constraints so it can sweep 10&nbsp;MB sequences without extra copies.

## Usage

```bash
# memory-mapped scan (fastest when the sequence fits comfortably in RAM)
cargo run --release -- <min_tetrads> <path-to-sequence> --mmap

# streaming scan (bounded RAM, parses FASTA headers line-by-line)
cargo run --release -- <min_tetrads> <path-to-sequence> --stream
```

Each qualifying substring is written to standard output as a tab-delimited line:

```
<start-index>\t<length>
```

`start-index` is zero-based. For a run of length `L`, the tool emits every possible substring whose length ranges from `min_tetrads` through `L`.

## Implementation highlights

- **Project layout**: `src/config.rs` encapsulates CLI parsing, `src/scanner/` hosts the mmap and streaming backends plus shared iterators, and `src/main.rs` wires everything together.
- `memmap2` exposes the input file through the OS page cache, avoiding the cost of copying ~10&nbsp;MB buffers.
- The streaming mode uses a `BufReader` that strips FASTA headers/whitespace and keeps only the current `G`-run state, so arbitrarily large inputs stay within a tiny memory budget.
- `memchr` provides a NEON-accelerated search primitive on Apple M1, skipping directly to the next `G` byte.
- A custom `GRunScanner` iterator walks the mapped slice once, so even worst-case inputs stay cache-friendly.
- Output streaming uses a 1&nbsp;MiB `BufWriter` to reduce syscalls when enumerating large numbers of substrings.

## Testing

```bash
cargo test
```

Use `cargo run --release` for production input sizes; this enables full compiler optimizations and targets the M1's vector units.
