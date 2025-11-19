# QGRS-Rust

Rust port of the original `qgrs.cpp` logic. The crate now focuses solely on discovering G-quadruplexes (G4s) and exporting the hits as CSV or Parquet. The legacy `find-repeat-G` binary is kept only as a stub; run the `qgrs` target instead.

## Usage

```bash
# build once (optimized)
cargo build --release --bin qgrs

# G-quadruplex search mirroring qgrs.cpp (CSV output by default)
cargo run --release --bin qgrs -- --sequence GGGGAGGGGAGGGGAGGGG --min-tetrads 4 --min-score 17 --format csv
# or feed a FASTA/plain sequence, pick mmap/stream input, and emit one file per chromosome
cargo run --release --bin qgrs -- --file path/to/sequence.fa --mode stream --min-tetrads 4 --format parquet --output-dir ./qgrs_out
```

### Recommended workflow

1. 构建：`cargo build --release --bin qgrs`
2. 运行
   - 处理内联序列：`target/release/qgrs --sequence ...`
   - 处理 FASTA/文本文件：`target/release/qgrs --file input.fa --mode mmap|stream --output-dir out_dir ...`
3. 选择输出
   - CSV（内联序列默认 stdout，或使用 `--output` 指定文件；文件输入会在 `--output-dir` 下按染色体生成多个 CSV）
   - Parquet（内联序列需 `--output`，文件输入在 `--output-dir` 下生成多个 Parquet）

## Implementation highlights

- `src/qgrs/` owns everything: sequence normalization, mmap/stream loaders, the translated `find` logic, and CSV/Parquet exporters.
- `src/bin/qgrs.rs` is the only shipping CLI. It accepts inline sequences or FASTA/plain files, splits FASTA inputs by chromosome header, and provides knobs for minimum tetrads, g-score, output format, input mode, and the destination (`--output` or `--output-dir`).
- `memchr` accelerates the internal `GRunScanner`, which seeds every candidate G-run without copying the input slice.
- `memmap2` and `BufReader` cover the two input strategies so huge genomes can be processed without buffering entire files into RAM.
- `src/qgrs/stream.rs` implements the true streaming pipeline: FASTA bytes are parsed incrementally, candidates are expanded directly from a sliding window, and no chromosome-sized `String` is ever materialized when `--mode stream` is selected.

## Testing

```bash
cargo test
```

Use `cargo run --release` for production input sizes; this enables full compiler optimizations and targets the M1's vector units.
