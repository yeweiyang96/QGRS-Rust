# QGRS-Rust <img alt="DNA icon" align="right" width="90" src="https://img.shields.io/badge/G4-scan-blueviolet?logo=databricks&logoColor=white" />

<div align="center">
 
![Rust 1.75+](https://img.shields.io/badge/Rust-1.75%2B-orange?logo=rust&logoColor=white)
![Platforms macOS/Linux](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-lightgrey?logo=apple&logoColor=white)
![Outputs CSV / Parquet](https://img.shields.io/badge/output-CSV%20%7C%20Parquet-0e83cd)
![Mode mmap / stream](https://img.shields.io/badge/mode-mmap%20%7C%20stream-6f42c1)

</div>

QGRS-Rust is a ground-up Rust rewrite of the [freezer333/qgrs-cpp](https://github.com/freezer333/qgrs-cpp) project. The refactor preserves the legacy scoring math while unlocking modern performance tricks: zero-copy SIMD-friendly scanners, Rayon-powered chromosome fan-out, and dual ingestion paths (`mmap` or streaming) that keep throughput high and memory usage predictable. New functionality includes Parquet export, strict CLI validation, and a true streaming FASTA parser capable of processing multi-hundred-gigabyte `.fa` archives without ever materializing an entire chromosome in RAM. The single `qgrs` CLI now covers inline sequences, mmap batches, and the ultra-large streaming mode in one coherent interface.

## Table of contents

- [✨ Features](#-features)
- [🧰 Requirements](#-requirements)
- [⚙️ Build](#️-build)
- [🧪 Usage](#-usage)
  - [Quick recipes](#quick-recipes)
  - [CLI reference](#cli-reference)
  - [Output schema](#output-schema)
- [✅ Testing & QA](#-testing--qa)
- [📊 Benchmarking tips](#-benchmarking-tips)

## ✨ Features

- Rust scanner mirrors the legacy scoring heuristics while benefiting from zero-copy iterators and Rayon parallelism.
- Memory-mapped (`mmap`) and streaming (`stream`) readers let you pick the best strategy per dataset.
- CSV/Parquet exporters always report 1-based, inclusive coordinates for genome-browser compatibility.
- CLI validation enforces sane tetrad, loop, and window settings to avoid silent misconfiguration.

## 🧰 Requirements

- Rust toolchain 1.75 or newer (install via [`rustup`](https://rustup.rs)).
- `cargo` is required for building, running, and testing.
- macOS and Linux are tested; Windows should work via WSL2.

## 🧱 架构

QGRS-Rust 的核心位于 `src/qgrs/`，每个模块对应搜索管线中的独立阶段：

- `data.rs`：定义 `ChromSequence`、`SequenceData`、`ScanLimits` 等零拷贝数据容器。
- `search.rs`：实现 G-run 扫描、BFS 候选扩展、计分与原始 `G4` 构造。
- `chunks.rs`：根据 `ScanLimits` 计算窗口与重叠，调度 `find_raw_*` 并合并 Rayon 结果。
- `consolidation.rs`：对原始命中去重、聚类，保留每个重叠家族的最高 `gscore`。
- `stream.rs`：实现 `StreamChromosome`/`StreamChunkScheduler`，增量解析超大 FASTA。
- `loaders.rs`：封装 mmap 与普通文件加载器，供 CLI 在批量模式下复用。
- `export.rs`：提供 CSV/Parquet 渲染器及错误类型，统一 1-based 坐标输出。
- `tests/`：集中单元与集成测试，确保 chunk/stream 模式结果一致。

`src/lib.rs` 只负责 re-export 公共 API，`src/bin/qgrs.rs` 将 CLI 选项映射到上述模块，保持入口简洁。

## ⚙️ Build

```bash
# clone and enter the workspace
git clone https://github.com/<your-org>/QGRS-Rust.git
cd QGRS-Rust

# dev build (fast iteration)
cargo build --bin qgrs

# optimized binary for large genomes
cargo build --release --bin qgrs

# the optimized binary lives here after a release build
target/release/qgrs --help
```

> Tip: use `cargo install --path . --bin qgrs` if you want the binary on your `~/.cargo/bin` for reuse across projects.

## 🧪 Usage

`qgrs` accepts either an inline sequence (`--sequence`) or an input file (`--file`). FASTA inputs are split per chromosome header, and each slice is processed independently. If you provide a file, choose either the memory-mapped (`mmap`) or buffered streaming (`stream`) pipeline with `--mode`. All examples below assume you already built the release binary (`target/release/qgrs`) or installed it as `qgrs`; use `cargo run --release --bin qgrs -- …` only when iterating locally. The banner below comes straight from `src/bin/qgrs.rs` so it always matches the binary.

```
Usage: qgrs -- [--sequence <SEQ> | --file <PATH>] [options]
Options:
   --sequence <SEQ>       Inline DNA/RNA sequence to scan
   --file <PATH>          Read sequences from FASTA (chromosomes split independently)
   --min-tetrads <N>      Minimum tetrads to seed (default 2)
   --min-score <S>        Minimum g-score (default 17)
   --max-g-run <N>        Maximum allowed G-run length (default 10)
   --max-g4-length <N>    Maximum allowed G4 length in bp (default 45)
   --format <csv|parquet> Output format (default csv)
   --output <PATH>        Destination file when using --sequence (required for parquet)
   --output-dir <DIR>     Directory for per-chromosome exports when using --file
   --mode <mmap|stream>   Input mode when using --file (default mmap)
   --help                 Show this message
```

### Quick recipes

```bash
# 1. Quick sanity check against a short inline sequence
target/release/qgrs \
   --sequence GGGGAGGGGAGGGGAGGGG \
   --min-tetrads 4 \
   --min-score 17 \
   --format csv

# 2. Process a FASTA file using mmap and emit CSV files per chromosome
target/release/qgrs \
   --file data/genome.fa \
   --mode mmap \
   --min-tetrads 4 \
   --min-score 20 \
   --output-dir ./qgrs_csv

# 3. Stream extremely large FASTA with Parquet output
target/release/qgrs \
   --file hg38.fa \
   --mode stream \
   --min-tetrads 3 \
   --max-g-run 12 \
   --max-g4-length 60 \
   --format parquet \
   --output-dir ./qgrs_parquet

# 4. Clamp loop and run lengths for custom heuristics while writing to stdout
target/release/qgrs \
   --sequence GGGGTTTTGGGGTTTTGGGGTTTTGGGG \
   --max-g-run 8 \
   --max-g4-length 45 \
   --output -
```

### CLI reference

| Flag                      | Description                                                                                | Default                  |
| ------------------------- | ------------------------------------------------------------------------------------------ | ------------------------ |
| `--sequence <SEQ>`        | Inline DNA sequence to scan (mutually exclusive with `--file`).                            | _none_                   |
| `--file <PATH>`           | FASTA or plain-text file containing one or more sequences.                                 | _none_                   |
| `--mode <mmap\|stream>`   | File ingestion strategy; `mmap` favors fast disks, `stream` lowers RAM.                    | `mmap`                   |
| `--min-tetrads <INT>`     | Minimum number of stacked tetrads required for a hit.                                      | `2`                      |
| `--min-score <INT>`       | Minimum legacy G-score threshold.                                                          | `17`                     |
| `--max-g-run <INT>`       | Upper bound for contiguous G-run length (must be ≥ `min-tetrads`).                         | `10`                     |
| `--max-g4-length <INT>`   | Upper bound for the full quadruplex length (must be ≥ `4 * min_tetrads`).                  | `45`                     |
| `--format <csv\|parquet>` | Output encoding. CSV defaults to stdout for inline sequences; Parquet requires a file/dir. | `csv`                    |
| `--output <FILE\|- >`     | Single output file (or `-` for stdout) when scanning inline sequences.                     | stdout for CSV           |
| `--output-dir <DIR>`      | Directory for per-chromosome files when reading FASTA/plain inputs.                        | _required with `--file`_ |

The CLI aborts with a descriptive error if incompatible parameters are provided (e.g., `--mode stream` without `--file`, or `--max-g-run < min-tetrads`). When scanning files you must pass `--output-dir`; when scanning inline sequences `--output` is optional for CSV but required for Parquet.

### Output schema

Both exporters emit the same fields (see `render_csv` and `write_parquet_from_results` in `src/qgrs/mod.rs`):

| Column           | Meaning                                                                                 |
| ---------------- | --------------------------------------------------------------------------------------- |
| `start`          | 1-based inclusive start coordinate of the hit within the processed sequence/chromosome. |
| `end`            | 1-based inclusive end coordinate (`start + length - 1`).                                |
| `length`         | Total number of bases spanned by the quadruplex.                                        |
| `tetrads`        | Count of stacked tetrads contributing to the hit.                                       |
| `y1`, `y2`, `y3` | Loop lengths between successive G-runs (0 means no spacer).                             |
| `gscore`         | Legacy G-score used for filtering and ranking candidates.                               |
| `sequence`       | Exact G4 motif sequence extracted from the input.                                       |

CSV output always includes the header `start,end,length,tetrads,y1,y2,y3,gscore,sequence`. When scanning FASTA inputs, each chromosome is written to its own file (so the filename, not a column, captures the chromosome name). Parquet exports contain the same columns using Arrow types (`UInt64` for coordinates/lengths, `Int32` for loop lengths and score, and UTF-8 for sequences).

## Testing & QA

```bash
# unit tests
cargo test

# lint + formatting (optional but recommended before sending patches)
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
```

## Benchmarking tips

```bash
time target/release/qgrs --file aaa.fa --mode mmap   --max-g4-length 32 --max-g-run 8  --output-dir out-mmap
time target/release/qgrs --file aaa.fa --mode stream --max-g4-length 45 --max-g-run 10 --output-dir out-stream
```

Track `real` time, CPU%, and RSS with your preferred profiler to decide whether `mmap` or `stream` is better for your environment. Always benchmark with `--release` builds to enable full optimizations.

### compare_modes consistency tester

`target/release/compare_modes` (defined in `src/bin/compare_modes.rs`) benchmarks and cross-checks the two ingestion pipelines against the same FASTA input. It scans every chromosome once with the mmap batch loader and once with the streaming reader, reports per-mode timings and hit counts, then diff-checks every field (`start`, `end`, `length`, loops, tetrads, gscore, sequence) to ensure both paths stay bit-for-bit aligned. The process exits with code `0` on success and `1` with detailed mismatch logs when discrepancies are detected.

**Before you run it**

- Build the binary at least once so that `target/release/compare_modes` exists (typically from a prior release build).
- Provide a readable FASTA file containing one or more chromosomes; the tool operates read-only and prints summaries to stdout.

**CLI arguments**

| Positional argument | Description                         | Default |
| ------------------- | ----------------------------------- | ------- |
| `<fasta>`           | Path to the FASTA file to compare.  | _none_  |
| `[min_tetrads]`     | Minimum tetrads threshold per scan. | `2`     |
| `[min_gscore]`      | Minimum g-score threshold per scan. | `17`    |

**Example commands** (using the compiled release binary; adjust the path if you install it elsewhere):

```bash
# Compare dme.fa with default thresholds
target/release/compare_modes dme.fa

# Tighten heuristics to 3 tetrads / gscore 20
target/release/compare_modes dme.fa 3 20

# Point at a chromosome subset file
target/release/compare_modes output/chromosome-2L.fa 4 30
```

During a run you will see individual sections for the mmap phase, the stream phase, a speed comparison, and the final consistency verdict. An error summary (up to 10 detailed mismatches) is printed before the program returns a non-zero exit status, which makes the tool suitable for automated regression checks.

