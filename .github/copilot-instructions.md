# QGRS-Rust AI Coding Agent Instructions

## 项目概览
QGRS-Rust 是 freezer333/qgrs-cpp 的 Rust 版本，通过原始的 G-score 公式扫描 DNA/RNA 序列。特点：

- `mmap` 与 `stream` 两种读入模式，覆盖常规和超大 FASTA；
- 完整零拷贝：内部统一小写 `Arc<Vec<u8>>`，输出时再 uppercase；
- BFS 扩展 + 去重/家族合并，确保 chunk 与 stream 完全一致；
- CLI 支持 inline 序列、批量 FASTA、CSV/Parquet 导出，并附带差异/基准工具。

## 模块与 API 边界
- `src/lib.rs` 只暴露 `pub mod qgrs;`，crate 根不再 re-export 任何函数。
- 对外 API 通过 `qgrs_rust::qgrs::*`，当前仅为 CLI/内部工具服务，默认不承诺稳定库接口。
- `ChromSequence` 字段为 `pub(crate)`，外部仅可通过 `name()`, `sequence()`, `into_parts()` 访问；`SequenceData`、`SequenceSlice` 等辅助类型保持 crate-private。

## 核心搜索流程
1. 原始 FASTA 字节归一化为 `SequenceData`（小写 `Arc<Vec<u8>>`）。
2. `GRunScanner` 使用 `memchr2` 查找 G-run，生成 BFS 种子。
3. `G4Candidate::expand` 穷举三个 loop 长度（`find_loop_lengths_from`），确保 `partial_length() <= maximum_length()`。
4. `G4Candidate::score()` 保留与 C++ 一致的 `floor(gmax - gavg + gmax*(tetrads-2))`。
5. `find_raw_on_window_bytes` / `find_raw_bytes_no_chunking` 生成 raw 命中；`find_owned_bytes*` 仅负责组装这些 raw hits，需由调用方显式传给 `consolidate_g4s` 做去重/家族合并。
6. 结果通过 `render_csv_results` 或 `write_parquet_results` 输出。

### 分块策略
- `chunk_size_for_limits()` 基于 `ScanLimits.max_g4_length` + padding，范围固定在 32~64bp。
- `compute_chunk_overlap()` 始终返回 `max_g4_length`，避免窗口边缘截断。
- 大序列拆成 `(start, primary_end, window_end)`，使用 Rayon `into_par_iter().flat_map_iter()` 聚合；短序列直接调用 `find_with_sequence()`。

### 家族合并 (`consolidation.rs`)
- `DedupKey = (start, end, SequenceSlice)`，保留最高 gscore；
- 去重结果按 `(start, end)` 排序后线性扫描，`belongs_in` 判断是否重叠；
- 每个家族输出最高 gscore 的成员，以确保 deterministic 输出。

## Streaming 管线 (`src/qgrs/stream.rs`)
`StreamChromosome` 管理每条染色体，`StreamChunkScheduler` 维护缓冲区、offset 与 worker 信道：
- FASTA 逐行读取，非序列字符跳过并转小写；
- 缓冲长度到达 `chunk_size + overlap` 即调度 worker，worker 直接调用 `find_raw_bytes_no_chunking()`；
- 染色体结束后调用 `finish()` 返回去重后的结果。

## CLI (`src/bin/qgrs.rs`)
- Inline `--sequence`：调用 `find_owned_bytes_with_limits()` 获取 raw hits，再交给 `consolidate_g4s()`，CSV 默认写 stdout，Parquet 需 `--output`。
- `--file` + `--mode mmap`：`load_sequences_from_path()` → Rayon 并行 → CSV/Parquet 分染色体写入。
- `--file` + `--mode stream`：`stream::process_fasta_stream_with_limits()` 回调中写文件。
- 线程数固定为 `num_cpus::get()`，避免依赖 `RAYON_NUM_THREADS`。

## `.rs` 文件速查表
| 文件 | 说明 |
| --- | --- |
| `src/lib.rs` | 仅 `pub mod qgrs;`，避免 crate 根 API 泄露。 |
| `src/qgrs/mod.rs` | 模块入口：声明 `pub mod stream;`，`pub use` 暴露搜索/导出 API。 |
| `src/qgrs/data.rs` | `InputMode`, `ChromSequence`, `ScanLimits`, `SequenceData`, `SequenceSlice`, `is_g` 等基础数据结构。 |
| `src/qgrs/search.rs` | BFS 搜索实现：`G4`, `G4Candidate`, `GRunScanner`, `find_raw_*` 等核心算法。 |
| `src/qgrs/chunks.rs` | `find_owned_bytes*`, chunk/overlap 计算、`find_with_sequence`、`shift_g4`。负责串联 Rayon 窗口并返回 raw hits（调用者自行 `consolidate_g4s`）。 |
| `src/qgrs/consolidation.rs` | `consolidate_g4s`, `DedupKey`, overlap 判断、家族 winner 逻辑。 |
| `src/qgrs/export.rs` | `render_csv_results`, `write_parquet_results`, `ExportError`，包含 CSV 转义与 Arrow/Parquet 写入。 |
| `src/qgrs/loaders.rs` | `load_sequences_from_path`, `load_sequences_mmap`, `load_sequences_stream`, `parse_chrom_name` 以及内部 push helper。 |
| `src/qgrs/stream.rs` | Streaming FASTA 处理：`process_fasta_stream(_with_limits)`, `process_reader`, `StreamChromosome`, `StreamChunkScheduler`。 |
| `src/qgrs/tests/helpers.rs` | 测试辅助（例如 `arc_from_sequence`）。 |
| `src/qgrs/tests/unit.rs` | 单测：命中检测、CSV/Parquet 输出、loader 行为。 |
| `src/qgrs/tests/integration_chunk.rs` | Chunk 与非 chunk 结果一致性、边界情况。 |
| `src/qgrs/tests/integration_stream.rs` | Stream 与 batch 结果等价验证。 |
| `src/qgrs/tests/mod.rs` | 组织 helper + 子模块。 |
| `src/bin/qgrs.rs` | 主 CLI，解析参数、调度 `qgrs` 模块并写 CSV/Parquet。 |
| `src/bin/compare_modes.rs` | 开发用工具：比较 mmap vs stream 的性能与命中差异。 |
| `src/bin/compare_csv_outputs.rs` | 对比两个输出目录的 CSV 是否逐行一致，打印差异详情。 |

## 关键约定与陷阱
1. **坐标体系**：内部 0-based 半开区间 `[start,end)`；输出 start+1、end inclusive，CSV/Parquet 必须保持一致。
2. **排序后合并**：`consolidate_g4s` 在排序后的 raw hits 上工作，确保家族顺序 deterministic。
3. **Stream worker 禁止再分块**：`StreamChunkScheduler` 已提供 overlap，worker 只能直接跑 `find_raw_bytes_no_chunking()`。
4. **避免调试输出**：库函数不得打印 stdout；如需调试请加 feature flag 或日志。
5. **Arc clone 便宜**：`Arc<Vec<u8>>` clone 仅递增引用计数，可在 Rayon worker 中放心复用。
6. **ChromSequence 访问**：外部如需名字或序列，使用 `chrom.name()` / `chrom.sequence()` / `chrom.into_parts()`。

## 构建、测试与常用脚本
```bash
# Release 构建 CLI
cargo build --release --bin qgrs

# 单元 + 集成测试
cargo test --lib

# mmap vs stream 性能/一致性
target/release/compare_modes data.fa 3 20

# CSV 目录对比
target/release/compare_csv_outputs out_mmap out_stream
```
> `compare_modes` 既汇报性能也会在发现差异时列出前若干条，适合作为轻量基准脚本。

## 修改核心逻辑的 Checklist
- 任何涉及搜索、分块或去重的改动都要跑 `cargo test --lib`，重点关注 `big_sequence_internal_equals_chunked` 和 `stream_pipeline_matches_batch_results`；
- CLI/输出格式改动：同步更新 README、`--help`，并使用 `compare_csv_outputs`/`compare_modes` 做回归；
- 新增输出格式：扩展 `OutputFormat`、实现渲染器、更新 CLI 两条路径以及测试；
- 调整 `ScanLimits`/chunk 行为时记得刷新 README 架构描述并记录默认参数；
- Rayon 线程数固定在 `main()`，避免依赖环境变量导致不确定行为。
