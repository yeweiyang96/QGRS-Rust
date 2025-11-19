**行动框架**

1. 现有依赖：`memmap2` 负责零拷贝访问、`memchr` 在 ARM64 上自带 NEON SIMD；它们已经自动利用 M1 的指令集与预取机制，所以不必改动。
2. 识别瓶颈：当 FASTA 只有一个染色体（例如 aaa.fa）时，所有工作都落在单线程的 `seed_queue`/`G4Candidate` 扩展上，CPU 核心无法同时工作。
3. 设计并行策略：保持单个映射视图不变的前提下拆分工作负载，让多个核心在不同窗口上并行扫描，再合并候选。

**加速建议**

- **Release 编译 + M1 优化开关**：确保执行 `cargo run --release`，并在 `.cargo/config.toml` 中加入

  ```toml
  [build]
  rustflags = ["-C", "target-cpu=native"]
  ```

  这样 LLVM 会针对 Apple M1 的 big.LITTLE 架构生成 NEON 指令，提高 `memchr`、分支预测等底层性能。

- **按染色体并行**：如果 FASTA 含多条染色体，直接用 `rayon::par_iter()` 包裹 `load_sequences_from_path` 返回的 `Vec<ChromSequence>`，每个条目独立运行 `find_owned` 并写出结果即可——数据完全隔离，无需同步。当前 CLI 可以在 `process_fasta_file` 的 mmap 分支中加入：

  ```rust
  sequences.par_iter().try_for_each(|chrom| { ... })
  ```

  这样多核在多个染色体上同时工作。

- **单染色体拆分窗口**：对于像 aaa.fa 这种只有一条序列的情况，可以把序列按固定大小（例如 2–4 MB）划成块，并在块之间加入 `(max_tetrads * 4 + max_loop)` 的重叠，以免跨块 G-run 被截断。每个块用 `rayon::scope` 派发到线程池运行 `find_slice(chunk_range)`, 最后把所有候选合并再执行 `consolidate_g4s` 去重。实现要点：

  1. 预先计算块边界并保留 `start_offset`，传入 `find_with_sequence` 的变体，使候选位置保持全局坐标。
  2. 合并阶段按 `start` 排序，复用现有 `consolidate_g4s`。
  3. 大序列 memory-map 后仍旧是同一个 `&[u8]` 视图，不会额外占内存。

- **流式模式的并行**：`stream` 目前是单线程滑动窗口。若要并行，可在读取器层面把 FASTA 分块（同样带重叠）后用 `rayon::ThreadPool` 为每个块启动一个新的 `StreamDetector`，再归并结果。由于窗口依赖前缀，下一个块必须携带前一块的尾巴才能保持准确度。

- **多进程方案**：若 Rust 端修改成本太高，可以在 shell 层面把超长 FASTA 切片（`split -l`）后并行启动多个 `qgrs` 实例，每个实例处理不同区段并写出临时结果，再在外部脚本中归并。M1 有 4 颗性能核，可同时跑 4 个进程。

总之，`memmap2`/`memchr` 已自动利用 M1 的硬件；要进一步提速，需要显式在 G4 搜索阶段引入 `rayon` 或多进程，把单染色体的扫描任务拆成可并发的窗口。你可以先尝试最简单的“多染色体并行 + release + target-cpu=native”，再视需要实现单染色体窗口化。
