## 计划：实时染色体进度展示

内置一条进度通道，让 mmap 与 stream 两种路径都能输出“染色体名称 + 当前 seed 进度百分比”，并保证不会污染 CSV/Parquet 结果。实现方式是在 `src/bin/qgrs.rs` 中创建轻量的 reporter 回调，向下传递到 `src/qgrs/mod.rs` 与 `src/qgrs/stream.rs` 的扫描流程中：每当某条染色体开始处理或 seed 起始位置前进时就触发事件。所有进度信息统一打印到 stderr，并通过 `\r` 回写与节流实现终端动态刷新，避免 Rayon 并行线程产生噪声。

### 实施步骤

1. 在 `src/bin/qgrs.rs` 内实现默认启用的 reporter，使用 `\r` 重写同一行向 stderr 输出 `name: current/total (percent)`，实现终端动态刷新；同时更新 `usage` 文案说明实时进度特性。
2. 在 `src/lib.rs` 定义 `ProgressReporter` trait（`start`、`seed_progress`、`finish`），并提供 `NoopReporter`；对外接口使用 `Arc<dyn ProgressReporter + Send + Sync>` 便于并发共享。
3. 调整 mmap 流程（`process_inline_sequence`、`process_fasta_file`、`qgrs::find_csv_owned_with_limits` 等）接受 reporter：扫描开始前调用 `start`，枚举 seed 时调用 `seed_progress(name, seed_start, total_len)`，结束后调用 `finish`。
4. 调整 stream 流程（`qgrs::stream::process_fasta_stream_with_limits`），将 reporter 贯穿 chunk 调度与回调，使得每条染色体在 chunk 起止推进时能报告总长度上的进度。
5. 通过 `Arc` + 互斥锁或 `crossbeam-channel` 等方式安全地从 Rayon 线程发送进度事件，并确保所有日志只写入 stderr，避免与 CSV/Parquet 输出混合。

### 额外注意

1. 确保进度信息永远与结果输出分离（仅写 stderr），并为极短序列设置最小刷新间隔，防止日志频率过高。
