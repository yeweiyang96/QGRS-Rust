**优化目标**

1. 统一“最大 G-run 长度 = 10 / 最大 G4 长度 = 45”的约束，并允许通过 CLI 参数 `--max-g-run`、`--max-g4-length` 自定义，默认 10/45。
2. 把 mmap 与 stream 两条管线都限制在 ≤64bp 的小窗口内，以提升多核利用率、减少无效扫描。
3. 在种子、候选扩展与合并阶段全面裁剪超出约束的情况，确保计算仅集中在潜在 G4 区域。
4. 对 CLI 输入做动态校验：`max_g4_length ≥ 4 * min_tetrads` 且 `max_g_run ≥ min_tetrads`，否则提前报错并提示用户。

**实施步骤**

1. 常量与参数结构

- 在 `src/qgrs/mod.rs` 定义 `DEFAULT_MAX_G4_LENGTH = 45`、`DEFAULT_MAX_G_RUN = 10`。
- 为 `find`/`find_owned`、chunk/stream 入口添加 `(max_g4_length, max_g_run)` 参数；默认使用常量，CLI 可覆盖。

2. Clamp Chunk Windows（≤64bp）

- 将 `CHUNK_SIZE` 调整为 `max_g4_length + safety`（建议 64）。
- `compute_chunk_overlap`、stream scheduler overlap 统一采用 `max_g4_length + padding`，确保跨 chunk 的候选不会被截断。

3. Seed Filter

- 修改 `GRunScanner` / `seed_queue`：忽略连续 G-run 长度 > `max_g_run` 的片段。
- 只枚举满足 `4 * tetrads ≤ max_g4_length` 的 tetrad 组合，过长的直接跳过。

4. Candidate Guardrails

- `G4Candidate::find_loop_lengths_from` / `expand` 在 `partial_length > max_g4_length` 时立即停止扩展；loop 搜索区间基于剩余预算裁剪。
- 若某个 tetrad 区段长度会超过 `max_g_run`，提前放弃该分支。

5. Scoring & Consolidation Adjustments

- `find_owned_chunked`、stream chunk merge 全部使用新的窗口/overlap 参数。
- 复核 `consolidate_g4s` 在小窗口场景仍可正确归并；必要时在合并前按 `start` 排序。

6. CLI 校验与文档

- `src/bin/qgrs.rs` 新增 `--max-g4-length`、`--max-g-run`；解析后校验 `max_g4_length ≥ 4 * min_tetrads`、`max_g_run ≥ min_tetrads`，否则返回 usage。
- README/usage 补充新参数说明及推荐范围（如 32–64bp / 8–12bp）。

7. 测试与基准

- 单测覆盖：默认参数、自定义合法参数、非法输入报错、chunk/stream 结果一致性。
- `cargo test` 全量回归。
- 基准建议：对 `aaa.fa` 分别运行 `--max-g4-length 32|48|64`、`--max-g-run 8|10|12`，记录 mmap/stream 的 real/user/RSS 指标，为 M1 Pro (64GB) 选出最佳配置。

**验证指标**

- CLI 解析：错误输入能即时报错；默认值与旧版行为一致。
- 扫描性能：stream 模式相对旧版至少提升 ~20%，内存占用下降；mmap 模式在长染色体场景下接近 stream 性能。
- 结果一致性：不同参数组合下，`stream` 与 `find_owned` 仍给出相同的 G4 集合。
