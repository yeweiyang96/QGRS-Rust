在不丢失 hits 的前提下，通过缓存、剪枝以及线程级复用降低 find_raw_on_window_bytes / find_raw_with_sequence 的分配与扩展开销；允许在 search.rs 内引入 Rayon 或线程本地缓存以支撑更复杂的优化。

## Steps
1. 在 search.rs 顶部补充约束注释（chunk/overlap、0-based 坐标、禁止 worker 再分块、consolidate_g4s 责任），确保后续复杂改动不破坏既定规则。
2. 修改 qgrs/chunks.rs, qgrs/stream.rs, search.rs 之间的接口：缓存 Arc<SequenceData> 并传递引用，让 find_raw_on_window_bytes 只处理切片与 offset，避免每个窗口重复创建 SequenceData::from_bytes。
3. 在 search.rs 中增强候选剪枝：基于 primary_end、max_length 和 window 尺寸限制 tetrads/max_offset，在种子阶段即丢弃必超长或重复区域的候选，减少 BFS 队列规模。
4. 为 G4Candidate::expand 与 BFS 队列引入线程本地缓存：使用 Rayon thread-local (或 std::thread::LocalKey) 复用 Vec<i32>、对象池或 slab，减少 VecDeque push/pop 的分配，并评估改用 Vec + 索引环的可行性。
5. 延迟 G4::from_candidate 的昂贵部分：仅在 cand.viable 通过后构建 SequenceSlice；必要时引入轻量结构暂存 start/length 并在最终阶段创建视图，保持 find_raw_with_sequence、chunk、stream 路径的命中一致性。