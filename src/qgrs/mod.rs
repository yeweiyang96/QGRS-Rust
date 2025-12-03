use std::collections::VecDeque;
use std::fmt;
use std::fs::File;
use std::hash::Hash;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, OnceLock};

use arrow_array::{ArrayRef, Int32Array, RecordBatch, StringArray, UInt64Array};
use arrow_schema::{DataType, Field, Schema};
use memchr::memchr2;
use memmap2::MmapOptions;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::errors::ParquetError;
use rayon::prelude::*;
// debug helpers removed for cleanup

pub mod stream;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InputMode {
    Mmap,
    Stream,
}

#[derive(Clone, Debug)]
pub struct ChromSequence {
    pub name: String,
    pub sequence: Arc<Vec<u8>>,
}

impl ChromSequence {
    pub fn as_uppercase_string(&self) -> String {
        let mut seq = unsafe { String::from_utf8_unchecked(self.sequence.as_ref().clone()) };
        seq.make_ascii_uppercase();
        seq
    }
}

pub const DEFAULT_MAX_G4_LENGTH: usize = 45;
pub const DEFAULT_MAX_G_RUN: usize = 10;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ScanLimits {
    pub max_g4_length: usize,
    pub max_g_run: usize,
}

impl ScanLimits {
    pub const fn new(max_g4_length: usize, max_g_run: usize) -> Self {
        Self {
            max_g4_length,
            max_g_run,
        }
    }
}

impl Default for ScanLimits {
    fn default() -> Self {
        ScanLimits::new(DEFAULT_MAX_G4_LENGTH, DEFAULT_MAX_G_RUN)
    }
}

const WINDOW_MIN_BP: usize = 32;
const WINDOW_MAX_BP: usize = 64;
const WINDOW_PADDING_BP: usize = 27;
// Increased overlap margin for stream chunking experiments (was 16).
// Larger overlap helps ensure families spanning chunk boundaries are scanned
// within a single worker's window, reducing divergent merged boundaries.

// Debug targets: sequences (lowercase) or start positions we want to trace.
// debug helpers removed

#[derive(Clone, Debug)]
struct SequenceData {
    normalized: Arc<Vec<u8>>,
}

impl SequenceData {
    #[cfg_attr(not(test), allow(dead_code))]
    fn new(sequence: &str) -> Self {
        let normalized = Arc::new(sequence.to_ascii_lowercase().into_bytes());
        Self { normalized }
    }

    fn from_bytes(normalized: Arc<Vec<u8>>) -> Self {
        Self { normalized }
    }
}

#[derive(Clone, Debug)]
struct SequenceSlice {
    normalized: Arc<Vec<u8>>,
    start: usize,
    length: usize,
}

impl SequenceSlice {
    fn new(normalized: Arc<Vec<u8>>, start: usize, length: usize) -> Self {
        Self {
            normalized,
            start,
            length,
        }
    }

    fn bytes(&self) -> &[u8] {
        let end = self.start + self.length;
        &self.normalized[self.start..end]
    }

    fn to_uppercase_string(&self) -> String {
        let mut sequence = unsafe { String::from_utf8_unchecked(self.bytes().to_vec()) };
        sequence.make_ascii_uppercase();
        sequence
    }
}

impl PartialEq for SequenceSlice {
    fn eq(&self, other: &Self) -> bool {
        self.length == other.length && self.bytes() == other.bytes()
    }
}

impl Eq for SequenceSlice {}

impl std::hash::Hash for SequenceSlice {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.length.hash(state);
        state.write(self.bytes());
    }
}

#[derive(Eq, PartialEq, Hash)]
struct DedupKey {
    start: usize,
    end: usize,
    slice: SequenceSlice,
}

impl DedupKey {
    fn new(g4: &G4) -> Self {
        Self {
            start: g4.start,
            end: g4.end,
            slice: g4.sequence_view.clone(),
        }
    }
}

struct GRunScanner<'a> {
    data: &'a [u8],
    cursor: usize,
    min_tetrads: usize,
    max_g_run: usize,
}

impl<'a> GRunScanner<'a> {
    fn new(data: &'a [u8], min_tetrads: usize, max_g_run: usize) -> Self {
        Self {
            data,
            cursor: 0,
            min_tetrads,
            max_g_run,
        }
    }
}

impl<'a> Iterator for GRunScanner<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let len = self.data.len();
        while self.cursor < len {
            let search_slice = &self.data[self.cursor..];
            let relative = memchr2(b'g', b'G', search_slice)?;
            let run_start = self.cursor + relative;
            let mut run_end = run_start;
            while run_end < len && is_g(unsafe { *self.data.get_unchecked(run_end) }) {
                run_end += 1;
            }
            self.cursor = if run_end < len { run_end + 1 } else { len };
            let run_len = run_end - run_start;
            if run_len >= self.min_tetrads && run_len <= self.max_g_run {
                return Some((run_start, run_len));
            }
        }
        None
    }
}

#[inline(always)]
pub(super) fn is_g(byte: u8) -> bool {
    byte == b'G' || byte == b'g'
}

#[derive(Clone, Debug)]
struct G4Candidate {
    seq: Arc<SequenceData>,
    num_tetrads: usize,
    start: usize,
    y1: i32,
    y2: i32,
    y3: i32,
    max_length: usize,
}

impl G4Candidate {
    fn new(seq: Arc<SequenceData>, num_tetrads: usize, start: usize, limits: ScanLimits) -> Self {
        Self {
            seq,
            num_tetrads,
            start,
            y1: -1,
            y2: -1,
            y3: -1,
            max_length: maximum_length(num_tetrads, limits),
        }
    }

    fn score(&self) -> i32 {
        let gavg = (f64::from((self.y1 - self.y2).abs())
            + f64::from((self.y2 - self.y3).abs())
            + f64::from((self.y1 - self.y3).abs()))
            / 3.0;
        let gmax = (self.max_length as i32 - (self.num_tetrads as i32 * 4 + 1)) as f64;
        let bonus = gmax * ((self.num_tetrads as i32 - 2) as f64);
        (gmax - gavg + bonus).floor() as i32
    }

    fn length(&self) -> usize {
        (4 * self.num_tetrads)
            + self.y1.max(0) as usize
            + self.y2.max(0) as usize
            + self.y3.max(0) as usize
    }

    fn t1(&self) -> usize {
        self.start
    }

    fn t2(&self) -> usize {
        self.t1() + self.num_tetrads + self.y1.max(0) as usize
    }

    fn t3(&self) -> usize {
        self.t2() + self.num_tetrads + self.y2.max(0) as usize
    }

    fn t4(&self) -> usize {
        self.t3() + self.num_tetrads + self.y3.max(0) as usize
    }

    fn cursor(&self) -> Option<usize> {
        if self.y1 < 0 {
            Some(self.t1() + self.num_tetrads)
        } else if self.y2 < 0 {
            Some(self.t2() + self.num_tetrads)
        } else if self.y3 < 0 {
            Some(self.t3() + self.num_tetrads)
        } else {
            None
        }
    }

    fn partial_length(&self) -> i32 {
        let mut length = (self.num_tetrads * 4) as i32;
        if self.y1 >= 0 && self.y2 < 0 {
            length += if self.y1 == 0 { 2 } else { 1 };
        } else if self.y2 >= 0 && self.y3 < 0 {
            length += if self.y1 == 0 || self.y2 == 0 { 1 } else { 0 };
        }
        if self.y1 > 0 {
            length += self.y1;
        }
        if self.y2 > 0 {
            length += self.y2;
        }
        if self.y3 > 0 {
            length += self.y3;
        }
        length
    }

    fn min_acceptable_loop_length(&self) -> i32 {
        if self.y1 == 0 || self.y2 == 0 || self.y3 == 0 {
            1
        } else {
            0
        }
    }

    fn complete(&self) -> bool {
        self.y1 >= 0 && self.y2 >= 0 && self.y3 >= 0
    }

    fn viable(&self, min_score: i32) -> bool {
        if self.score() < min_score {
            return false;
        }
        if self.length() > self.max_length {
            return false;
        }
        let mut zero_loops = 0;
        if self.y1 < 1 {
            zero_loops += 1;
        }
        if self.y2 < 1 {
            zero_loops += 1;
        }
        if self.y3 < 1 {
            zero_loops += 1;
        }
        zero_loops < 2
    }

    fn find_loop_lengths_from(&self, ys: &mut Vec<i32>, cursor: usize) {
        let mut p = cursor;
        let seq = &self.seq.normalized;
        let max_pos = self.start + self.max_length + 1;
        let target_len = self.num_tetrads;
        let min_loop = self.min_acceptable_loop_length();

        while p + target_len <= seq.len() {
            if p >= max_pos {
                break;
            }
            if seq[p..p + target_len].iter().all(|&b| b == b'g') {
                let y = (p - cursor) as i32;
                if y >= min_loop && (p - self.start + target_len - 1) < self.max_length {
                    ys.push(y);
                } else {
                    break;
                }
            }
            p += 1;
        }
    }

    fn expand(&self) -> Vec<G4Candidate> {
        let mut results = Vec::new();
        if let Some(cursor) = self.cursor() {
            let mut ys = Vec::new();
            self.find_loop_lengths_from(&mut ys, cursor);
            for y in ys {
                let mut next = self.clone();
                if next.y1 < 0 {
                    next.y1 = y;
                } else if next.y2 < 0 {
                    next.y2 = y;
                } else if next.y3 < 0 {
                    next.y3 = y;
                }
                if next.partial_length() <= next.max_length as i32 {
                    results.push(next);
                }
            }
        }
        results
    }
}

fn seed_queue(
    cands: &mut VecDeque<G4Candidate>,
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    limits: ScanLimits,
) {
    let mut max_tetrads_allowed = limits.max_g_run;
    if limits.max_g4_length >= 4 {
        max_tetrads_allowed = max_tetrads_allowed.min(limits.max_g4_length / 4);
    }
    if max_tetrads_allowed < min_tetrads {
        return;
    }
    for (run_start, run_len) in GRunScanner::new(&seq.normalized, min_tetrads, limits.max_g_run) {
        let capped_run_len = run_len.min(max_tetrads_allowed);
        let mut tetrads = min_tetrads;
        while tetrads <= capped_run_len {
            if tetrads * 4 > limits.max_g4_length {
                break;
            }
            let max_offset = capped_run_len - tetrads;
            for offset in 0..=max_offset {
                let start = run_start + offset;
                cands.push_back(G4Candidate::new(seq.clone(), tetrads, start, limits));
            }
            tetrads += 1;
        }
    }
}

#[derive(Debug)]
pub struct G4 {
    pub start: usize,
    pub end: usize,
    pub tetrad1: usize,
    pub tetrad2: usize,
    pub tetrad3: usize,
    pub tetrad4: usize,
    pub y1: i32,
    pub y2: i32,
    pub y3: i32,
    pub tetrads: usize,
    pub length: usize,
    pub gscore: i32,
    sequence_view: SequenceSlice,
    sequence_cache: OnceLock<String>,
}

impl G4 {
    fn from_candidate(candidate: &G4Candidate) -> Self {
        let length = candidate.length();
        let end = candidate.start + length;
        let normalized = candidate.seq.normalized.clone();
        let sequence_view = SequenceSlice::new(normalized, candidate.start, length);
        Self {
            start: candidate.start + 1,
            end,
            tetrad1: candidate.t1() + 1,
            tetrad2: candidate.t2() + 1,
            tetrad3: candidate.t3() + 1,
            tetrad4: candidate.t4() + 1,
            y1: candidate.y1,
            y2: candidate.y2,
            y3: candidate.y3,
            tetrads: candidate.num_tetrads,
            length,
            gscore: candidate.score(),
            sequence_view,
            sequence_cache: OnceLock::new(),
        }
    }

    pub fn sequence(&self) -> &str {
        self.sequence_cache
            .get_or_init(|| self.sequence_view.to_uppercase_string())
    }
}

impl Clone for G4 {
    fn clone(&self) -> Self {
        Self {
            start: self.start,
            end: self.end,
            tetrad1: self.tetrad1,
            tetrad2: self.tetrad2,
            tetrad3: self.tetrad3,
            tetrad4: self.tetrad4,
            y1: self.y1,
            y2: self.y2,
            y3: self.y3,
            tetrads: self.tetrads,
            length: self.length,
            gscore: self.gscore,
            sequence_view: self.sequence_view.clone(),
            sequence_cache: OnceLock::new(),
        }
    }
}

pub(super) fn maximum_length(num_tetrads: usize, limits: ScanLimits) -> usize {
    let base = if num_tetrads < 3 { 30 } else { 45 };
    base.min(limits.max_g4_length)
}

fn overlapped(a: &G4, b: &G4) -> bool {
    let a_start = a.start as isize;
    let a_end = (a.start + a.length) as isize;
    let b_start = b.start as isize;
    let b_end = (b.start + b.length) as isize;
    (a_start >= b_start && a_start <= b_end)
        || (a_end >= b_start && a_end <= b_end)
        || (b_start >= a_start && b_start <= a_end)
        || (b_end >= a_start && b_end <= a_end)
}

fn belongs_in(g4: &G4, family: &[G4]) -> bool {
    family.iter().any(|member| overlapped(g4, member))
}

pub(super) fn consolidate_g4s(mut raw_g4s: Vec<G4>) -> Vec<G4> {
    // First, sort by start so grouping is stable.
    raw_g4s.sort_by(|a, b| a.start.cmp(&b.start));

    // Remove exact-duplicate candidates produced from overlapping chunks，
    // 同一 key (start/end/序列) 只保留最高 gscore 的成员，避免早期插入
    // 的低分候选把高分版本顶掉。
    use std::collections::HashMap;
    use std::collections::hash_map::Entry;
    let mut best_by_key: HashMap<DedupKey, G4> = HashMap::new();
    // 消费 raw_g4s，转移所有权（避免克隆大量数据）
    for g in raw_g4s.into_iter() {
        let key = DedupKey::new(&g);
        match best_by_key.entry(key) {
            Entry::Vacant(slot) => {
                slot.insert(g);
            }
            Entry::Occupied(mut slot) => {
                if g.gscore > slot.get().gscore {
                    slot.insert(g);
                }
            }
        }
    }
    let mut deduped: Vec<G4> = best_by_key.into_values().collect();
    // HashMap 的迭代顺序是不确定的（取决于哈希值和内部扩容）
    // 后续家族分组（belongs_in）是顺序敏感的：候选遇到第一个重叠家族就加入并停止
    // 不同顺序会导致同一批候选被分到不同的家族
    deduped.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));

    let mut families: Vec<Vec<G4>> = Vec::new();

    for g4 in deduped.into_iter() {
        let mut inserted = false;
        for family in &mut families {
            if belongs_in(&g4, family) {
                family.push(g4.clone());
                inserted = true;
                break;
            }
        }
        if !inserted {
            families.push(vec![g4]);
        }
    }

    let mut results = Vec::new();
    for family in families.into_iter() {
        if family.is_empty() {
            continue;
        }
        // Normal family consolidation: choose best by gscore
        let mut best = family[0].clone();
        for member in &family {
            if member.gscore > best.gscore {
                best = member.clone();
            }
        }
        // end per-family debug removed
        results.push(best);
    }

    results
}

pub fn find_owned_bytes(sequence: Arc<Vec<u8>>, min_tetrads: usize, min_score: i32) -> Vec<G4> {
    find_owned_bytes_with_limits(sequence, min_tetrads, min_score, ScanLimits::default())
}

pub fn find_owned_bytes_with_limits(
    sequence: Arc<Vec<u8>>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    let chunk_size = chunk_size_for_limits(limits);
    if sequence.len() > chunk_size {
        // 零拷贝窗口策略：窗口仅决定种子扫描范围，候选扩展沿用
        // 全局 Arc<Vec<u8>>，保证 family 的最佳成员不会因为 chunk 边界
        // 被截断；overlap 仅提供延伸上下文。
        let len = sequence.len();
        let overlap = compute_chunk_overlap(min_tetrads, limits);
        let mut start = 0usize;
        let seq_arc = sequence.clone();
        let windows: Vec<(usize, usize, usize)> = {
            let mut v = Vec::new();
            while start < len {
                let primary_end = (start + chunk_size).min(len);
                let window_end = (primary_end + overlap).min(len);
                v.push((start, primary_end, window_end));
                start = primary_end;
            }
            v
        };
        let merged_raw: Vec<G4> = windows
            .into_par_iter()
            .flat_map_iter(|(offset, primary_end, window_end)| {
                let window_slice = &seq_arc[offset..window_end];
                let hits = find_raw_on_window_bytes(
                    seq_arc.clone(),
                    offset,
                    primary_end,
                    window_slice,
                    min_tetrads,
                    min_score,
                    limits,
                );
                hits.into_iter()
            })
            .collect();
        return consolidate_g4s(merged_raw);
    }
    // Non-chunked path: build SequenceData from bytes and call find_with_sequence.
    let seq = Arc::new(SequenceData::from_bytes(sequence));
    find_with_sequence(seq, min_tetrads, min_score, limits)
}

pub(super) fn find_raw_bytes_no_chunking(
    mut sequence: Vec<u8>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    sequence.make_ascii_lowercase();
    let seq = Arc::new(SequenceData::from_bytes(Arc::new(sequence)));
    find_raw_with_sequence(seq, min_tetrads, min_score, limits)
}

// Zero-copy windowed raw search: scan `window` (a slice of the full sequence)
// but build `G4Candidate` entries that reference the full underlying
// sequence via `seq_bytes_arc` so final `G4` outputs have global coordinates
// and can reuse `G4::from_candidate` without copying the whole input.
#[allow(dead_code)]
fn find_raw_on_window_bytes(
    seq_bytes_arc: Arc<Vec<u8>>,
    base_offset: usize,
    primary_end: usize,
    window: &[u8],
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    // Create a SequenceData that references the full normalized bytes (no copy).
    let seq = Arc::new(SequenceData::from_bytes(seq_bytes_arc));
    // Seed queue using GRunScanner over the provided window slice。种子
    // 只在主 chunk 范围内生成（run_start < primary_end），而 loop 扩展
    // 可访问完整序列，避免窗口字符串长度限制。
    let mut cands = VecDeque::new();
    let mut max_tetrads_allowed = limits.max_g_run;
    if limits.max_g4_length >= 4 {
        max_tetrads_allowed = max_tetrads_allowed.min(limits.max_g4_length / 4);
    }
    if max_tetrads_allowed >= min_tetrads {
        for (run_start_rel, run_len) in GRunScanner::new(window, min_tetrads, limits.max_g_run) {
            let run_start = base_offset + run_start_rel;
            if run_start >= primary_end {
                continue;
            }
            let capped_run_len = run_len.min(max_tetrads_allowed);
            let mut tetrads = min_tetrads;
            while tetrads <= capped_run_len {
                if tetrads * 4 > limits.max_g4_length {
                    break;
                }
                let max_offset = capped_run_len - tetrads;
                for offset in 0..=max_offset {
                    let start = run_start + offset;
                    cands.push_back(G4Candidate::new(seq.clone(), tetrads, start, limits));
                }
                tetrads += 1;
            }
        }
    }

    let mut raw_g4s = Vec::new();
    while let Some(cand) = cands.pop_front() {
        if cand.complete() {
            if cand.viable(min_score) {
                let g4 = G4::from_candidate(&cand);
                raw_g4s.push(g4);
            }
        } else {
            for expanded in cand.expand() {
                cands.push_back(expanded);
            }
        }
    }
    raw_g4s
}

pub(crate) fn chunk_size_for_limits(limits: ScanLimits) -> usize {
    let desired = limits.max_g4_length.saturating_add(WINDOW_PADDING_BP);
    desired.clamp(WINDOW_MIN_BP, WINDOW_MAX_BP)
}

pub(crate) fn compute_chunk_overlap(_min_tetrads: usize, limits: ScanLimits) -> usize {
    // overlap 必须至少覆盖一个完整的最大 G4，确保跨 chunk 的候选
    // 能在单次窗口中完成扩展与评分。
    limits.max_g4_length.max(1)
}

pub(crate) fn shift_g4(g4: &mut G4, offset: usize) {
    g4.start += offset;
    g4.end += offset;
    g4.tetrad1 += offset;
    g4.tetrad2 += offset;
    g4.tetrad3 += offset;
    g4.tetrad4 += offset;
}

fn find_with_sequence(
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    // Find raw candidates first
    let raw = find_raw_with_sequence(seq, min_tetrads, min_score, limits);
    print!("{}", raw.len());
    // Then consolidate families
    consolidate_g4s(raw)
}

fn find_raw_with_sequence(
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    let mut cands = VecDeque::new();
    seed_queue(&mut cands, seq.clone(), min_tetrads, limits);
    let mut raw_g4s = Vec::new();
    while let Some(cand) = cands.pop_front() {
        if cand.complete() {
            if cand.viable(min_score) {
                let g4 = G4::from_candidate(&cand);
                raw_g4s.push(g4);
            }
        } else {
            for expanded in cand.expand() {
                cands.push_back(expanded);
            }
        }
    }
    raw_g4s
}

pub fn render_csv_results(g4s: &[G4]) -> String {
    render_csv(g4s)
}

fn render_csv(g4s: &[G4]) -> String {
    let mut out = String::from("start,end,length,tetrads,y1,y2,y3,gscore,sequence\n");
    for g4 in g4s {
        let sequence_field = escape_csv_field(g4.sequence());
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{}\n",
            g4.start, g4.end, g4.length, g4.tetrads, g4.y1, g4.y2, g4.y3, g4.gscore, sequence_field
        ));
    }
    out
}

fn escape_csv_field(value: &str) -> String {
    if value.is_empty() {
        return String::new();
    }
    let needs_quotes = value.contains([',', '"', '\n']);
    if !needs_quotes {
        return value.to_string();
    }
    let mut escaped = String::from('"');
    for ch in value.chars() {
        if ch == '"' {
            escaped.push_str("\"\"");
        } else {
            escaped.push(ch);
        }
    }
    escaped.push('"');
    escaped
}

#[derive(Debug)]
pub enum ExportError {
    Arrow(arrow_schema::ArrowError),
    Parquet(ParquetError),
}

impl From<arrow_schema::ArrowError> for ExportError {
    fn from(value: arrow_schema::ArrowError) -> Self {
        ExportError::Arrow(value)
    }
}

impl From<ParquetError> for ExportError {
    fn from(value: ParquetError) -> Self {
        ExportError::Parquet(value)
    }
}

impl fmt::Display for ExportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExportError::Arrow(err) => write!(f, "arrow error: {err}"),
            ExportError::Parquet(err) => write!(f, "parquet error: {err}"),
        }
    }
}

impl std::error::Error for ExportError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ExportError::Arrow(err) => Some(err),
            ExportError::Parquet(err) => Some(err),
        }
    }
}

pub fn write_parquet_results<W: Write + Send + 'static>(
    g4s: &[G4],
    writer: W,
) -> Result<(), ExportError> {
    write_parquet_from_results(g4s, writer)
}

fn write_parquet_from_results<W: Write + Send + 'static>(
    g4s: &[G4],
    writer: W,
) -> Result<(), ExportError> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("start", DataType::UInt64, false),
        Field::new("end", DataType::UInt64, false),
        Field::new("length", DataType::UInt64, false),
        Field::new("tetrads", DataType::UInt64, false),
        Field::new("y1", DataType::Int32, false),
        Field::new("y2", DataType::Int32, false),
        Field::new("y3", DataType::Int32, false),
        Field::new("gscore", DataType::Int32, false),
        Field::new("sequence", DataType::Utf8, false),
    ]));

    let starts: Vec<u64> = g4s.iter().map(|g| g.start as u64).collect();
    let ends: Vec<u64> = g4s.iter().map(|g| g.end as u64).collect();
    let lengths: Vec<u64> = g4s.iter().map(|g| g.length as u64).collect();
    let tetrads: Vec<u64> = g4s.iter().map(|g| g.tetrads as u64).collect();
    let y1s: Vec<i32> = g4s.iter().map(|g| g.y1).collect();
    let y2s: Vec<i32> = g4s.iter().map(|g| g.y2).collect();
    let y3s: Vec<i32> = g4s.iter().map(|g| g.y3).collect();
    let gscores: Vec<i32> = g4s.iter().map(|g| g.gscore).collect();
    let sequences: Vec<String> = g4s.iter().map(|g| g.sequence().to_string()).collect();

    let columns: Vec<ArrayRef> = vec![
        Arc::new(UInt64Array::from(starts)),
        Arc::new(UInt64Array::from(ends)),
        Arc::new(UInt64Array::from(lengths)),
        Arc::new(UInt64Array::from(tetrads)),
        Arc::new(Int32Array::from(y1s)),
        Arc::new(Int32Array::from(y2s)),
        Arc::new(Int32Array::from(y3s)),
        Arc::new(Int32Array::from(gscores)),
        Arc::new(StringArray::from(sequences)),
    ];

    let batch = RecordBatch::try_new(schema.clone(), columns)?;
    let mut arrow_writer = ArrowWriter::try_new(writer, schema, None)?;
    arrow_writer.write(&batch)?;
    arrow_writer.close()?;
    Ok(())
}

pub fn load_sequences_from_path(path: &Path, mode: InputMode) -> io::Result<Vec<ChromSequence>> {
    match mode {
        InputMode::Mmap => load_sequences_mmap(path),
        InputMode::Stream => load_sequences_stream(path),
    }
}

fn load_sequences_stream(path: &Path) -> io::Result<Vec<ChromSequence>> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut sequences = Vec::new();
    let mut current_name: Option<String> = None;
    let mut sequence: Vec<u8> = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
            current_name = Some(parse_chrom_name(&line, sequences.len() + 1));
            continue;
        }
        for byte in line.bytes() {
            if byte.is_ascii_whitespace() {
                continue;
            }
            sequence.push(byte.to_ascii_lowercase());
        }
    }
    finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
    if !sequence.is_empty() {
        sequences.push(ChromSequence {
            name: format!("chromosome_{}", sequences.len() + 1),
            sequence: Arc::new(std::mem::take(&mut sequence)),
        });
    }
    Ok(sequences)
}

fn load_sequences_mmap(path: &Path) -> io::Result<Vec<ChromSequence>> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let mut sequences = Vec::new();
    let mut sequence = Vec::with_capacity(mmap.len());
    let mut current_name: Option<String> = None;
    let mut at_line_start = true;
    let mut i = 0;
    while i < mmap.len() {
        let byte = mmap[i];
        if byte == b'\n' || byte == b'\r' {
            at_line_start = true;
            i += 1;
            continue;
        }
        if at_line_start && byte == b'>' {
            finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
            i += 1;
            let header_start = i;
            while i < mmap.len() && mmap[i] != b'\n' && mmap[i] != b'\r' {
                i += 1;
            }
            let header = &mmap[header_start..i];
            current_name = Some(parse_chrom_name_bytes(header, sequences.len() + 1));
            at_line_start = true;
            continue;
        }
        at_line_start = false;
        if byte.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        sequence.push(byte.to_ascii_lowercase());
        i += 1;
    }
    finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
    if !sequence.is_empty() {
        let fallback =
            current_name.unwrap_or_else(|| format!("chromosome_{}", sequences.len() + 1));
        sequences.push(ChromSequence {
            name: fallback,
            sequence: Arc::new(std::mem::take(&mut sequence)),
        });
    }
    Ok(sequences)
}

fn finalize_sequence(
    current_name: &mut Option<String>,
    sequence: &mut Vec<u8>,
    sequences: &mut Vec<ChromSequence>,
) {
    if let Some(name) = current_name.take()
        && !sequence.is_empty()
    {
        sequences.push(ChromSequence {
            name,
            sequence: Arc::new(std::mem::take(sequence)),
        });
    }
}

pub(super) fn parse_chrom_name(line: &str, index: usize) -> String {
    let header = line.trim_start_matches('>');
    header
        .split_whitespace()
        .next()
        .map(|s| s.to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| format!("chromosome_{index}"))
}

pub(super) fn parse_chrom_name_bytes(header: &[u8], index: usize) -> String {
    let header_str = std::str::from_utf8(header).unwrap_or("");
    parse_chrom_name(header_str, index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{collections::HashMap, env, fs, sync::Arc};

    type G4Signature = (usize, usize, usize, usize, i32, i32, i32, i32, String);

    fn arc_from_sequence(seq: &str) -> Arc<Vec<u8>> {
        Arc::new(seq.bytes().map(|b| b.to_ascii_lowercase()).collect())
    }

    fn g4_signatures(g4s: &[G4]) -> Vec<G4Signature> {
        let mut sigs: Vec<_> = g4s
            .iter()
            .map(|g| {
                (
                    g.start,
                    g.end,
                    g.length,
                    g.tetrads,
                    g.y1,
                    g.y2,
                    g.y3,
                    g.gscore,
                    g.sequence().to_string(),
                )
            })
            .collect();
        sigs.sort();
        sigs
    }

    #[test]
    fn finds_single_g4() {
        let sequence = "GGGGAGGGGAGGGGAGGGG";
        let results = find_owned_bytes(arc_from_sequence(sequence), 4, 17);
        assert_eq!(results.len(), 1);
        let g = &results[0];
        assert_eq!(g.start, 1);
        assert_eq!(g.tetrads, 4);
        assert_eq!(g.y1, 1);
        assert_eq!(g.y2, 1);
        assert_eq!(g.y3, 1);
        assert_eq!(g.sequence(), sequence);
    }

    #[test]
    fn empty_sequence_has_no_hits() {
        let results = find_owned_bytes(arc_from_sequence("ACACAC"), 4, 17);
        assert!(results.is_empty());
    }

    #[test]
    fn csv_output_includes_header_and_rows() {
        let sequence = "GGGGAGGGGAGGGGAGGGG";
        let results = find_owned_bytes(arc_from_sequence(sequence), 4, 17);
        let csv = render_csv_results(&results);
        assert!(csv.starts_with("start,end,length"));
        assert!(csv.contains("GGGGAGGGGAGGGGAGGGG"));
    }

    #[test]
    fn parquet_writer_emits_bytes() {
        let sequence = "GGGGAGGGGAGGGGAGGGG";
        let path = env::temp_dir().join("qgrs_parquet_test.parquet");
        let file = fs::File::create(&path).expect("temp parquet file");
        let results = find_owned_bytes(arc_from_sequence(sequence), 4, 17);
        write_parquet_results(&results, file).expect("parquet export");
        let metadata = fs::metadata(&path).expect("metadata");
        assert!(metadata.len() > 0);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn load_sequences_stream_mode_splits_chromosomes() {
        let path = env::temp_dir().join("qgrs_stream_input.fa");
        fs::write(&path, b">chr1 description\nGGGG\nAC\n>chr2\nTTTT\n").unwrap();
        let seqs = load_sequences_from_path(&path, InputMode::Stream).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "chr1");
        assert_eq!(seqs[0].as_uppercase_string(), "GGGGAC");
        assert_eq!(seqs[1].name, "chr2");
        assert_eq!(seqs[1].as_uppercase_string(), "TTTT");
        fs::remove_file(&path).unwrap();
    }

    #[test]
    fn load_sequences_mmap_mode_splits_chromosomes() {
        let path = env::temp_dir().join("qgrs_mmap_input.fa");
        fs::write(&path, b">chr1\r\nGGGG\r\nAC\r\n>chrX\r\nCCCC\r\n").unwrap();
        let seqs = load_sequences_from_path(&path, InputMode::Mmap).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "chr1");
        assert_eq!(seqs[0].as_uppercase_string(), "GGGGAC");
        assert_eq!(seqs[1].name, "chrX");
        assert_eq!(seqs[1].as_uppercase_string(), "CCCC");
        fs::remove_file(&path).unwrap();
    }

    #[test]
    fn stream_pipeline_matches_batch_results() {
        let path = env::temp_dir().join("qgrs_stream_pipeline.fa");
        let fasta = b">chr1 desc\nGGGGAGGGGTTTTGGGG\n>chr2\nACACGGGGACACGGGG\n";
        fs::write(&path, fasta).unwrap();
        let sequences = load_sequences_from_path(&path, InputMode::Stream).unwrap();
        let mut expected: HashMap<String, Vec<G4>> = HashMap::new();
        for chrom in &sequences {
            let hits = find_owned_bytes(chrom.sequence.clone(), 2, 17);
            expected.insert(chrom.name.clone(), hits);
        }
        let mut actual: HashMap<String, Vec<G4>> = HashMap::new();
        stream::process_fasta_stream(&path, 2, 17, |name, results| {
            actual.insert(name, results);
            Ok(())
        })
        .unwrap();
        assert_eq!(expected.len(), actual.len());
        for (name, expected_hits) in expected {
            let observed = actual.get(&name).expect("missing chromosome");
            assert_eq!(expected_hits.len(), observed.len());
            for (lhs, rhs) in expected_hits.iter().zip(observed.iter()) {
                assert_eq!(lhs.start, rhs.start);
                assert_eq!(lhs.end, rhs.end);
                assert_eq!(lhs.sequence(), rhs.sequence());
                assert_eq!(lhs.tetrads, rhs.tetrads);
                assert_eq!(lhs.gscore, rhs.gscore);
            }
        }
        fs::remove_file(&path).unwrap();
    }

    #[test]
    fn chunked_search_matches_internal_results() {
        let limits = ScanLimits::default();
        let chunk_size = chunk_size_for_limits(limits);
        assert!(chunk_size < 100);
        let mut sequence = String::new();
        sequence.push_str(&"A".repeat(chunk_size - 5));
        sequence.push_str("GGGGAGGGGAGGGGAGGGG");
        sequence.push_str(&"A".repeat(chunk_size / 2));
        sequence.push_str("GGGGAGGGGAGGGGAGGGG");
        sequence.push_str(&"T".repeat(10));

        let chunked = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 4, 17, limits);

        let starts: Vec<_> = chunked.iter().map(|g| g.start).collect();
        assert_eq!(starts.len(), 2);
        let expected_starts = vec![
            chunk_size - 5 + 1,
            (chunk_size - 5) + "GGGGAGGGGAGGGGAGGGG".len() + chunk_size / 2 + 1,
        ];
        assert_eq!(starts, expected_starts);
    }

    #[test]
    fn chunked_bytes_matches_full_scan_on_boundary() {
        let limits = ScanLimits::default();
        let chunk_size = chunk_size_for_limits(limits);
        let mut sequence = String::new();
        sequence.push_str(&"A".repeat(chunk_size - 5));
        sequence.push_str("GGGGAGGGGAGGGGAGGGG");
        sequence.push_str(&"T".repeat(32));

        let chunked = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 4, 17, limits);
        let reference = find_with_sequence(Arc::new(SequenceData::new(&sequence)), 4, 17, limits);

        assert_eq!(g4_signatures(&chunked), g4_signatures(&reference));
    }

    #[test]
    fn chunked_bytes_handles_adjacent_cross_boundary_families() {
        let limits = ScanLimits::default();
        let chunk_size = chunk_size_for_limits(limits);
        let mut sequence = String::new();
        sequence.push_str(&"C".repeat(chunk_size - 8));
        sequence.push_str("GGGGAGGGGAGGGGAGGGG");
        sequence.push_str("AA");
        sequence.push_str("GGGGTTGGGGTTGGGGTTGGGG");
        sequence.push_str(&"C".repeat(24));

        let chunked = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 4, 17, limits);
        let reference = find_with_sequence(Arc::new(SequenceData::new(&sequence)), 4, 17, limits);

        assert_eq!(g4_signatures(&chunked), g4_signatures(&reference));
    }

    fn load_big_sequence() -> String {
        let fasta = include_str!("../../big.txt");
        fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect::<String>()
    }

    #[test]
    fn big_sequence_internal_equals_chunked() {
        let sequence = load_big_sequence();
        let limits = ScanLimits::default();
        let chunked = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 2, 17, limits);
        let internal = find_with_sequence(Arc::new(SequenceData::new(&sequence)), 2, 17, limits);
        assert_eq!(g4_signatures(&chunked), g4_signatures(&internal));
    }
}
