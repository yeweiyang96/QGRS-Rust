use std::cell::RefCell;
use std::collections::VecDeque;
use std::sync::{Arc, OnceLock};

use memchr::memchr2;

use crate::qgrs::data::{ScanLimits, SequenceData, SequenceSlice, is_g};

// Invariants for the raw-search layer:
// 1. All coordinates remain 0-based half-open internally. `G4::start` is adjusted
//    to 1-based when materializing results, so upstream logic must not re-shift.
// 2. Chunked scans provide windows shaped as (primary, primary + overlap). This
//    module must clamp every emission to `primary_end` so overlap regions do not
//    double-count. Stream workers follow the same rule and must never re-chunk
//    the window they receive from the scheduler.
// 3. Raw finders never perform deduplication. Callers are responsible for passing
//    concatenated hits to `consolidate_g4s` to preserve parity across mmap and
//    streaming paths.

thread_local! {
    static LOOP_BUFFER: RefCell<Vec<i32>> = RefCell::new(Vec::with_capacity(16));
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
    slice_start: usize,
    sequence_data: Arc<Vec<u8>>,
    slice_cache: OnceLock<SequenceSlice>,
    sequence_cache: OnceLock<String>,
}

impl G4 {
    fn from_candidate(candidate: &G4Candidate) -> Self {
        let length = candidate.length();
        let end = candidate.start + length;
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
            slice_start: candidate.start,
            sequence_data: candidate.seq.normalized.clone(),
            slice_cache: OnceLock::new(),
            sequence_cache: OnceLock::new(),
        }
    }

    pub fn sequence(&self) -> &str {
        self.sequence_cache
            .get_or_init(|| self.sequence_slice().to_uppercase_string())
    }

    pub(crate) fn sequence_slice(&self) -> SequenceSlice {
        self.slice_cache
            .get_or_init(|| {
                SequenceSlice::new(self.sequence_data.clone(), self.slice_start, self.length)
            })
            .clone()
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
            slice_start: self.slice_start,
            sequence_data: self.sequence_data.clone(),
            slice_cache: OnceLock::new(),
            sequence_cache: OnceLock::new(),
        }
    }
}

#[derive(Clone)]
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
            LOOP_BUFFER.with(|slot| {
                let mut ys = slot.borrow_mut();
                ys.clear();
                self.find_loop_lengths_from(&mut ys, cursor);
                for &y in ys.iter() {
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
            });
        }
        results
    }
}

pub(crate) struct GRunScanner<'a> {
    data: &'a [u8],
    cursor: usize,
    min_tetrads: usize,
    max_g_run: usize,
}

impl<'a> GRunScanner<'a> {
    pub(crate) fn new(data: &'a [u8], min_tetrads: usize, max_g_run: usize) -> Self {
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

pub(crate) fn maximum_length(num_tetrads: usize, limits: ScanLimits) -> usize {
    let base = if num_tetrads < 3 { 30 } else { 45 };
    base.min(limits.max_g4_length)
}

pub(crate) fn find_raw_bytes_no_chunking(
    mut sequence: Vec<u8>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    sequence.make_ascii_lowercase();
    let seq = Arc::new(SequenceData::from_bytes(Arc::new(sequence)));
    find_raw_with_sequence(seq, min_tetrads, min_score, limits)
}

pub(crate) fn find_raw_on_window_bytes(
    seq: Arc<SequenceData>,
    base_offset: usize,
    primary_end: usize,
    window_end: usize,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    // Called by chunked batch scans only. The chunk scheduler already shapes
    // windows as (primary, primary+overlap) and expects this function to avoid
    // emitting hits whose start ≥ primary_end so that overlap regions don't
    // double-count.
    let window = &seq.normalized[base_offset..window_end];
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
                let base_max_offset = capped_run_len - tetrads;
                // Ensure each raw hit is emitted by at most one window by
                // clamping offsets to the primary section of the chunk.
                let boundary_offset = primary_end.saturating_sub(run_start + 1);
                let allowed_offset = base_max_offset.min(boundary_offset);
                for offset in 0..=allowed_offset {
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
    raw_g4s.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));
    raw_g4s
}

pub(crate) fn find_raw_with_sequence(
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
    raw_g4s.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));
    raw_g4s
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
