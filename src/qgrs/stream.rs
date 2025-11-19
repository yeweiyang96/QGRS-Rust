use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, Sender};

use rayon::spawn;

use super::{
    G4,
    CHUNK_SIZE,
    consolidate_g4s,
    compute_chunk_overlap,
    is_g,
    maximum_length,
    parse_chrom_name,
    shift_g4,
};

const WINDOW_CAPACITY: usize = 512;
const MAX_TETRADS_ALLOWED: usize = 11;

pub fn process_fasta_stream<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1 << 20, file);
    process_reader(reader, min_tetrads, min_score, &mut on_chromosome)
}

pub fn process_reader<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results) = chrom.finish();
                on_chromosome(name, results)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(name, min_tetrads, min_score));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new(fallback, min_tetrads, min_score));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte);
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results) = chrom.finish();
        on_chromosome(name, results)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

struct StreamChromosome {
    name: String,
    scheduler: StreamChunkScheduler,
}

impl StreamChromosome {
    fn new(name: String, min_tetrads: usize, min_score: i32) -> Self {
        Self {
            name,
            scheduler: StreamChunkScheduler::new(min_tetrads, min_score),
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.scheduler.push_byte(byte);
    }

    fn finish(self) -> (String, Vec<G4>) {
        let results = self.scheduler.finish();
        (self.name, results)
    }
}

struct StreamChunkScheduler {
    min_tetrads: usize,
    min_score: i32,
    chunk_size: usize,
    overlap: usize,
    buffer: VecDeque<u8>,
    offset: usize,
    tx: Sender<Vec<G4>>,
    rx: Receiver<Vec<G4>>,
    inflight: usize,
}

impl StreamChunkScheduler {
    fn new(min_tetrads: usize, min_score: i32) -> Self {
        let (tx, rx) = mpsc::channel();
        Self {
            min_tetrads,
            min_score,
            chunk_size: CHUNK_SIZE,
            overlap: compute_chunk_overlap(min_tetrads),
            buffer: VecDeque::with_capacity(CHUNK_SIZE),
            offset: 0,
            tx,
            rx,
            inflight: 0,
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.buffer.push_back(byte);
        self.flush_ready_chunks(false);
    }

    fn flush_ready_chunks(&mut self, finishing: bool) {
        let threshold = self.chunk_size + self.overlap;
        while self.buffer.len() >= threshold {
            self.dispatch_chunk(false, threshold);
        }
        if finishing && !self.buffer.is_empty() {
            self.dispatch_chunk(true, self.buffer.len());
        }
    }

    fn dispatch_chunk(&mut self, is_last: bool, window_len: usize) {
        if self.buffer.is_empty() {
            return;
        }
        let primary_len = if is_last {
            self.buffer.len()
        } else {
            self.chunk_size.min(self.buffer.len())
        };
        if primary_len == 0 {
            return;
        }
        let take = window_len.min(self.buffer.len());
        let mut chunk = Vec::with_capacity(take);
        for &byte in self.buffer.iter().take(take) {
            chunk.push(byte);
        }
        for _ in 0..primary_len {
            self.buffer.pop_front();
        }
        let offset = self.offset;
        let cutoff = offset + primary_len;
        self.offset += primary_len;
        let min_tetrads = self.min_tetrads;
        let min_score = self.min_score;
        let tx = self.tx.clone();
        self.inflight += 1;
        spawn(move || {
            let mut detector = StreamDetector::new(min_tetrads, min_score);
            for byte in chunk {
                detector.push_byte(byte);
            }
            let mut hits = detector.finish();
            for g4 in &mut hits {
                shift_g4(g4, offset);
            }
            if !is_last {
                hits.retain(|g4| g4.start < cutoff);
            }
            let _ = tx.send(hits);
        });
    }

    fn finish(mut self) -> Vec<G4> {
        self.flush_ready_chunks(true);
        drop(self.tx);
        let mut combined = Vec::new();
        for _ in 0..self.inflight {
            if let Ok(mut chunk) = self.rx.recv() {
                combined.append(&mut chunk);
            }
        }
        consolidate_g4s(combined)
    }
}

struct StreamDetector {
    min_tetrads: usize,
    min_score: i32,
    position: usize,
    window: SlidingWindow,
    run: Option<GRunTracker>,
    candidates: Vec<StreamCandidate>,
    raw_g4s: Vec<G4>,
}

impl StreamDetector {
    fn new(min_tetrads: usize, min_score: i32) -> Self {
        Self {
            min_tetrads,
            min_score,
            position: 0,
            window: SlidingWindow::new(WINDOW_CAPACITY),
            run: None,
            candidates: Vec::new(),
            raw_g4s: Vec::new(),
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.window.push(byte);
        let idx = self.position;
        self.position += 1;
        if is_g(byte) {
            self.extend_run(idx);
        } else {
            self.run = None;
        }
        self.advance_candidates();
        let retain_from = self.position.saturating_sub(self.window.capacity());
        self.window.trim(retain_from);
    }

    fn extend_run(&mut self, idx: usize) {
        let mut pending = Vec::new();
        if let Some(run) = self.run.as_mut() {
            run.push(&mut |start, tetrads| pending.push((start, tetrads)));
        } else {
            let mut tracker = GRunTracker::new(idx, self.min_tetrads);
            tracker.push(&mut |start, tetrads| pending.push((start, tetrads)));
            self.run = Some(tracker);
        }
        for (start, tetrads) in pending {
            self.spawn_candidate(start, tetrads);
        }
    }

    fn spawn_candidate(&mut self, start: usize, tetrads: usize) {
        if !tetrad_supported(tetrads) {
            return;
        }
        self.candidates.push(StreamCandidate::new(start, tetrads));
    }

    fn advance_candidates(&mut self) {
        let available = self.position;
        let window = &self.window;
        let min_score = self.min_score;
        let raw = &mut self.raw_g4s;
        let mut i = 0;
        while i < self.candidates.len() {
            let mut new_children = Vec::new();
            let state = {
                let cand = &mut self.candidates[i];
                cand.advance(available, window, min_score, raw, &mut new_children)
            };
            match state {
                CandidateState::Pending => i += 1,
                _ => {
                    self.candidates.swap_remove(i);
                }
            }
            self.candidates.extend(new_children);
        }
    }

    fn finish(self) -> Vec<G4> {
        consolidate_g4s(self.raw_g4s)
    }
}

struct SlidingWindow {
    original: VecDeque<u8>,
    normalized: VecDeque<u8>,
    head_index: usize,
    capacity: usize,
}

impl SlidingWindow {
    fn new(capacity: usize) -> Self {
        Self {
            original: VecDeque::with_capacity(capacity),
            normalized: VecDeque::with_capacity(capacity),
            head_index: 0,
            capacity,
        }
    }

    fn capacity(&self) -> usize {
        self.capacity
    }

    fn push(&mut self, byte: u8) {
        self.original.push_back(byte);
        self.normalized.push_back(byte.to_ascii_lowercase());
    }

    fn trim(&mut self, retain_from: usize) {
        while self.head_index < retain_from {
            if self.original.pop_front().is_none() {
                break;
            }
            self.normalized.pop_front();
            self.head_index += 1;
        }
    }

    fn is_all_g(&self, start: usize, len: usize) -> Option<bool> {
        let offset = start.checked_sub(self.head_index)?;
        if offset + len > self.normalized.len() {
            return None;
        }
        for i in 0..len {
            if *self.normalized.get(offset + i)? != b'g' {
                return Some(false);
            }
        }
        Some(true)
    }

    fn extract_original(&self, start: usize, len: usize) -> Option<String> {
        let offset = start.checked_sub(self.head_index)?;
        if offset + len > self.original.len() {
            return None;
        }
        let mut buf = Vec::with_capacity(len);
        for i in 0..len {
            buf.push(*self.original.get(offset + i)?);
        }
        String::from_utf8(buf).ok()
    }
}

struct GRunTracker {
    start: usize,
    length: usize,
    tetrad_lengths: Vec<usize>,
    emitted: Vec<usize>,
}

impl GRunTracker {
    fn new(start: usize, min_tetrads: usize) -> Self {
        let tetrad_lengths = supported_tetrads(min_tetrads);
        let emitted = vec![0; tetrad_lengths.len()];
        Self {
            start,
            length: 0,
            tetrad_lengths,
            emitted,
        }
    }

    fn push<F>(&mut self, emit: &mut F)
    where
        F: FnMut(usize, usize),
    {
        self.length += 1;
        self.emit_available(emit);
    }

    fn emit_available<F>(&mut self, emit: &mut F)
    where
        F: FnMut(usize, usize),
    {
        for (idx, &tetrad) in self.tetrad_lengths.iter().enumerate() {
            if self.length < tetrad {
                break;
            }
            let total = self.length - tetrad + 1;
            while self.emitted[idx] < total {
                let start = self.start + self.emitted[idx];
                emit(start, tetrad);
                self.emitted[idx] += 1;
            }
        }
    }
}

fn supported_tetrads(min_tetrads: usize) -> Vec<usize> {
    let mut lengths = Vec::new();
    for tetrad in min_tetrads..=MAX_TETRADS_ALLOWED {
        if tetrad_supported(tetrad) {
            lengths.push(tetrad);
        }
    }
    lengths
}

fn tetrad_supported(num_tetrads: usize) -> bool {
    4 * num_tetrads <= maximum_length(num_tetrads)
}

#[derive(Clone)]
struct StreamCandidate {
    start: usize,
    num_tetrads: usize,
    y: [i32; 3],
    stage: usize,
    cursor: usize,
    search_pos: usize,
    max_length: usize,
}

impl StreamCandidate {
    fn new(start: usize, num_tetrads: usize) -> Self {
        let max_length = maximum_length(num_tetrads);
        let cursor = start + num_tetrads;
        Self {
            start,
            num_tetrads,
            y: [-1; 3],
            stage: 0,
            cursor,
            search_pos: cursor,
            max_length,
        }
    }

    fn advance(
        &mut self,
        available_pos: usize,
        window: &SlidingWindow,
        min_score: i32,
        results: &mut Vec<G4>,
        expansions: &mut Vec<StreamCandidate>,
    ) -> CandidateState {
        if available_pos < self.num_tetrads {
            return CandidateState::Pending;
        }
        let max_scan_start = available_pos - self.num_tetrads;
        let mut position = self.search_pos;
        while position <= max_scan_start {
            if position >= self.start + self.max_length {
                self.search_pos = position;
                return CandidateState::Consume;
            }
            match window.is_all_g(position, self.num_tetrads) {
                Some(true) => {
                    let y = (position - self.cursor) as i32;
                    if y < self.min_acceptable_loop_length()
                        || (position - self.start + self.num_tetrads - 1) >= self.max_length
                    {
                        self.search_pos = position;
                        return CandidateState::Consume;
                    }
                    if let Some(next) = self.spawn_with_loop(y) {
                        if next.partial_length() <= next.max_length as i32 {
                            if next.complete() {
                                if next.viable(min_score) {
                                    if let Some(seq) =
                                        window.extract_original(next.start, next.length())
                                    {
                                        results.push(next.into_g4(seq));
                                    }
                                }
                            } else {
                                expansions.push(next);
                            }
                        }
                    }
                }
                Some(false) => {}
                None => {
                    self.search_pos = position;
                    return CandidateState::Consume;
                }
            }
            position += 1;
        }
        self.search_pos = position;
        CandidateState::Pending
    }

    fn spawn_with_loop(&self, value: i32) -> Option<StreamCandidate> {
        if self.stage >= self.y.len() {
            return None;
        }
        let mut next = self.clone();
        next.y[next.stage] = value;
        next.stage += 1;
        next.cursor = match next.stage {
            0 => next.start + next.num_tetrads,
            1 => next.t2() + next.num_tetrads,
            2 => next.t3() + next.num_tetrads,
            _ => next.t4(),
        };
        next.search_pos = next.cursor;
        Some(next)
    }

    fn length(&self) -> usize {
        (4 * self.num_tetrads)
            + self.y[0].max(0) as usize
            + self.y[1].max(0) as usize
            + self.y[2].max(0) as usize
    }

    fn t1(&self) -> usize {
        self.start
    }

    fn t2(&self) -> usize {
        self.t1() + self.num_tetrads + self.y[0].max(0) as usize
    }

    fn t3(&self) -> usize {
        self.t2() + self.num_tetrads + self.y[1].max(0) as usize
    }

    fn t4(&self) -> usize {
        self.t3() + self.num_tetrads + self.y[2].max(0) as usize
    }

    fn complete(&self) -> bool {
        self.stage >= 3
    }

    fn min_acceptable_loop_length(&self) -> i32 {
        if self.y.iter().any(|&y| y == 0) { 1 } else { 0 }
    }

    fn partial_length(&self) -> i32 {
        let mut length = (self.num_tetrads * 4) as i32;
        if self.y[0] >= 0 && self.y[1] < 0 {
            length += if self.y[0] == 0 { 2 } else { 1 };
        } else if self.y[1] >= 0 && self.y[2] < 0 {
            if self.y[0] == 0 || self.y[1] == 0 {
                length += 1;
            }
        }
        for &y in &self.y {
            if y > 0 {
                length += y;
            }
        }
        length
    }

    fn score(&self) -> i32 {
        let gavg = (f64::from((self.y[0] - self.y[1]).abs())
            + f64::from((self.y[1] - self.y[2]).abs())
            + f64::from((self.y[0] - self.y[2]).abs()))
            / 3.0;
        let gmax = (self.max_length as i32 - (self.num_tetrads as i32 * 4 + 1)) as f64;
        let bonus = gmax * ((self.num_tetrads as i32 - 2) as f64);
        (gmax - gavg + bonus).floor() as i32
    }

    fn viable(&self, min_score: i32) -> bool {
        if self.score() < min_score {
            return false;
        }
        if self.length() > self.max_length {
            return false;
        }
        let mut zero_loops = 0;
        for &y in &self.y {
            if y < 1 {
                zero_loops += 1;
            }
        }
        zero_loops < 2
    }

    fn into_g4(&self, sequence: String) -> G4 {
        let length = self.length();
        let end = self.start + length;
        G4 {
            start: self.start,
            end,
            tetrad1: self.t1(),
            tetrad2: self.t2(),
            tetrad3: self.t3(),
            tetrad4: self.t4(),
            y1: self.y[0],
            y2: self.y[1],
            y3: self.y[2],
            tetrads: self.num_tetrads,
            length,
            gscore: self.score(),
            sequence,
        }
    }
}

enum CandidateState {
    Pending,
    Consume,
}
