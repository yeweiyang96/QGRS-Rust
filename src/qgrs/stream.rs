use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, Sender};

use rayon::spawn;

use super::{
    G4, ScanLimits, SearchResults, chunk_size_for_limits, compute_chunk_overlap,
    finalize_search_results, find_raw_bytes_no_chunking, parse_chrom_name, shift_g4,
};

// debug helpers removed

pub fn process_fasta_stream<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_fasta_stream_with_limits(
        path,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1 << 20, file);
    process_reader_with_limits(reader, min_tetrads, min_score, limits, &mut on_chromosome)
}

pub fn process_fasta_stream_with_limits_ex<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    collect_families: bool,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(SearchResults) -> io::Result<()>,
{
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1 << 20, file);
    process_reader_with_limits_ex(
        reader,
        min_tetrads,
        min_score,
        limits,
        collect_families,
        &mut on_chromosome,
    )
}

pub fn process_reader<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_reader_with_limits(
        reader,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        on_chromosome,
    )
}

pub fn process_reader_with_limits<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let mut adapter = |result: SearchResults| {
        on_chromosome(result.chrom_name, result.g4s)?;
        Ok(())
    };
    process_reader_with_limits_ex(reader, min_tetrads, min_score, limits, false, &mut adapter)
}

pub fn process_reader_with_limits_ex<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    collect_families: bool,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(SearchResults) -> io::Result<()>,
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
                let results = chrom.finish();
                on_chromosome(results)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                collect_families,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new(
                fallback,
                min_tetrads,
                min_score,
                limits,
                collect_families,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let results = chrom.finish();
        on_chromosome(results)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

struct StreamChromosome {
    name: String,
    scheduler: StreamChunkScheduler,
    collect_families: bool,
}

impl StreamChromosome {
    fn new(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        collect_families: bool,
    ) -> Self {
        Self {
            name,
            scheduler: StreamChunkScheduler::new(min_tetrads, min_score, limits),
            collect_families,
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.scheduler.push_byte(byte);
    }

    fn finish(self) -> SearchResults {
        let raw_hits = self.scheduler.finish();
        finalize_search_results(raw_hits, &self.name, self.collect_families)
    }
}

struct StreamChunkScheduler {
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    chunk_size: usize,
    overlap: usize,
    buffer: VecDeque<u8>,
    offset: usize,
    tx: Sender<Vec<G4>>,
    rx: Receiver<Vec<G4>>,
    inflight: usize,
}

impl StreamChunkScheduler {
    fn new(min_tetrads: usize, min_score: i32, limits: ScanLimits) -> Self {
        let (tx, rx) = mpsc::channel();
        let chunk_size = chunk_size_for_limits(limits);
        let overlap = compute_chunk_overlap(min_tetrads, limits);
        let capacity = chunk_size + overlap;
        Self {
            min_tetrads,
            min_score,
            limits,
            chunk_size,
            overlap,
            buffer: VecDeque::with_capacity(capacity),
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
        let (front, back) = self.buffer.as_slices();
        if take <= front.len() {
            chunk.extend_from_slice(&front[..take]);
        } else {
            chunk.extend_from_slice(front);
            let remaining = take - front.len();
            chunk.extend_from_slice(&back[..remaining]);
        }
        // Efficiently remove the primary_len elements from the front.
        self.buffer.drain(..primary_len);
        let offset = self.offset;
        let _cutoff = offset + primary_len;
        self.offset += primary_len;
        let min_tetrads = self.min_tetrads;
        let min_score = self.min_score;
        let limits = self.limits;
        let tx = self.tx.clone();
        self.inflight += 1;
        spawn(move || {
            // Use the no-chunking variant here: the scheduler already supplied
            // a window (primary + overlap) and we must not re-chunk it.
            let mut hits = find_raw_bytes_no_chunking(chunk, min_tetrads, min_score, limits);
            for g4 in &mut hits {
                shift_g4(g4, offset);
            }
            // worker-local dedup is disabled; send raw hits to consolidator
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
        combined
    }
}
