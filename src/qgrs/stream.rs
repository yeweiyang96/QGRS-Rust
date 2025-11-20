use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, Sender};

use rayon::spawn;

use crate::ProgressReporterHandle;

use super::{
    G4, ProgressScope, ScanLimits, chunk_size_for_limits, compute_chunk_overlap, consolidate_g4s,
    find_owned_with_limits, find_owned_with_limits_with_scope, parse_chrom_name, shift_g4,
};

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

pub fn process_fasta_stream_with_progress<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    reporter: ProgressReporterHandle,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_fasta_stream_with_limits_with_progress(
        path,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        reporter,
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
    process_reader_with_limits_internal(
        reader,
        min_tetrads,
        min_score,
        limits,
        None,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_with_progress<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    reporter: ProgressReporterHandle,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1 << 20, file);
    process_reader_with_limits_internal(
        reader,
        min_tetrads,
        min_score,
        limits,
        Some(reporter),
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
    process_reader_with_limits_internal(
        reader,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        None,
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
    process_reader_with_limits_internal(reader, min_tetrads, min_score, limits, None, on_chromosome)
}

pub fn process_reader_with_limits_with_progress<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    reporter: ProgressReporterHandle,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_reader_with_limits_internal(
        reader,
        min_tetrads,
        min_score,
        limits,
        Some(reporter),
        on_chromosome,
    )
}

fn process_reader_with_limits_internal<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    reporter: Option<ProgressReporterHandle>,
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
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                reporter.as_ref().map(|handle| handle.clone()),
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
                reporter.as_ref().map(|handle| handle.clone()),
            ));
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
    length: usize,
    progress: Option<ProgressScope>,
}

impl StreamChromosome {
    fn new(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        reporter: Option<ProgressReporterHandle>,
    ) -> Self {
        let scope = reporter.map(|handle| ProgressScope::new(&name, 1, handle));
        let scheduler = StreamChunkScheduler::new(min_tetrads, min_score, limits, scope.clone());
        Self {
            name,
            scheduler,
            length: 0,
            progress: scope,
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.length += 1;
        if let Some(scope) = &self.progress {
            scope.update_total(self.length);
            if self.length == 1 {
                scope.start();
            }
        }
        self.scheduler.push_byte(byte);
    }

    fn finish(self) -> (String, Vec<G4>) {
        if self.length == 0 {
            if let Some(scope) = &self.progress {
                scope.start();
            }
        }
        let results = self.scheduler.finish();
        if let Some(scope) = &self.progress {
            scope.finish();
        }
        (self.name, results)
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
    progress: Option<ProgressScope>,
}

impl StreamChunkScheduler {
    fn new(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        progress: Option<ProgressScope>,
    ) -> Self {
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
            progress,
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
        let limits = self.limits;
        let tx = self.tx.clone();
        let scope = self
            .progress
            .as_ref()
            .map(|progress| progress.with_offset(offset));
        self.inflight += 1;
        spawn(move || {
            let sequence: String = chunk.into_iter().map(|b| b as char).collect();
            let mut hits = if let Some(progress_scope) = scope {
                find_owned_with_limits_with_scope(
                    sequence,
                    min_tetrads,
                    min_score,
                    limits,
                    Some(progress_scope),
                )
            } else {
                find_owned_with_limits(sequence, min_tetrads, min_score, limits)
            };
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
