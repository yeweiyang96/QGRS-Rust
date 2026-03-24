use std::collections::VecDeque;
use std::io::{self, BufRead};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, Sender};

use rayon::spawn;

use super::{
    G4, ScanLimits, SequenceTopology, chunk_size_for_limits, compute_chunk_overlap,
    consolidate_g4s_with_topology, find_raw_bytes_no_chunking, input::open_input_reader,
    parse_chrom_name, retain_circular_raw_hits, shift_g4,
};

pub struct StreamChromosomeResults {
    pub hits: Vec<G4>,
    pub family_ranges: Vec<(usize, usize)>,
    pub raw_hits: Option<Vec<G4>>,
}

pub fn process_fasta_stream<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_fasta_stream_with_limits_topology(
        path,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_overlap<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_fasta_stream_with_limits_overlap_topology(
        path,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_fasta_stream_with_limits_topology(
        path,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology_and_len<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology_and_len(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology_and_sequence<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>, Vec<u8>) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology_and_sequence(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_fasta_stream_with_limits_overlap_topology(
        path,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology_and_len<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology_and_len(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology_and_sequence<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults, Vec<u8>) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology_and_sequence(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
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
    process_reader_with_limits_topology(
        reader,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_overlap<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_reader_with_limits_overlap_topology(
        reader,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
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
    process_reader_with_limits_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_limits_topology<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
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
                topology,
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
                topology,
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
        let (name, results) = chrom.finish();
        on_chromosome(name, results)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_topology_and_len<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>, usize) -> io::Result<()>,
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
                let (name, results, sequence_len) = chrom.finish_with_sequence_len();
                on_chromosome(name, results, sequence_len)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
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
                topology,
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
        let (name, results, sequence_len) = chrom.finish_with_sequence_len();
        on_chromosome(name, results, sequence_len)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_topology_and_sequence<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>, Vec<u8>) -> io::Result<()>,
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
                let (name, results, sequence) = chrom.finish_with_sequence();
                on_chromosome(name, results, sequence)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
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
        let (name, results, sequence) = chrom.finish_with_sequence();
        on_chromosome(name, results, sequence)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

pub fn process_reader_with_limits_overlap<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_reader_with_limits_overlap_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_limits_overlap_topology<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
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
                let (name, results) = chrom.finish_with_overlap();
                on_chromosome(name, results)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
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
                topology,
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
        let (name, results) = chrom.finish_with_overlap();
        on_chromosome(name, results)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_overlap_topology_and_len<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults, usize) -> io::Result<()>,
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
                let (name, results, sequence_len) = chrom.finish_with_overlap_and_sequence_len();
                on_chromosome(name, results, sequence_len)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
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
                topology,
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
        let (name, results, sequence_len) = chrom.finish_with_overlap_and_sequence_len();
        on_chromosome(name, results, sequence_len)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_overlap_topology_and_sequence<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults, Vec<u8>) -> io::Result<()>,
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
                let (name, results, sequence) = chrom.finish_with_overlap_and_sequence();
                on_chromosome(name, results, sequence)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
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
        let (name, results, sequence) = chrom.finish_with_overlap_and_sequence();
        on_chromosome(name, results, sequence)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

struct StreamChromosome {
    name: String,
    scheduler: StreamChunkScheduler,
    captured_sequence: Option<Vec<u8>>,
}

impl StreamChromosome {
    fn new(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
    ) -> Self {
        Self::new_with_sequence_capture(name, min_tetrads, min_score, limits, topology, false)
    }

    fn new_with_sequence_capture(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        capture_sequence: bool,
    ) -> Self {
        Self {
            name,
            scheduler: StreamChunkScheduler::new(min_tetrads, min_score, limits, topology),
            captured_sequence: capture_sequence.then(Vec::new),
        }
    }

    fn push_byte(&mut self, byte: u8) {
        if let Some(sequence) = self.captured_sequence.as_mut() {
            sequence.push(byte);
        }
        self.scheduler.push_byte(byte);
    }

    fn finish(self) -> (String, Vec<G4>) {
        let results = self.scheduler.finish();
        (self.name, results)
    }

    fn finish_with_sequence_len(self) -> (String, Vec<G4>, usize) {
        let sequence_len = self.scheduler.sequence_len();
        let results = self.scheduler.finish();
        (self.name, results, sequence_len)
    }

    fn finish_with_sequence(self) -> (String, Vec<G4>, Vec<u8>) {
        let sequence = self.captured_sequence.unwrap_or_default();
        let results = self.scheduler.finish();
        (self.name, results, sequence)
    }

    fn finish_with_overlap(self) -> (String, StreamChromosomeResults) {
        let (hits, ranges, raw_hits) = self.scheduler.finish_with_overlap();
        (
            self.name,
            StreamChromosomeResults {
                hits,
                family_ranges: ranges,
                raw_hits: Some(raw_hits),
            },
        )
    }

    fn finish_with_overlap_and_sequence_len(self) -> (String, StreamChromosomeResults, usize) {
        let sequence_len = self.scheduler.sequence_len();
        let (name, results) = self.finish_with_overlap();
        (name, results, sequence_len)
    }

    fn finish_with_overlap_and_sequence(self) -> (String, StreamChromosomeResults, Vec<u8>) {
        let sequence = self.captured_sequence.unwrap_or_default();
        let (hits, ranges, raw_hits) = self.scheduler.finish_with_overlap();
        (
            self.name,
            StreamChromosomeResults {
                hits,
                family_ranges: ranges,
                raw_hits: Some(raw_hits),
            },
            sequence,
        )
    }
}

struct StreamChunkScheduler {
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    chunk_size: usize,
    overlap: usize,
    buffer: VecDeque<u8>,
    offset: usize,
    sequence_len: usize,
    circular_boundary_bp: usize,
    circular_head: VecDeque<u8>,
    circular_tail: VecDeque<u8>,
    tx: Sender<Vec<G4>>,
    rx: Receiver<Vec<G4>>,
    inflight: usize,
}

type FinishParts = (Vec<G4>, Vec<(usize, usize)>, Option<Vec<G4>>);

impl StreamChunkScheduler {
    fn new(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
    ) -> Self {
        let (tx, rx) = mpsc::channel();
        let chunk_size = chunk_size_for_limits(limits);
        let overlap = compute_chunk_overlap(min_tetrads, limits);
        let capacity = chunk_size + overlap;
        let circular_boundary_bp = if topology.is_circular() {
            limits.max_g4_length.saturating_sub(1)
        } else {
            0
        };
        Self {
            min_tetrads,
            min_score,
            limits,
            topology,
            chunk_size,
            overlap,
            buffer: VecDeque::with_capacity(capacity),
            offset: 0,
            sequence_len: 0,
            circular_boundary_bp,
            circular_head: VecDeque::with_capacity(circular_boundary_bp),
            circular_tail: VecDeque::with_capacity(circular_boundary_bp),
            tx,
            rx,
            inflight: 0,
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.sequence_len += 1;
        if self.circular_boundary_bp > 0 {
            if self.circular_head.len() < self.circular_boundary_bp {
                self.circular_head.push_back(byte);
            }
            self.circular_tail.push_back(byte);
            if self.circular_tail.len() > self.circular_boundary_bp {
                self.circular_tail.pop_front();
            }
        }
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

    fn finish(self) -> Vec<G4> {
        let (hits, _, _) = self.finish_internal(false);
        hits
    }

    fn finish_with_overlap(self) -> (Vec<G4>, Vec<(usize, usize)>, Vec<G4>) {
        let (hits, ranges, raw) = self.finish_internal(true);
        (
            hits,
            ranges,
            raw.expect("raw hits must be captured when capture_raw is true"),
        )
    }

    fn finish_internal(mut self, capture_raw: bool) -> FinishParts {
        self.flush_ready_chunks(true);
        let mut combined = Vec::new();
        for _ in 0..self.inflight {
            if let Ok(mut chunk) = self.rx.recv() {
                combined.append(&mut chunk);
            }
        }
        if self.topology.is_circular() {
            self.append_wraparound_hits(&mut combined);
            retain_circular_raw_hits(&mut combined, self.sequence_len);
        } else {
            combined.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));
        }
        let raw_hits = if capture_raw {
            Some(combined.clone())
        } else {
            None
        };
        let (hits, ranges) =
            consolidate_g4s_with_topology(combined, self.topology, self.sequence_len);
        (hits, ranges, raw_hits)
    }

    fn sequence_len(&self) -> usize {
        self.sequence_len
    }

    fn append_wraparound_hits(&self, combined: &mut Vec<G4>) {
        if self.sequence_len == 0
            || self.circular_boundary_bp == 0
            || self.circular_head.is_empty()
            || self.circular_tail.is_empty()
        {
            return;
        }
        let mut boundary = Vec::with_capacity(self.circular_tail.len() + self.circular_head.len());
        boundary.extend(self.circular_tail.iter().copied());
        boundary.extend(self.circular_head.iter().copied());
        let mut hits =
            find_raw_bytes_no_chunking(boundary, self.min_tetrads, self.min_score, self.limits);
        let offset = self.sequence_len.saturating_sub(self.circular_tail.len());
        for g4 in &mut hits {
            shift_g4(g4, offset);
        }
        hits.retain(|g4| g4.end > self.sequence_len);
        combined.extend(hits);
    }
}
