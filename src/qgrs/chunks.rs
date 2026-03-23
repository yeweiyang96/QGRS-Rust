use std::sync::Arc;

use rayon::prelude::*;

use crate::qgrs::data::{ScanLimits, SequenceData, SequenceTopology};
use crate::qgrs::search::{G4, find_raw_on_window_bytes, find_raw_with_sequence};

const WINDOW_MIN_BP: usize = 32;
const WINDOW_MAX_BP: usize = 64;
const WINDOW_PADDING_BP: usize = 27;

pub fn find_owned_bytes(sequence: Arc<Vec<u8>>, min_tetrads: usize, min_score: i32) -> Vec<G4> {
    find_owned_bytes_with_topology(
        sequence,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
    )
}

pub fn find_owned_bytes_with_limits(
    sequence: Arc<Vec<u8>>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    find_owned_bytes_with_topology(
        sequence,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
    )
}

pub fn find_owned_bytes_with_topology(
    sequence: Arc<Vec<u8>>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
) -> Vec<G4> {
    if topology.is_circular() {
        return find_owned_bytes_circular(sequence, min_tetrads, min_score, limits);
    }
    find_owned_bytes_linear(sequence, min_tetrads, min_score, limits)
}

fn find_owned_bytes_linear(
    sequence: Arc<Vec<u8>>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    let chunk_size = chunk_size_for_limits(limits);
    if sequence.len() > chunk_size {
        let len = sequence.len();
        let overlap = compute_chunk_overlap(min_tetrads, limits);
        let mut start = 0usize;
        let seq_data = Arc::new(SequenceData::from_bytes(sequence.clone()));
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
                let hits = find_raw_on_window_bytes(
                    seq_data.clone(),
                    offset,
                    primary_end,
                    window_end,
                    min_tetrads,
                    min_score,
                    limits,
                );
                hits.into_iter()
            })
            .collect();

        return merged_raw;
    }
    let seq = Arc::new(SequenceData::from_bytes(sequence));
    find_with_sequence(seq, min_tetrads, min_score, limits)
}

fn find_owned_bytes_circular(
    sequence: Arc<Vec<u8>>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    let sequence_len = sequence.len();
    if sequence_len == 0 {
        return Vec::new();
    }
    let prefix_len = circular_prefix_len(sequence_len, limits);
    let mut extended = Vec::with_capacity(sequence_len + prefix_len);
    extended.extend_from_slice(sequence.as_slice());
    if prefix_len > 0 {
        extended.extend_from_slice(&sequence[..prefix_len]);
    }
    let mut hits = find_owned_bytes_linear(Arc::new(extended), min_tetrads, min_score, limits);
    retain_circular_raw_hits(&mut hits, sequence_len);
    hits
}

fn circular_prefix_len(sequence_len: usize, limits: ScanLimits) -> usize {
    if sequence_len <= 1 {
        return 0;
    }
    limits
        .max_g4_length
        .saturating_sub(1)
        .min(sequence_len.saturating_sub(1))
}

pub(crate) fn retain_circular_raw_hits(raw_hits: &mut Vec<G4>, sequence_len: usize) {
    if sequence_len == 0 {
        raw_hits.clear();
        return;
    }
    raw_hits.retain(|g4| g4.start <= sequence_len && g4.length <= sequence_len);
    raw_hits.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));
}

pub(crate) fn chunk_size_for_limits(limits: ScanLimits) -> usize {
    let desired = limits.max_g4_length.saturating_add(WINDOW_PADDING_BP);
    desired.clamp(WINDOW_MIN_BP, WINDOW_MAX_BP)
}

pub(crate) fn compute_chunk_overlap(_min_tetrads: usize, limits: ScanLimits) -> usize {
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

pub(crate) fn find_with_sequence(
    seq: Arc<SequenceData>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    find_raw_with_sequence(seq, min_tetrads, min_score, limits)
}
