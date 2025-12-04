use std::sync::Arc;

use rayon::prelude::*;

use crate::qgrs::data::{ScanLimits, SequenceData};
use crate::qgrs::search::{G4, find_raw_on_window_bytes, find_raw_with_sequence};

const WINDOW_MIN_BP: usize = 32;
const WINDOW_MAX_BP: usize = 64;
const WINDOW_PADDING_BP: usize = 27;

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

        return merged_raw;
    }
    let seq = Arc::new(SequenceData::from_bytes(sequence));
    find_with_sequence(seq, min_tetrads, min_score, limits)
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
