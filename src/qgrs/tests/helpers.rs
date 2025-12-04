use std::sync::Arc;

use crate::qgrs::data::SequenceData;
use crate::qgrs::{G4, ScanLimits, consolidate_g4s, find_with_sequence};

pub(super) type G4Signature = (usize, usize, usize, usize, i32, i32, i32, i32, String);

pub(super) fn arc_from_sequence(seq: &str) -> Arc<Vec<u8>> {
    Arc::new(seq.bytes().map(|b| b.to_ascii_lowercase()).collect())
}

pub(super) fn g4_signatures(g4s: &[G4]) -> Vec<G4Signature> {
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

pub(super) fn load_big_sequence() -> String {
    let fasta = include_str!("../../../big.txt");
    fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>()
}

pub(super) fn run_internal_scan(
    sequence: &str,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
) -> Vec<G4> {
    let seq = Arc::new(SequenceData::new(sequence));
    let raw = find_with_sequence(seq, min_tetrads, min_score, limits);
    consolidate_g4s(raw)
}
