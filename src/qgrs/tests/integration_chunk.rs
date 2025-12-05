use crate::qgrs::{
    ScanLimits, chunk_size_for_limits, consolidate_g4s, find_owned_bytes_with_limits,
};

use super::helpers::{arc_from_sequence, g4_signatures, load_big_sequence, run_internal_scan};

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

    let chunked_raw = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 4, 17, limits);
    let (chunked, _ranges) = consolidate_g4s(chunked_raw);

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

    let chunked_raw = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 4, 17, limits);
    let (chunked, _ranges) = consolidate_g4s(chunked_raw);
    let reference = run_internal_scan(&sequence, 4, 17, limits);

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

    let chunked_raw = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 4, 17, limits);
    let (chunked, _ranges) = consolidate_g4s(chunked_raw);
    let reference = run_internal_scan(&sequence, 4, 17, limits);

    assert_eq!(g4_signatures(&chunked), g4_signatures(&reference));
}

#[test]
fn big_sequence_internal_equals_chunked() {
    let sequence = load_big_sequence();
    let limits = ScanLimits::default();
    let chunked_raw = find_owned_bytes_with_limits(arc_from_sequence(&sequence), 2, 17, limits);
    let (chunked, _ranges) = consolidate_g4s(chunked_raw);
    let internal = run_internal_scan(&sequence, 2, 17, limits);
    assert_eq!(g4_signatures(&chunked), g4_signatures(&internal));
}
