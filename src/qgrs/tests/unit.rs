use std::env;
use std::fs;

use crate::qgrs::{
    InputMode, ScanLimits, SequenceTopology, consolidate_g4s, consolidate_g4s_with_topology,
    find_owned_bytes, find_owned_bytes_with_topology, load_sequences_from_path, render_csv_results,
    render_csv_results_with_projection, render_family_ranges_csv_with_projection,
    write_parquet_results,
};

use super::helpers::arc_from_sequence;

#[test]
fn finds_single_g4() {
    let sequence = "GGGGAGGGGAGGGGAGGGG";
    let raw = find_owned_bytes(arc_from_sequence(sequence), 4, 17);
    let (results, _ranges) = consolidate_g4s(raw);
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
    let raw = find_owned_bytes(arc_from_sequence("ACACAC"), 4, 17);
    let (results, _ranges) = consolidate_g4s(raw);
    assert!(results.is_empty());
}

#[test]
fn csv_output_includes_header_and_rows() {
    let sequence = "GGGGAGGGGAGGGGAGGGG";
    let raw = find_owned_bytes(arc_from_sequence(sequence), 4, 17);
    let (results, _ranges) = consolidate_g4s(raw);
    let csv = render_csv_results(&results);
    assert!(csv.starts_with("start,end,length"));
    assert!(csv.contains("GGGGAGGGGAGGGGAGGGG"));
}

#[test]
fn parquet_writer_emits_bytes() {
    let sequence = "GGGGAGGGGAGGGGAGGGG";
    let path = env::temp_dir().join("qgrs_parquet_test.parquet");
    let file = fs::File::create(&path).expect("temp parquet file");
    let raw = find_owned_bytes(arc_from_sequence(sequence), 4, 17);
    let (results, _ranges) = consolidate_g4s(raw);
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
    assert_eq!(seqs[0].name(), "chr1");
    assert_eq!(seqs[0].as_uppercase_string(), "GGGGAC");
    assert_eq!(seqs[1].name(), "chr2");
    assert_eq!(seqs[1].as_uppercase_string(), "TTTT");
    fs::remove_file(&path).unwrap();
}

#[test]
fn load_sequences_mmap_mode_splits_chromosomes() {
    let path = env::temp_dir().join("qgrs_mmap_input.fa");
    fs::write(&path, b">chr1\r\nGGGG\r\nAC\r\n>chrX\r\nCCCC\r\n").unwrap();
    let seqs = load_sequences_from_path(&path, InputMode::Mmap).unwrap();
    assert_eq!(seqs.len(), 2);
    assert_eq!(seqs[0].name(), "chr1");
    assert_eq!(seqs[0].as_uppercase_string(), "GGGGAC");
    assert_eq!(seqs[1].name(), "chrX");
    assert_eq!(seqs[1].as_uppercase_string(), "CCCC");
    fs::remove_file(&path).unwrap();
}

#[test]
fn circular_mode_finds_wraparound_hit_when_linear_does_not() {
    let sequence = "GAGGGGAGGGGAGGGGGGG";
    let arc = arc_from_sequence(sequence);
    let limits = ScanLimits::default();

    let linear_raw =
        find_owned_bytes_with_topology(arc.clone(), 4, 17, limits, SequenceTopology::Linear);
    let (linear_hits, _ranges) =
        consolidate_g4s_with_topology(linear_raw, SequenceTopology::Linear, sequence.len());
    assert!(linear_hits.is_empty());

    let circular_raw =
        find_owned_bytes_with_topology(arc, 4, 17, limits, SequenceTopology::Circular);
    let (circular_hits, family_ranges) =
        consolidate_g4s_with_topology(circular_raw, SequenceTopology::Circular, sequence.len());
    assert_eq!(circular_hits.len(), 1);
    assert_eq!(family_ranges.len(), 1);

    let hit = &circular_hits[0];
    assert!(hit.start > 1);
    assert!(hit.end > sequence.len());
    assert_eq!(hit.sequence(), "GGGGAGGGGAGGGGAGGGG");
}

#[test]
fn circular_consolidation_merges_wraparound_family() {
    let sequence = "GAGGGGAGGGGAGGGGGGG";
    let limits = ScanLimits::default();
    let raw = find_owned_bytes_with_topology(
        arc_from_sequence(sequence),
        4,
        17,
        limits,
        SequenceTopology::Circular,
    );
    assert!(raw.len() > 1);
    let wrap_count = raw.iter().filter(|g4| g4.end > sequence.len()).count();
    assert!(wrap_count >= 2);

    let (hits, ranges) =
        consolidate_g4s_with_topology(raw, SequenceTopology::Circular, sequence.len());
    assert_eq!(hits.len(), 1);
    assert_eq!(ranges.len(), 1);
    assert!(ranges[0].1 > sequence.len());
}

#[test]
fn circular_projected_csv_maps_end_back_to_ring() {
    let sequence = "GAGGGGAGGGGAGGGGGGG";
    let limits = ScanLimits::default();
    let raw = find_owned_bytes_with_topology(
        arc_from_sequence(sequence),
        4,
        17,
        limits,
        SequenceTopology::Circular,
    );
    let (hits, ranges) =
        consolidate_g4s_with_topology(raw, SequenceTopology::Circular, sequence.len());

    let csv = render_csv_results_with_projection(&hits, SequenceTopology::Circular, sequence.len());
    assert!(csv.contains("\n17,16,19,4,1,1,1,84,GGGGAGGGGAGGGGAGGGG\n"));

    let family_csv = render_family_ranges_csv_with_projection(
        &ranges,
        SequenceTopology::Circular,
        sequence.len(),
    );
    let family_line = family_csv.lines().nth(1).expect("family row");
    let mut cols = family_line.split(',');
    assert_eq!(cols.next(), Some("1"));
    let start = cols.next().unwrap().parse::<usize>().unwrap();
    let end = cols.next().unwrap().parse::<usize>().unwrap();
    assert!(start <= sequence.len());
    assert!(end <= sequence.len());
    assert!(end < start);
}
