use std::env;
use std::fs;

use crate::qgrs::{
    InputMode, consolidate_g4s, find_owned_bytes, load_sequences_from_path, render_csv_results,
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
