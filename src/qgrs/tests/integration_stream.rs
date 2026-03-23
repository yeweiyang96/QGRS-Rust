use std::collections::HashMap;
use std::fs;

use crate::qgrs::stream;
use crate::qgrs::{
    InputMode, ScanLimits, SequenceTopology, consolidate_g4s, consolidate_g4s_with_topology,
    find_owned_bytes, find_owned_bytes_with_topology,
};

#[test]
fn stream_pipeline_matches_batch_results() {
    let path = std::env::temp_dir().join("qgrs_stream_pipeline.fa");
    let fasta = b">chr1 desc\nGGGGAGGGGTTTTGGGG\n>chr2\nACACGGGGACACGGGG\n";
    fs::write(&path, fasta).unwrap();
    let sequences = crate::qgrs::load_sequences_from_path(&path, InputMode::Stream).unwrap();
    let mut expected: HashMap<String, Vec<_>> = HashMap::new();
    for chrom in &sequences {
        let raw = find_owned_bytes(chrom.sequence(), 2, 17);
        let (hits, _ranges) = consolidate_g4s(raw);
        expected.insert(chrom.name().to_string(), hits);
    }
    let mut actual: HashMap<String, Vec<_>> = HashMap::new();
    stream::process_fasta_stream(&path, 2, 17, |name, results| {
        actual.insert(name, results);
        Ok(())
    })
    .unwrap();
    assert_eq!(expected.len(), actual.len());
    for (name, expected_hits) in expected {
        let observed = actual.get(&name).expect("missing chromosome");
        assert_eq!(expected_hits.len(), observed.len());
        for (lhs, rhs) in expected_hits.iter().zip(observed.iter()) {
            assert_eq!(lhs.start, rhs.start);
            assert_eq!(lhs.end, rhs.end);
            assert_eq!(lhs.sequence(), rhs.sequence());
            assert_eq!(lhs.tetrads, rhs.tetrads);
            assert_eq!(lhs.gscore, rhs.gscore);
        }
    }
    fs::remove_file(&path).unwrap();
}

#[test]
fn stream_overlap_exposes_metadata() {
    let path = std::env::temp_dir().join("qgrs_stream_overlap.fa");
    let fasta = b">chr1 desc\nGGGGAGGGGTTTTGGGG\n";
    fs::write(&path, fasta).unwrap();
    let mut observed = Vec::new();
    stream::process_fasta_stream_with_limits_overlap(
        &path,
        2,
        17,
        ScanLimits::default(),
        |name, results| {
            assert_eq!(name, "chr1");
            assert!(!results.hits.is_empty());
            assert!(
                results
                    .raw_hits
                    .as_ref()
                    .expect("overlap results missing raw hits")
                    .len()
                    >= results.hits.len()
            );
            assert_eq!(results.family_ranges.len(), results.hits.len());
            observed.push(results.hits.len());
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(observed.len(), 1);
    fs::remove_file(&path).unwrap();
}

#[test]
fn stream_pipeline_matches_batch_results_in_circular_mode() {
    let path = std::env::temp_dir().join("qgrs_stream_pipeline_circular.fa");
    let fasta = b">chr1\nGAGGGGAGGGGAGGGGGGG\n>chr2\nACACGGGGAGGGGAGGGGGGGAC\n";
    fs::write(&path, fasta).unwrap();
    let limits = ScanLimits::default();

    let sequences = crate::qgrs::load_sequences_from_path(&path, InputMode::Stream).unwrap();
    let mut expected: HashMap<String, Vec<_>> = HashMap::new();
    for chrom in &sequences {
        let seq_len = chrom.sequence().len();
        let raw = find_owned_bytes_with_topology(
            chrom.sequence(),
            4,
            17,
            limits,
            SequenceTopology::Circular,
        );
        let (hits, _ranges) =
            consolidate_g4s_with_topology(raw, SequenceTopology::Circular, seq_len);
        expected.insert(chrom.name().to_string(), hits);
    }

    let mut actual: HashMap<String, Vec<_>> = HashMap::new();
    stream::process_fasta_stream_with_limits_topology(
        &path,
        4,
        17,
        limits,
        SequenceTopology::Circular,
        |name, results| {
            actual.insert(name, results);
            Ok(())
        },
    )
    .unwrap();

    assert_eq!(expected.len(), actual.len());
    for (name, expected_hits) in expected {
        let observed = actual.get(&name).expect("missing chromosome");
        assert_eq!(expected_hits.len(), observed.len());
        for (lhs, rhs) in expected_hits.iter().zip(observed.iter()) {
            assert_eq!(lhs.start, rhs.start);
            assert_eq!(lhs.end, rhs.end);
            assert_eq!(lhs.length, rhs.length);
            assert_eq!(lhs.sequence(), rhs.sequence());
            assert_eq!(lhs.tetrads, rhs.tetrads);
            assert_eq!(lhs.gscore, rhs.gscore);
            assert_eq!(lhs.y1, rhs.y1);
            assert_eq!(lhs.y2, rhs.y2);
            assert_eq!(lhs.y3, rhs.y3);
        }
    }
    fs::remove_file(&path).unwrap();
}
