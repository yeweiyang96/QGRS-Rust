use std::collections::HashMap;
use std::fs;

use crate::qgrs::stream;
use crate::qgrs::{InputMode, consolidate_g4s, find_owned_bytes};

#[test]
fn stream_pipeline_matches_batch_results() {
    let path = std::env::temp_dir().join("qgrs_stream_pipeline.fa");
    let fasta = b">chr1 desc\nGGGGAGGGGTTTTGGGG\n>chr2\nACACGGGGACACGGGG\n";
    fs::write(&path, fasta).unwrap();
    let sequences = crate::qgrs::load_sequences_from_path(&path, InputMode::Stream).unwrap();
    let mut expected: HashMap<String, Vec<_>> = HashMap::new();
    for chrom in &sequences {
        let raw = find_owned_bytes(chrom.sequence(), 2, 17);
        let hits = consolidate_g4s(raw);
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
