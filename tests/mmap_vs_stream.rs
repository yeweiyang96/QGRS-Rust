use std::path::Path;
use std::sync::{Arc, Mutex};

use qgrs_rust::qgrs::stream as qstream;
use qgrs_rust::qgrs::{self, find_owned_with_limits, load_sequences_from_path, render_csv_results};
use qgrs_rust::{InputMode, ScanLimits};

#[test]
fn mmap_and_stream_produce_same_results_on_dme() {
    let path = Path::new("dme.fa");
    let min_tetrads = 2usize;
    let min_score = 17i32;
    let limits = ScanLimits::default();

    // Collect mmap results by loading sequences and scanning each
    let mut mmap_lines: Vec<String> = Vec::new();
    let sequences =
        load_sequences_from_path(path, InputMode::Mmap).expect("failed to load sequences via mmap");
    for chrom in sequences {
        let results = find_owned_with_limits(chrom.sequence, min_tetrads, min_score, limits);
        let csv = render_csv_results(&results);
        for (i, line) in csv.lines().enumerate() {
            if i == 0 {
                continue;
            } // skip header
            mmap_lines.push(line.to_string());
        }
    }

    // Collect stream results by processing the file and capturing callback outputs
    let stream_lines: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let collector = stream_lines.clone();
    qstream::process_fasta_stream_with_limits(
        path,
        min_tetrads,
        min_score,
        limits,
        move |_name, results| {
            let csv = render_csv_results(&results);
            let mut guard = collector.lock().unwrap();
            for (i, line) in csv.lines().enumerate() {
                if i == 0 {
                    continue;
                }
                guard.push(line.to_string());
            }
            Ok(())
        },
    )
    .expect("stream processing failed");

    let mut mmap_lines_sorted = mmap_lines;
    mmap_lines_sorted.sort();
    let mut stream_lines_sorted = stream_lines.lock().unwrap().clone();
    stream_lines_sorted.sort();

    assert_eq!(
        mmap_lines_sorted, stream_lines_sorted,
        "mmap vs stream outputs differ"
    );
}
