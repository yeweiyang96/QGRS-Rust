use qgrs_rust::qgrs::{
    self, InputMode, ScanLimits, find_owned_bytes_with_limits, load_sequences_from_path, stream,
};
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_path> [min_tetrads] [min_gscore]", args[0]);
        eprintln!("\nExamples:");
        eprintln!("  {} dme.fa", args[0]);
        eprintln!("  {} dme.fa 3 17", args[0]);
        std::process::exit(1);
    }

    let path = PathBuf::from(&args[1]);
    let min_tetrads = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2);
    let min_gscore = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(17);

    if !path.exists() {
        eprintln!("❌ File does not exist: {:?}", path);
        std::process::exit(1);
    }

    println!("════════════════════════════════════════════════════════");
    println!("🔬 QGRS Stream vs Mmap Mode Performance Comparison");
    println!("════════════════════════════════════════════════════════");
    println!("File: {}", path.display());
    println!(
        "Parameters: min_tetrads={}, min_gscore={}",
        min_tetrads, min_gscore
    );
    println!("════════════════════════════════════════════════════════\n");

    // ========== Test Batch/Mmap mode ==========
    println!("⏳ Running Batch/Mmap mode...");
    let start = Instant::now();

    let sequences = match load_sequences_from_path(&path, InputMode::Mmap) {
        Ok(seqs) => seqs,
        Err(e) => {
            eprintln!("❌ Failed to load sequences: {}", e);
            std::process::exit(1);
        }
    };

    let load_time = start.elapsed();
    println!("  ✓ Sequence loading complete: {:?}", load_time);

    let process_start = Instant::now();
    let mut batch_results: HashMap<String, Vec<_>> = HashMap::new();
    let limits = ScanLimits::default();
    for chrom in &sequences {
        let hits = find_owned_bytes_with_limits(chrom.sequence(), min_tetrads, min_gscore, limits);
        batch_results.insert(chrom.name().to_string(), hits);
    }
    let process_time = process_start.elapsed();
    let batch_total_time = start.elapsed();

    let batch_total_hits: usize = batch_results.values().map(|v| v.len()).sum();

    println!("  ✓ Sequence processing complete: {:?}", process_time);
    println!("\n📊 Batch/Mmap mode results:");
    println!("  Chromosome count: {}", batch_results.len());
    println!("  Total G4s: {}", batch_total_hits);
    println!("  Loading time: {:?}", load_time);
    println!("  Processing time: {:?}", process_time);
    println!("  Total time: {:?}", batch_total_time);

    // Display detailed info for each chromosome
    let mut chrom_list: Vec<_> = batch_results.iter().collect();
    chrom_list.sort_by_key(|(name, _)| name.as_str());
    println!("\n  Detailed results:");
    for (name, hits) in &chrom_list {
        println!("    {}: {} G4s", name, hits.len());
    }

    println!("\n════════════════════════════════════════════════════════\n");

    // ========== Test Stream mode ==========
    println!("⏳ Running Stream mode...");
    let start = Instant::now();

    let mut stream_results: HashMap<String, Vec<_>> = HashMap::new();
    if let Err(e) = stream::process_fasta_stream(&path, min_tetrads, min_gscore, |name, results| {
        stream_results.insert(name, results);
        Ok(())
    }) {
        eprintln!("❌ Stream processing failed: {}", e);
        std::process::exit(1);
    }

    let stream_total_time = start.elapsed();
    let stream_total_hits: usize = stream_results.values().map(|v| v.len()).sum();

    println!("  ✓ Processing complete");
    println!("\n📊 Stream mode results:");
    println!("  Chromosome count: {}", stream_results.len());
    println!("  Total G4s: {}", stream_total_hits);
    println!("  Total time: {:?}", stream_total_time);

    // Display detailed info for each chromosome
    let mut stream_chrom_list: Vec<_> = stream_results.iter().collect();
    stream_chrom_list.sort_by_key(|(name, _)| name.as_str());
    println!("\n  Detailed results:");
    for (name, hits) in &stream_chrom_list {
        println!("    {}: {} G4s", name, hits.len());
    }

    println!("\n════════════════════════════════════════════════════════\n");

    // ========== Performance comparison ==========
    println!("⚡ Performance comparison:");
    let speedup = batch_total_time.as_secs_f64() / stream_total_time.as_secs_f64();
    println!("  Batch/Mmap: {:?}", batch_total_time);
    println!("  Stream:     {:?}", stream_total_time);
    if speedup > 1.0 {
        println!("  Stream is {:.2}x faster", speedup);
    } else {
        println!("  Batch/Mmap is {:.2}x faster", 1.0 / speedup);
    }

    println!("\n════════════════════════════════════════════════════════\n");

    // ========== Consistency verification ==========
    println!("🔍 Verifying result consistency...");

    let mut mismatches = 0;
    let mut details = Vec::new();

    // Check chromosome count
    if batch_results.len() != stream_results.len() {
        details.push(format!(
            "  ⚠️  Chromosome count mismatch: Batch={}, Stream={}",
            batch_results.len(),
            stream_results.len()
        ));
        mismatches += 1;
    }

    // Check total G4 count
    if batch_total_hits != stream_total_hits {
        details.push(format!(
            "  ⚠️  Total G4 count mismatch: Batch={}, Stream={}",
            batch_total_hits, stream_total_hits
        ));
        mismatches += 1;
    }

    // Check each chromosome
    for (name, batch_hits) in &batch_results {
        if let Some(stream_hits) = stream_results.get(name) {
            if batch_hits.len() != stream_hits.len() {
                details.push(format!(
                    "  ⚠️  G4 count mismatch for chromosome {}: Batch={}, Stream={}",
                    name,
                    batch_hits.len(),
                    stream_hits.len()
                ));
                mismatches += 1;
            } else {
                // Compare G4 details one by one
                for (i, (batch_g4, stream_g4)) in
                    batch_hits.iter().zip(stream_hits.iter()).enumerate()
                {
                    // Compare all fields to ensure complete consistency
                    if batch_g4.start != stream_g4.start
                        || batch_g4.end != stream_g4.end
                        || batch_g4.sequence() != stream_g4.sequence()
                        || batch_g4.tetrads != stream_g4.tetrads
                        || batch_g4.gscore != stream_g4.gscore
                        || batch_g4.length != stream_g4.length
                        || batch_g4.y1 != stream_g4.y1
                        || batch_g4.y2 != stream_g4.y2
                        || batch_g4.y3 != stream_g4.y3
                    {
                        details.push(format!(
                            "  ⚠️  G4 #{} mismatch in chromosome {}:",
                            i + 1,
                            name
                        ));
                        details.push(format!(
                            "      Batch:  pos={}..{}, len={}, seq={}, tetrads={}, y=({},{},{}), gscore={}",
                            batch_g4.start,
                            batch_g4.end,
                            batch_g4.length,
                            batch_g4.sequence(),
                            batch_g4.tetrads,
                            batch_g4.y1,
                            batch_g4.y2,
                            batch_g4.y3,
                            batch_g4.gscore
                        ));
                        details.push(format!(
                            "      Stream: pos={}..{}, len={}, seq={}, tetrads={}, y=({},{},{}), gscore={}",
                            stream_g4.start,
                            stream_g4.end,
                            stream_g4.length,
                            stream_g4.sequence(),
                            stream_g4.tetrads,
                            stream_g4.y1,
                            stream_g4.y2,
                            stream_g4.y3,
                            stream_g4.gscore
                        ));
                        mismatches += 1;
                        if mismatches >= 10 {
                            details.push("  ... (additional mismatches omitted)".to_string());
                            break;
                        }
                    }
                }
            }
        } else {
            details.push(format!("  ⚠️  Stream mode missing chromosome: {}", name));
            mismatches += 1;
        }

        if mismatches >= 10 {
            break;
        }
    }

    // Check if Stream has additional chromosomes
    if mismatches < 10 {
        for name in stream_results.keys() {
            if !batch_results.contains_key(name) {
                details.push(format!("  ⚠️  Batch mode missing chromosome: {}", name));
                mismatches += 1;
                if mismatches >= 10 {
                    break;
                }
            }
        }
    }

    if mismatches == 0 {
        println!("  ✅ All results are completely consistent!");
        println!("     - Chromosome count: {}", batch_results.len());
        println!("     - Total G4s: {}", batch_total_hits);
        println!("     - All G4 fields (position, length, sequence, tetrads, loops, gscore) match");
    } else {
        println!("  ❌ Found {} mismatch(es):", mismatches);
        for detail in details {
            println!("{}", detail);
        }
        println!("\n════════════════════════════════════════════════════════");
        std::process::exit(1);
    }

    println!("\n════════════════════════════════════════════════════════");
    println!("✅ Test completed!");
    println!("════════════════════════════════════════════════════════");
}
