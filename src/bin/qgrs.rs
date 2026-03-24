use std::collections::HashMap;
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use qgrs_rust::qgrs::{
    self, DEFAULT_MAX_G_RUN, DEFAULT_MAX_G4_LENGTH, G4, InputMode, ScanLimits, SequenceTopology,
};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

fn main() {
    // Initialize Rayon global thread pool to match machine CPU count.
    // This makes parallelism deterministic across runs and avoids relying on
    // the environment variable `RAYON_NUM_THREADS`.
    let threads = num_cpus::get();
    let _ = ThreadPoolBuilder::new().num_threads(threads).build_global();

    if let Err(err) = run_env(env::args().skip(1)) {
        eprintln!("Error: {err}");
        std::process::exit(1);
    }
}

fn run_env<I>(mut args: I) -> Result<(), String>
where
    I: Iterator<Item = String>,
{
    let mut sequence_arg: Option<String> = None;
    let mut file_arg: Option<PathBuf> = None;
    let mut min_tetrads: usize = 2;
    let mut min_score: i32 = 17;
    let mut max_g_run: usize = DEFAULT_MAX_G_RUN;
    let mut max_g4_length: usize = DEFAULT_MAX_G4_LENGTH;
    let mut format = OutputFormat::Csv;
    let mut output_path: Option<PathBuf> = None;
    let mut output_dir: Option<PathBuf> = None;
    let mut mode = InputMode::Mmap;
    let mut include_overlap = false;
    let mut include_revcomp = false;
    let mut circular = false;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--sequence" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --sequence"))?;
                sequence_arg = Some(value);
            }
            "--file" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --file"))?;
                file_arg = Some(PathBuf::from(value));
            }
            "--min-tetrads" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --min-tetrads"))?
                    .parse::<usize>()
                    .map_err(|_| usage("--min-tetrads must be a positive integer"))?;
                if value == 0 {
                    return Err(usage("--min-tetrads must be > 0"));
                }
                min_tetrads = value;
            }
            "--min-score" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --min-score"))?
                    .parse::<i32>()
                    .map_err(|_| usage("--min-score must be an integer"))?;
                min_score = value;
            }
            "--format" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --format"))?;
                format = value.try_into()?;
            }
            "--mode" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --mode"))?;
                mode = parse_mode(&value)?;
            }
            "--max-g-run" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --max-g-run"))?
                    .parse::<usize>()
                    .map_err(|_| usage("--max-g-run must be a positive integer"))?;
                if value == 0 {
                    return Err(usage("--max-g-run must be > 0"));
                }
                max_g_run = value;
            }
            "--max-g4-length" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --max-g4-length"))?
                    .parse::<usize>()
                    .map_err(|_| usage("--max-g4-length must be a positive integer"))?;
                if value == 0 {
                    return Err(usage("--max-g4-length must be > 0"));
                }
                max_g4_length = value;
            }
            "--output" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --output"))?;
                output_path = Some(PathBuf::from(value));
            }
            "--output-dir" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --output-dir"))?;
                output_dir = Some(PathBuf::from(value));
            }
            "--overlap" => {
                include_overlap = true;
            }
            "--revcomp" => {
                include_revcomp = true;
            }
            "--circular" => {
                circular = true;
            }
            "--help" | "-h" => return Err(usage("")),
            other => {
                return Err(usage(&format!("unknown argument '{other}'")));
            }
        }
    }

    let input = match (sequence_arg, file_arg) {
        (Some(_), Some(_)) => {
            return Err(usage("cannot provide both --sequence and --file"));
        }
        (Some(seq), None) => InputSpec::Inline(seq),
        (None, Some(path)) => InputSpec::File(path),
        (None, None) => return Err(usage("must provide --sequence or --file")),
    };

    let min_required_length = min_tetrads
        .checked_mul(4)
        .ok_or_else(|| usage("--min-tetrads is too large"))?;
    if max_g_run < min_tetrads {
        return Err(usage("--max-g-run must be ≥ --min-tetrads"));
    }
    if max_g4_length < min_required_length {
        return Err(usage("--max-g4-length must be ≥ 4 * --min-tetrads"));
    }

    let limits = ScanLimits::new(max_g4_length, max_g_run);
    let topology = if circular {
        SequenceTopology::Circular
    } else {
        SequenceTopology::Linear
    };
    let scan = ScanConfig::new(min_tetrads, min_score, limits, topology);

    match input {
        InputSpec::Inline(seq) => {
            if output_dir.is_some() {
                return Err(usage("--output-dir can only be used with --file"));
            }
            process_inline_sequence(
                seq,
                format,
                output_path,
                scan,
                include_overlap,
                include_revcomp,
            )?;
        }
        InputSpec::File(path) => {
            if output_path.is_some() {
                return Err(usage(
                    "--output is only valid with --sequence; use --output-dir for --file",
                ));
            }
            process_fasta_file(
                path,
                mode,
                format,
                scan,
                output_dir,
                include_overlap,
                include_revcomp,
            )?;
        }
    }
    Ok(())
}

fn usage(reason: &str) -> String {
    let mut msg = String::new();
    if !reason.is_empty() {
        msg.push_str(reason);
        msg.push('\n');
    }
    msg.push_str("Usage: cargo run --bin qgrs -- [--sequence <SEQ> | --file <PATH>] [options]\n");
    msg.push_str("Options:\n");
    msg.push_str("  --sequence <SEQ>     Inline DNA/RNA sequence to scan\n");
    msg.push_str(
        "  --file <PATH>        Read sequences from FASTA/FASTA.gz (chromosomes split independently)\n",
    );
    msg.push_str("  --min-tetrads <N>    Minimum tetrads to seed (default 2)\n");
    msg.push_str("  --min-score <S>      Minimum g-score (default 17)\n");
    msg.push_str("  --max-g-run <N>      Maximum allowed G-run length (default 10)\n");
    msg.push_str("  --max-g4-length <N>  Maximum allowed G4 length in bp (default 45)\n");
    msg.push_str("  --format <csv|parquet>  Output format (default csv)\n");
    msg.push_str(
        "  --output <PATH>     Destination file when using --sequence (required for parquet)\n",
    );
    msg.push_str("  --output-dir <DIR>  Directory for per-chromosome exports when using --file\n");
    msg.push_str("  --mode <mmap|stream> Input mode when using --file (default mmap)\n");
    msg.push_str(
        "  --overlap            Emit raw hits (.overlap.<format>) and family ranges (.family.<format>)\n",
    );
    msg.push_str(
        "  --revcomp            Also scan reverse-complement and emit .revcomp.<format>\n",
    );
    msg.push_str("  --circular           Treat each sequence/chromosome as circular\n");
    msg.push_str("  --help               Show this message\n");
    msg
}

fn parse_mode(value: &str) -> Result<InputMode, String> {
    match value {
        "mmap" => Ok(InputMode::Mmap),
        "stream" => Ok(InputMode::Stream),
        _ => Err(usage("--mode must be either 'mmap' or 'stream'")),
    }
}

enum InputSpec {
    Inline(String),
    File(PathBuf),
}

#[derive(Clone, Copy)]
struct ScanConfig {
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
}

impl ScanConfig {
    fn new(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
    ) -> Self {
        Self {
            min_tetrads,
            min_score,
            limits,
            topology,
        }
    }

    fn min_tetrads(self) -> usize {
        self.min_tetrads
    }

    fn min_score(self) -> i32 {
        self.min_score
    }

    fn limits(self) -> ScanLimits {
        self.limits
    }

    fn topology(self) -> SequenceTopology {
        self.topology
    }
}

fn process_inline_sequence(
    sequence: String,
    format: OutputFormat,
    output_path: Option<PathBuf>,
    scan: ScanConfig,
    include_overlap: bool,
    include_revcomp: bool,
) -> Result<(), String> {
    let mut normalized = sequence.into_bytes();
    normalized.make_ascii_lowercase();
    let sequence_len = normalized.len();
    if include_overlap && output_path.is_none() {
        return Err(usage("--overlap requires --output when using --sequence"));
    }
    if include_revcomp && output_path.is_none() {
        return Err(usage("--revcomp requires --output when using --sequence"));
    }

    let (results, family_ranges, raw_hits) = run_scan_for_export(
        Arc::new(normalized.clone()),
        scan,
        include_overlap,
        sequence_len,
    );
    write_primary_output(
        output_path.as_deref(),
        format,
        &results,
        scan.topology(),
        sequence_len,
    )?;

    if include_overlap {
        let base = output_path
            .as_ref()
            .expect("overlap outputs require an explicit --output path");
        write_overlap_exports(
            base,
            format,
            raw_hits.as_ref().unwrap(),
            &family_ranges,
            scan.topology(),
            sequence_len,
        )?;
    }

    if include_revcomp {
        let base = output_path
            .as_ref()
            .expect("revcomp outputs require an explicit --output path");
        let (revcomp_results, revcomp_family_ranges, revcomp_raw_hits) =
            run_revcomp_scan_for_export(&normalized, scan, include_overlap);
        let revcomp_path = revcomp_output_path(base, format);
        write_results_to_path(
            &revcomp_path,
            format,
            &revcomp_results,
            scan.topology(),
            sequence_len,
        )?;
        if include_overlap {
            let revcomp_raw_hits = revcomp_raw_hits
                .as_ref()
                .expect("revcomp raw hits must be captured when overlap is requested");
            write_overlap_exports(
                &revcomp_path,
                format,
                revcomp_raw_hits,
                &revcomp_family_ranges,
                scan.topology(),
                sequence_len,
            )?;
        }
    }
    Ok(())
}

fn process_fasta_file(
    path: PathBuf,
    mode: InputMode,
    format: OutputFormat,
    scan: ScanConfig,
    output_dir: Option<PathBuf>,
    include_overlap: bool,
    include_revcomp: bool,
) -> Result<(), String> {
    let dir = output_dir.ok_or_else(|| usage("--output-dir is required when --file is used"))?;
    fs::create_dir_all(&dir).map_err(|err| format!("failed to create {dir:?}: {err}"))?;
    let mut name_counts: HashMap<String, usize> = HashMap::new();
    match mode {
        InputMode::Mmap => {
            let sequences = qgrs::load_sequences_from_path(&path, InputMode::Mmap)
                .map_err(|err| format!("failed to read {path:?}: {err}"))?;
            if sequences.is_empty() {
                return Err(format!("no sequences found in {path:?}"));
            }
            let mut chrom_outputs = Vec::with_capacity(sequences.len());
            for chrom in sequences {
                let filename = next_output_filename(chrom.name(), format, &mut name_counts);
                chrom_outputs.push((chrom, dir.join(filename)));
            }
            chrom_outputs.into_par_iter().try_for_each(
                |(chrom, filepath)| -> Result<(), String> {
                    let (_name, sequence) = chrom.into_parts();
                    let sequence_len = sequence.len();
                    let (results, family_ranges, raw_hits) =
                        run_scan_for_export(sequence.clone(), scan, include_overlap, sequence_len);
                    write_results_to_path(
                        &filepath,
                        format,
                        &results,
                        scan.topology(),
                        sequence_len,
                    )?;
                    if include_overlap {
                        let raw_hits = raw_hits
                            .as_ref()
                            .expect("raw hits must be captured when overlap is requested");
                        write_overlap_exports(
                            &filepath,
                            format,
                            raw_hits,
                            &family_ranges,
                            scan.topology(),
                            sequence_len,
                        )?;
                    }
                    if include_revcomp {
                        let (revcomp_results, revcomp_ranges, revcomp_raw_hits) =
                            run_revcomp_scan_for_export(sequence.as_slice(), scan, include_overlap);
                        let revcomp_path = revcomp_output_path(&filepath, format);
                        write_results_to_path(
                            &revcomp_path,
                            format,
                            &revcomp_results,
                            scan.topology(),
                            sequence_len,
                        )?;
                        if include_overlap {
                            let revcomp_raw_hits = revcomp_raw_hits.as_ref().expect(
                                "revcomp raw hits must be captured when overlap is requested",
                            );
                            write_overlap_exports(
                                &revcomp_path,
                                format,
                                revcomp_raw_hits,
                                &revcomp_ranges,
                                scan.topology(),
                                sequence_len,
                            )?;
                        }
                    }
                    Ok(())
                },
            )?;
        }
        InputMode::Stream => {
            let mut processed = 0usize;
            if include_revcomp && include_overlap {
                qgrs::stream::process_fasta_stream_with_limits_overlap_topology_and_sequence(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    |name, mut stream_results, sequence| {
                        processed += 1;
                        let sequence_len = sequence.len();
                        let filename = next_output_filename(&name, format, &mut name_counts);
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &stream_results.hits,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        let raw_hits = stream_results
                            .raw_hits
                            .take()
                            .expect("raw hits missing from overlap stream results");
                        write_overlap_exports(
                            &filepath,
                            format,
                            &raw_hits,
                            &stream_results.family_ranges,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;

                        let (revcomp_results, revcomp_ranges, revcomp_raw_hits) =
                            run_revcomp_scan_for_export(&sequence, scan, true);
                        let revcomp_path = revcomp_output_path(&filepath, format);
                        write_results_to_path(
                            &revcomp_path,
                            format,
                            &revcomp_results,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        let revcomp_raw_hits = revcomp_raw_hits
                            .as_ref()
                            .expect("revcomp raw hits must be captured when overlap is requested");
                        write_overlap_exports(
                            &revcomp_path,
                            format,
                            revcomp_raw_hits,
                            &revcomp_ranges,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            } else if include_revcomp {
                qgrs::stream::process_fasta_stream_with_limits_topology_and_sequence(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    |name, results, sequence| {
                        processed += 1;
                        let sequence_len = sequence.len();
                        let filename = next_output_filename(&name, format, &mut name_counts);
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &results,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;

                        let (revcomp_results, _, _) =
                            run_revcomp_scan_for_export(&sequence, scan, false);
                        let revcomp_path = revcomp_output_path(&filepath, format);
                        write_results_to_path(
                            &revcomp_path,
                            format,
                            &revcomp_results,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            } else if include_overlap {
                qgrs::stream::process_fasta_stream_with_limits_overlap_topology_and_len(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    |name, mut stream_results, sequence_len| {
                        processed += 1;
                        let filename = next_output_filename(&name, format, &mut name_counts);
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &stream_results.hits,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        let raw_hits = stream_results
                            .raw_hits
                            .take()
                            .expect("raw hits missing from overlap stream results");

                        write_overlap_exports(
                            &filepath,
                            format,
                            &raw_hits,
                            &stream_results.family_ranges,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            } else {
                qgrs::stream::process_fasta_stream_with_limits_topology_and_len(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    |name, results, sequence_len| {
                        processed += 1;
                        let filename = next_output_filename(&name, format, &mut name_counts);
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &results,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            }
            if processed == 0 {
                return Err(format!("no sequences found in {path:?}"));
            }
        }
    }
    Ok(())
}

fn next_output_filename(
    name: &str,
    format: OutputFormat,
    counts: &mut HashMap<String, usize>,
) -> String {
    let sanitized = sanitize_name(name);
    // 处理同名染色体的重复输出(万一)
    let entry = counts.entry(sanitized.clone()).or_insert(0);
    let suffix = if *entry == 0 {
        String::new()
    } else {
        format!("_{}", entry)
    };
    *entry += 1;
    format!("{}{suffix}.{}", sanitized, format.extension())
}

fn sanitize_name(raw: &str) -> String {
    // let mut sanitized = String::new();
    // for ch in raw.chars() {
    //     if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_') {
    //         sanitized.push(ch);
    //     } else {
    //         sanitized.push('_');
    //     }
    // }
    // if sanitized.is_empty() {
    //     "chromosome".to_string()
    // } else {
    //     sanitized
    // }
    if raw.is_empty() {
        "chromosome".to_string()
    } else {
        raw.to_string()
    }
}

type ConsolidatedResults = (Vec<G4>, Vec<(usize, usize)>, Option<Vec<G4>>);

fn consolidate_for_export(
    raw: Vec<G4>,
    capture_raw: bool,
    topology: SequenceTopology,
    sequence_len: usize,
) -> ConsolidatedResults {
    if capture_raw {
        let raw_copy = raw.clone();
        let (hits, ranges) = qgrs::consolidate_g4s_with_topology(raw, topology, sequence_len);
        (hits, ranges, Some(raw_copy))
    } else {
        let (hits, ranges) = qgrs::consolidate_g4s_with_topology(raw, topology, sequence_len);
        (hits, ranges, None)
    }
}

fn run_scan_for_export(
    sequence: Arc<Vec<u8>>,
    scan: ScanConfig,
    capture_raw: bool,
    sequence_len: usize,
) -> ConsolidatedResults {
    let raw = qgrs::find_owned_bytes_with_topology(
        sequence,
        scan.min_tetrads(),
        scan.min_score(),
        scan.limits(),
        scan.topology(),
    );
    consolidate_for_export(raw, capture_raw, scan.topology(), sequence_len)
}

fn run_revcomp_scan_for_export(
    forward_sequence: &[u8],
    scan: ScanConfig,
    capture_raw: bool,
) -> ConsolidatedResults {
    let sequence_len = forward_sequence.len();
    let revcomp = Arc::new(reverse_complement_lowercase(forward_sequence));
    let (hits, ranges, raw_hits) = run_scan_for_export(revcomp, scan, capture_raw, sequence_len);
    let mapped_hits = map_revcomp_hits_to_forward(hits, scan.topology(), sequence_len);
    let mapped_ranges = map_revcomp_ranges_to_forward(ranges, scan.topology(), sequence_len);
    let mapped_raw_hits =
        raw_hits.map(|hits| map_revcomp_hits_to_forward(hits, scan.topology(), sequence_len));
    (mapped_hits, mapped_ranges, mapped_raw_hits)
}

fn reverse_complement_lowercase(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|byte| match *byte {
            b'a' | b'A' => b't',
            b't' | b'T' | b'u' | b'U' => b'a',
            b'c' | b'C' => b'g',
            b'g' | b'G' => b'c',
            b'n' | b'N' => b'n',
            other => other.to_ascii_lowercase(),
        })
        .collect()
}

fn map_revcomp_hits_to_forward(
    mut hits: Vec<G4>,
    topology: SequenceTopology,
    sequence_len: usize,
) -> Vec<G4> {
    if sequence_len == 0 {
        return hits;
    }
    for hit in &mut hits {
        let (start, end) = map_revcomp_interval(hit.start, hit.end, topology, sequence_len);
        hit.start = start;
        hit.end = end;
    }
    hits.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));
    hits
}

fn map_revcomp_ranges_to_forward(
    mut ranges: Vec<(usize, usize)>,
    topology: SequenceTopology,
    sequence_len: usize,
) -> Vec<(usize, usize)> {
    if sequence_len == 0 {
        return ranges;
    }
    for range in &mut ranges {
        *range = map_revcomp_interval(range.0, range.1, topology, sequence_len);
    }
    ranges.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    ranges
}

fn map_revcomp_interval(
    start_rc: usize,
    end_rc: usize,
    topology: SequenceTopology,
    sequence_len: usize,
) -> (usize, usize) {
    let end_rc_projected = projected_end(end_rc, topology, sequence_len);
    let start = sequence_len - end_rc_projected + 1;
    let end = sequence_len - start_rc + 1;
    (start, end)
}

fn projected_end(end: usize, topology: SequenceTopology, sequence_len: usize) -> usize {
    if !topology.is_circular() || sequence_len == 0 || end == 0 {
        return end;
    }
    ((end - 1) % sequence_len) + 1
}

fn write_primary_output(
    output_path: Option<&Path>,
    format: OutputFormat,
    results: &[G4],
    topology: SequenceTopology,
    sequence_len: usize,
) -> Result<(), String> {
    match format {
        OutputFormat::Csv => {
            let csv = qgrs::render_csv_results_with_projection(results, topology, sequence_len);
            if let Some(path) = output_path {
                fs::write(path, csv).map_err(|err| format!("failed to write {path:?}: {err}"))?;
            } else {
                print!("{csv}");
            }
            Ok(())
        }
        OutputFormat::Parquet => {
            let path =
                output_path.ok_or_else(|| usage("--output is required when --format parquet"))?;
            write_results_to_path(path, format, results, topology, sequence_len)
        }
    }
}

fn write_results_to_path(
    path: &Path,
    format: OutputFormat,
    results: &[G4],
    topology: SequenceTopology,
    sequence_len: usize,
) -> Result<(), String> {
    match format {
        OutputFormat::Csv => {
            let csv = qgrs::render_csv_results_with_projection(results, topology, sequence_len);
            fs::write(path, csv).map_err(|err| format!("failed to write {path:?}: {err}"))?;
        }
        OutputFormat::Parquet => {
            let file = fs::File::create(path)
                .map_err(|err| format!("failed to create {path:?}: {err}"))?;
            qgrs::write_parquet_results_with_projection(results, file, topology, sequence_len)
                .map_err(|err| format!("failed to write parquet {path:?}: {err}"))?;
        }
    }
    Ok(())
}

fn write_overlap_exports(
    base: &Path,
    format: OutputFormat,
    raw_hits: &[G4],
    family_ranges: &[(usize, usize)],
    topology: SequenceTopology,
    sequence_len: usize,
) -> Result<(), String> {
    let overlap_path = overlap_path(base, format);
    let family_path = family_path(base, format);
    match format {
        OutputFormat::Csv => {
            let overlap_csv =
                qgrs::render_csv_results_with_projection(raw_hits, topology, sequence_len);
            fs::write(&overlap_path, overlap_csv)
                .map_err(|err| format!("failed to write {overlap_path:?}: {err}"))?;

            let family_csv = qgrs::render_family_ranges_csv_with_projection(
                family_ranges,
                topology,
                sequence_len,
            );
            fs::write(&family_path, family_csv)
                .map_err(|err| format!("failed to write {family_path:?}: {err}"))?;
        }
        OutputFormat::Parquet => {
            let overlap_file = fs::File::create(&overlap_path)
                .map_err(|err| format!("failed to create {overlap_path:?}: {err}"))?;
            qgrs::write_parquet_results_with_projection(
                raw_hits,
                overlap_file,
                topology,
                sequence_len,
            )
            .map_err(|err| format!("failed to write parquet {overlap_path:?}: {err}"))?;

            let family_file = fs::File::create(&family_path)
                .map_err(|err| format!("failed to create {family_path:?}: {err}"))?;
            qgrs::write_parquet_family_ranges_with_projection(
                family_ranges,
                family_file,
                topology,
                sequence_len,
            )
            .map_err(|err| format!("failed to write parquet {family_path:?}: {err}"))?;
        }
    }
    Ok(())
}

fn revcomp_output_path(base: &Path, format: OutputFormat) -> PathBuf {
    append_output_suffix(base, ".revcomp", format)
}

fn overlap_path(base: &Path, format: OutputFormat) -> PathBuf {
    append_output_suffix(base, ".overlap", format)
}

fn family_path(base: &Path, format: OutputFormat) -> PathBuf {
    append_output_suffix(base, ".family", format)
}

fn append_output_suffix(path: &Path, suffix: &str, format: OutputFormat) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
        .unwrap_or("chromosome");
    let name = format!("{stem}{suffix}.{}", format.extension());
    parent.join(name)
}

#[derive(Clone, Copy)]
enum OutputFormat {
    Csv,
    Parquet,
}

impl TryFrom<String> for OutputFormat {
    type Error = String;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        match value.as_str() {
            "csv" => Ok(OutputFormat::Csv),
            "parquet" => Ok(OutputFormat::Parquet),
            _ => Err(usage("--format must be either 'csv' or 'parquet'")),
        }
    }
}

impl OutputFormat {
    fn extension(&self) -> &'static str {
        match self {
            OutputFormat::Csv => "csv",
            OutputFormat::Parquet => "parquet",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    use flate2::Compression;
    use flate2::write::GzEncoder;

    #[test]
    fn default_limits_are_valid() {
        let limits = ScanLimits::default();
        assert_eq!(limits.max_g4_length, DEFAULT_MAX_G4_LENGTH);
        assert_eq!(limits.max_g_run, DEFAULT_MAX_G_RUN);
        assert!(limits.max_g_run >= 2);
        assert!(limits.max_g4_length >= 8);
    }

    #[test]
    fn usage_fails_on_invalid_limits() {
        let err = run_with_args([
            "--sequence",
            "GGGG",
            "--min-tetrads",
            "4",
            "--max-g4-length",
            "12",
            "--max-g-run",
            "2",
        ]);
        assert!(err.is_err());
        let msg = err.unwrap_err().to_string();
        assert!(msg.contains("max-g4-length"));
    }

    #[test]
    fn overlap_requires_output_for_inline() {
        let err = run_with_args(["--sequence", "GGGG", "--overlap"]);
        assert!(err.is_err());
        let msg = err.unwrap_err();
        assert!(msg.contains("--overlap requires --output"));
    }

    #[test]
    fn revcomp_requires_output_for_inline() {
        let err = run_with_args(["--sequence", "GGGG", "--revcomp"]);
        assert!(err.is_err());
        let msg = err.unwrap_err();
        assert!(msg.contains("--revcomp requires --output"));
    }

    #[test]
    fn circular_flag_is_supported_for_inline_scan() {
        let result = run_with_args([
            "--sequence",
            "GAGGGGAGGGGAGGGGGGG",
            "--min-tetrads",
            "4",
            "--min-score",
            "17",
            "--circular",
        ]);
        assert!(result.is_ok());
    }

    #[test]
    fn circular_cli_outputs_map_coordinates_back_to_ring() {
        let base = unique_test_path("qgrs_circular_cli");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "GAGGGGAGGGGAGGGGGGG".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--circular".to_string(),
            "--overlap".to_string(),
            "--output".to_string(),
            output_str.clone(),
        ]);
        assert!(result.is_ok());

        let csv = fs::read_to_string(&output).expect("main output");
        assert!(csv.contains("\n17,16,19,4,1,1,1,84,GGGGAGGGGAGGGGAGGGG\n"));

        let overlap =
            fs::read_to_string(overlap_path(&output, OutputFormat::Csv)).expect("overlap output");
        for line in overlap.lines().skip(1) {
            let mut cols = line.split(',');
            let start = cols.next().unwrap().parse::<usize>().unwrap();
            let end = cols.next().unwrap().parse::<usize>().unwrap();
            assert!(start <= 19);
            assert!(end <= 19);
        }
        assert!(overlap.lines().skip(1).any(|line| {
            let mut cols = line.split(',');
            let start = cols.next().unwrap().parse::<usize>().unwrap();
            let end = cols.next().unwrap().parse::<usize>().unwrap();
            end < start
        }));

        let family =
            fs::read_to_string(family_path(&output, OutputFormat::Csv)).expect("family output");
        let family_line = family.lines().nth(1).expect("family row");
        let mut cols = family_line.split(',');
        assert_eq!(cols.next(), Some("1"));
        let start = cols.next().unwrap().parse::<usize>().unwrap();
        let end = cols.next().unwrap().parse::<usize>().unwrap();
        assert!(start <= 19);
        assert!(end <= 19);
        assert!(end < start);

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(overlap_path(&output, OutputFormat::Csv));
        let _ = fs::remove_file(family_path(&output, OutputFormat::Csv));
    }

    #[test]
    fn parquet_overlap_and_family_follow_format() {
        let base = unique_test_path("qgrs_parquet_sidecars");
        let output = base.with_extension("parquet");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "GGGGAGGGGAGGGGAGGGG".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--format".to_string(),
            "parquet".to_string(),
            "--overlap".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let overlap = overlap_path(&output, OutputFormat::Parquet);
        let family = family_path(&output, OutputFormat::Parquet);
        let overlap_meta = fs::metadata(&overlap).expect("overlap parquet output");
        let family_meta = fs::metadata(&family).expect("family parquet output");
        assert!(overlap_meta.len() > 0);
        assert!(family_meta.len() > 0);

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(overlap);
        let _ = fs::remove_file(family);
    }

    #[test]
    fn revcomp_inline_outputs_mapped_coordinates() {
        let base = unique_test_path("qgrs_revcomp_inline");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        // Build a sequence where only reverse-complement scan should find the canonical motif.
        let motif = b"ttggggaggggaggggaggggaaaaaaa";
        let sequence = String::from_utf8(reverse_complement_lowercase(motif)).unwrap();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            sequence,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--revcomp".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let revcomp_output = revcomp_output_path(&output, OutputFormat::Csv);
        let revcomp_csv = fs::read_to_string(&revcomp_output).expect("revcomp output");
        assert!(revcomp_csv.contains("\n8,26,19,4,1,1,1,84,GGGGAGGGGAGGGGAGGGG\n"));

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(revcomp_output);
    }

    #[test]
    fn revcomp_interval_mapping_for_circular_uses_projected_end() {
        let (start, end) = map_revcomp_interval(17, 35, SequenceTopology::Circular, 19);
        assert_eq!(start, 4);
        assert_eq!(end, 3);
        assert!(end < start);
    }

    #[test]
    fn revcomp_with_overlap_outputs_revcomp_sidecars() {
        let base = unique_test_path("qgrs_revcomp_overlap");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        let sequence = String::from_utf8(reverse_complement_lowercase(
            b"ttggggaggggaggggaggggaaaaaaa",
        ))
        .unwrap();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            sequence,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--revcomp".to_string(),
            "--overlap".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let revcomp_output = revcomp_output_path(&output, OutputFormat::Csv);
        assert!(fs::metadata(overlap_path(&revcomp_output, OutputFormat::Csv)).is_ok());
        assert!(fs::metadata(family_path(&revcomp_output, OutputFormat::Csv)).is_ok());

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(overlap_path(&output, OutputFormat::Csv));
        let _ = fs::remove_file(family_path(&output, OutputFormat::Csv));
        let _ = fs::remove_file(&revcomp_output);
        let _ = fs::remove_file(overlap_path(&revcomp_output, OutputFormat::Csv));
        let _ = fs::remove_file(family_path(&revcomp_output, OutputFormat::Csv));
    }

    #[test]
    fn circular_file_outputs_match_between_mmap_and_stream() {
        let fasta = unique_test_path("qgrs_circular_modes").with_extension("fa");
        fs::write(
            &fasta,
            b">chr1\nGAGGGGAGGGGAGGGGGGG\n>chr2\nGGGCGGGGAGGGGAGGGGAG\n",
        )
        .unwrap();
        let mmap_dir = unique_test_path("qgrs_mmap_out");
        let stream_dir = unique_test_path("qgrs_stream_out");
        fs::create_dir_all(&mmap_dir).unwrap();
        fs::create_dir_all(&stream_dir).unwrap();

        let fasta_str = fasta.to_string_lossy().into_owned();
        let mmap_dir_str = mmap_dir.to_string_lossy().into_owned();
        let stream_dir_str = stream_dir.to_string_lossy().into_owned();

        let mmap_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str.clone(),
            "--mode".to_string(),
            "mmap".to_string(),
            "--output-dir".to_string(),
            mmap_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--circular".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(mmap_result.is_ok());

        let stream_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str,
            "--mode".to_string(),
            "stream".to_string(),
            "--output-dir".to_string(),
            stream_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--circular".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(stream_result.is_ok());

        for filename in [
            "chr1.csv",
            "chr1.overlap.csv",
            "chr1.family.csv",
            "chr2.csv",
            "chr2.overlap.csv",
            "chr2.family.csv",
        ] {
            let mmap_contents = fs::read_to_string(mmap_dir.join(filename)).unwrap();
            let stream_contents = fs::read_to_string(stream_dir.join(filename)).unwrap();
            assert_eq!(mmap_contents, stream_contents, "mismatch for {filename}");
        }

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_dir_all(&mmap_dir);
        let _ = fs::remove_dir_all(&stream_dir);
    }

    #[test]
    fn revcomp_file_outputs_match_between_mmap_and_stream() {
        let fasta = unique_test_path("qgrs_revcomp_modes").with_extension("fa");
        fs::write(
            &fasta,
            b">chr1\nAAACCCTCCCCTCCCCTCCCCAAA\n>chr2\nTTTTCCCCCTCCCCTCCCCTCCCCGG\n",
        )
        .unwrap();
        let mmap_dir = unique_test_path("qgrs_revcomp_mmap_out");
        let stream_dir = unique_test_path("qgrs_revcomp_stream_out");
        fs::create_dir_all(&mmap_dir).unwrap();
        fs::create_dir_all(&stream_dir).unwrap();

        let fasta_str = fasta.to_string_lossy().into_owned();
        let mmap_dir_str = mmap_dir.to_string_lossy().into_owned();
        let stream_dir_str = stream_dir.to_string_lossy().into_owned();

        let mmap_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str.clone(),
            "--mode".to_string(),
            "mmap".to_string(),
            "--output-dir".to_string(),
            mmap_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--revcomp".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(mmap_result.is_ok());

        let stream_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str,
            "--mode".to_string(),
            "stream".to_string(),
            "--output-dir".to_string(),
            stream_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--revcomp".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(stream_result.is_ok());

        for filename in [
            "chr1.csv",
            "chr1.overlap.csv",
            "chr1.family.csv",
            "chr1.revcomp.csv",
            "chr1.revcomp.overlap.csv",
            "chr1.revcomp.family.csv",
            "chr2.csv",
            "chr2.overlap.csv",
            "chr2.family.csv",
            "chr2.revcomp.csv",
            "chr2.revcomp.overlap.csv",
            "chr2.revcomp.family.csv",
        ] {
            let mmap_contents = fs::read_to_string(mmap_dir.join(filename)).unwrap();
            let stream_contents = fs::read_to_string(stream_dir.join(filename)).unwrap();
            assert_eq!(mmap_contents, stream_contents, "mismatch for {filename}");
        }

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_dir_all(&mmap_dir);
        let _ = fs::remove_dir_all(&stream_dir);
    }

    #[test]
    fn gzip_file_outputs_match_between_default_mmap_and_stream() {
        let base = unique_test_path("qgrs_gzip_modes");
        let fasta = base.with_extension("fa");
        let fasta_gz = base.with_extension("fna.data");
        let fasta_bytes = b">chr1\nGAGGGGAGGGGAGGGGGGG\n>chr2\nGGGCGGGGAGGGGAGGGGAG\n";
        fs::write(&fasta, fasta_bytes).unwrap();
        write_gzip(&fasta_gz, fasta_bytes);

        let mmap_dir = unique_test_path("qgrs_gzip_mmap_out");
        let stream_dir = unique_test_path("qgrs_gzip_stream_out");
        fs::create_dir_all(&mmap_dir).unwrap();
        fs::create_dir_all(&stream_dir).unwrap();

        let fasta_gz_str = fasta_gz.to_string_lossy().into_owned();
        let mmap_dir_str = mmap_dir.to_string_lossy().into_owned();
        let stream_dir_str = stream_dir.to_string_lossy().into_owned();

        let mmap_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_gz_str.clone(),
            "--output-dir".to_string(),
            mmap_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(mmap_result.is_ok());

        let stream_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_gz_str,
            "--mode".to_string(),
            "stream".to_string(),
            "--output-dir".to_string(),
            stream_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(stream_result.is_ok());

        for filename in [
            "chr1.csv",
            "chr1.overlap.csv",
            "chr1.family.csv",
            "chr2.csv",
            "chr2.overlap.csv",
            "chr2.family.csv",
        ] {
            let mmap_contents = fs::read_to_string(mmap_dir.join(filename)).unwrap();
            let stream_contents = fs::read_to_string(stream_dir.join(filename)).unwrap();
            assert_eq!(mmap_contents, stream_contents, "mismatch for {filename}");
        }

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_file(&fasta_gz);
        let _ = fs::remove_dir_all(&mmap_dir);
        let _ = fs::remove_dir_all(&stream_dir);
    }

    fn run_with_args<const N: usize>(args: [&'static str; N]) -> Result<(), String> {
        let args = args.iter().map(|arg| arg.to_string()).collect::<Vec<_>>();
        run_with_owned_args(args)
    }

    fn run_with_owned_args(args: Vec<String>) -> Result<(), String> {
        let mut argv = vec![String::from("qgrs")];
        argv.extend(args);
        let original = env::args_os().collect::<Vec<_>>();
        let _ = original;
        run_env(argv.into_iter().skip(1))
    }

    fn unique_test_path(prefix: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time before unix epoch")
            .as_nanos();
        env::temp_dir().join(format!("{prefix}_{}_{}", std::process::id(), nonce))
    }

    fn write_gzip(path: &Path, bytes: &[u8]) {
        let file = fs::File::create(path).expect("create gzip file");
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(bytes).expect("write gzip data");
        encoder.finish().expect("finish gzip");
    }
}
