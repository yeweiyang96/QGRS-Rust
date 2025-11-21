use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use atty::Stream;
use qgrs_rust::{
    DEFAULT_MAX_G_RUN, DEFAULT_MAX_G4_LENGTH, InputMode, ProgressPhase, ProgressReporter,
    ProgressReporterHandle, ScanLimits, qgrs,
};
use rayon::current_num_threads;
use rayon::prelude::*;

fn main() {
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
    let progress = TerminalProgressReporter::new_handle();

    match input {
        InputSpec::Inline(seq) => {
            if output_dir.is_some() {
                return Err(usage("--output-dir can only be used with --file"));
            }
            process_inline_sequence(
                seq,
                format,
                output_path,
                min_tetrads,
                min_score,
                limits,
                progress.clone(),
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
                min_tetrads,
                min_score,
                limits,
                output_dir,
                progress.clone(),
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
        "  --file <PATH>        Read sequences from FASTA (chromosomes split independently)\n",
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

const INLINE_PROGRESS_NAME: &str = "inline-sequence";

fn process_inline_sequence(
    sequence: String,
    format: OutputFormat,
    output_path: Option<PathBuf>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    progress: ProgressReporterHandle,
) -> Result<(), String> {
    match format {
        OutputFormat::Csv => {
            let csv = qgrs::find_csv_owned_with_limits_with_progress(
                sequence,
                min_tetrads,
                min_score,
                limits,
                INLINE_PROGRESS_NAME,
                progress.clone(),
            );
            if let Some(path) = output_path {
                fs::write(&path, csv).map_err(|err| format!("failed to write {path:?}: {err}"))?
            } else {
                print!("{csv}");
            }
        }
        OutputFormat::Parquet => {
            let path =
                output_path.ok_or_else(|| usage("--output is required when --format parquet"))?;
            let file = fs::File::create(&path)
                .map_err(|err| format!("failed to create {path:?}: {err}"))?;
            qgrs::write_parquet_owned_with_limits_with_progress(
                sequence,
                min_tetrads,
                min_score,
                limits,
                INLINE_PROGRESS_NAME,
                progress.clone(),
                file,
            )
            .map_err(|err| format!("failed to write parquet: {err}"))?;
        }
    }
    Ok(())
}

fn process_fasta_file(
    path: PathBuf,
    mode: InputMode,
    format: OutputFormat,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    output_dir: Option<PathBuf>,
    progress: ProgressReporterHandle,
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
                let filename = next_output_filename(&chrom.name, format, &mut name_counts);
                chrom_outputs.push((chrom, dir.join(filename)));
            }
            let reporter = progress.clone();
            chrom_outputs.into_par_iter().try_for_each(
                |(chrom, filepath)| -> Result<(), String> {
                    let qgrs::ChromSequence { name, sequence } = chrom;
                    match format {
                        OutputFormat::Csv => {
                            let csv = qgrs::find_csv_owned_with_limits_with_progress(
                                sequence,
                                min_tetrads,
                                min_score,
                                limits,
                                &name,
                                reporter.clone(),
                            );
                            fs::write(&filepath, csv)
                                .map_err(|err| format!("failed to write {filepath:?}: {err}"))?
                        }
                        OutputFormat::Parquet => {
                            let file = fs::File::create(&filepath)
                                .map_err(|err| format!("failed to create {filepath:?}: {err}"))?;
                            qgrs::write_parquet_owned_with_limits_with_progress(
                                sequence,
                                min_tetrads,
                                min_score,
                                limits,
                                &name,
                                reporter.clone(),
                                file,
                            )
                            .map_err(|err| format!("failed to write parquet for {}: {err}", name))?
                        }
                    }
                    Ok(())
                },
            )?;
        }
        InputMode::Stream => {
            let mut processed = 0usize;
            qgrs::stream::process_fasta_stream_with_limits_with_progress(
                &path,
                min_tetrads,
                min_score,
                limits,
                progress.clone(),
                |name, results| {
                    processed += 1;
                    let filename = next_output_filename(&name, format, &mut name_counts);
                    let filepath = dir.join(&filename);
                    match format {
                        OutputFormat::Csv => {
                            let csv = qgrs::render_csv_results(&results);
                            fs::write(&filepath, csv)?;
                        }
                        OutputFormat::Parquet => {
                            let file = fs::File::create(&filepath)?;
                            qgrs::write_parquet_results(&results, file).map_err(|err| {
                                io::Error::new(io::ErrorKind::Other, err.to_string())
                            })?;
                        }
                    }
                    Ok(())
                },
            )
            .map_err(|err| format!("failed to process {path:?}: {err}"))?;
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
    let mut sanitized = String::new();
    for ch in raw.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_') {
            sanitized.push(ch);
        } else {
            sanitized.push('_');
        }
    }
    if sanitized.is_empty() {
        "chromosome".to_string()
    } else {
        sanitized
    }
}

#[derive(Clone, Copy)]
enum OutputFormat {
    Csv,
    Parquet,
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let msg = format!("{}", err.unwrap_err());
        assert!(msg.contains("max-g4-length"));
    }

    fn run_with_args<const N: usize>(args: [&'static str; N]) -> Result<(), String> {
        let mut argv = vec![String::from("qgrs")];
        argv.extend(args.iter().map(|a| a.to_string()));
        let original = env::args_os().collect::<Vec<_>>();
        let _ = original;
        run_env(argv.into_iter().skip(1))
    }
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

struct TerminalProgressReporter {
    enabled: bool,
    state: Mutex<PanelState>,
}

impl TerminalProgressReporter {
    fn new() -> Self {
        let enabled = atty::is(Stream::Stderr);
        let slots = current_num_threads().max(1);
        Self {
            enabled,
            state: Mutex::new(PanelState::new(slots, enabled)),
        }
    }

    fn new_handle() -> ProgressReporterHandle {
        Arc::new(Self::new()) as ProgressReporterHandle
    }

    fn update_line(&self, name: &str, phase: ProgressPhase, processed: usize, total: usize) {
        if !self.enabled {
            return;
        }
        let mut panel = self.state.lock().unwrap();
        panel.update_current_thread(name, phase, processed, total);
    }

    fn mark_idle(&self) {
        if !self.enabled {
            return;
        }
        let mut panel = self.state.lock().unwrap();
        panel.mark_current_thread_idle();
    }
}

impl ProgressReporter for TerminalProgressReporter {
    fn start_chromosome(&self, name: &str, total: usize) {
        self.update_line(name, ProgressPhase::Seeding, 0, total);
    }

    fn update_phase(&self, name: &str, phase: ProgressPhase, processed: usize, total: usize) {
        self.update_line(name, phase, processed, total);
    }

    fn finish_chromosome(&self, _name: &str) {
        self.mark_idle();
    }
}

struct PanelState {
    lines: Vec<LineState>,
    thread_map: HashMap<std::thread::ThreadId, usize>,
    initialized: bool,
    enabled: bool,
}

impl PanelState {
    fn new(slots: usize, enabled: bool) -> Self {
        let count = slots.max(1);
        Self {
            lines: vec![LineState::default(); count],
            thread_map: HashMap::new(),
            initialized: false,
            enabled,
        }
    }

    fn update_current_thread(
        &mut self,
        name: &str,
        phase: ProgressPhase,
        processed: usize,
        total: usize,
    ) {
        if !self.enabled {
            return;
        }
        if let Some(slot) = self.ensure_slot_for_current_thread() {
            self.lines[slot].content = LineContent::Active {
                chromosome: name.to_string(),
                phase,
                processed,
                total: total.max(1),
            };
            self.render();
        }
    }

    fn mark_current_thread_idle(&mut self) {
        if !self.enabled {
            return;
        }
        if let Some(slot) = self.slot_for_current_thread() {
            self.lines[slot].content = LineContent::Idle;
            self.render();
        }
    }

    fn ensure_slot_for_current_thread(&mut self) -> Option<usize> {
        let id = current_thread_id();
        if let Some(&slot) = self.thread_map.get(&id) {
            return Some(slot);
        }
        if let Some(slot) = self.lines.iter().position(|line| line.thread_id.is_none()) {
            self.lines[slot].thread_id = Some(id);
            self.thread_map.insert(id, slot);
            return Some(slot);
        }
        None
    }

    fn slot_for_current_thread(&mut self) -> Option<usize> {
        let id = current_thread_id();
        self.thread_map.get(&id).copied()
    }

    fn render(&mut self) {
        if !self.enabled {
            return;
        }
        let lines = self.lines.len();
        if lines == 0 {
            return;
        }
        let mut stderr = io::stderr();
        if !self.initialized {
            for _ in 0..lines {
                let _ = writeln!(stderr);
            }
            self.initialized = true;
        }
        let _ = write!(stderr, "\x1b[{}F", lines);
        for (idx, line) in self.lines.iter().enumerate() {
            let formatted = line.format(idx);
            let _ = write!(stderr, "\x1b[2K{formatted}\n");
        }
        let _ = stderr.flush();
    }
}

#[derive(Clone)]
struct LineState {
    thread_id: Option<std::thread::ThreadId>,
    content: LineContent,
}

impl Default for LineState {
    fn default() -> Self {
        Self {
            thread_id: None,
            content: LineContent::Idle,
        }
    }
}

impl LineState {
    fn format(&self, slot_index: usize) -> String {
        let label = format!("Thread {:02}", slot_index + 1);
        match &self.content {
            LineContent::Idle => format!("{label} | idle"),
            LineContent::Active {
                chromosome,
                phase,
                processed,
                total,
            } => {
                let percent = if *total == 0 {
                    0.0
                } else {
                    (*processed as f64 / *total as f64) * 100.0
                };
                format!(
                    "{label} | {:<20} | {:>6.2}% | {}",
                    chromosome,
                    percent.min(100.0),
                    phase_label(*phase)
                )
            }
        }
    }
}

#[derive(Clone)]
enum LineContent {
    Idle,
    Active {
        chromosome: String,
        phase: ProgressPhase,
        processed: usize,
        total: usize,
    },
}

fn phase_label(phase: ProgressPhase) -> &'static str {
    match phase {
        ProgressPhase::Seeding => "seeding",
        ProgressPhase::Expanding => "expanding",
        ProgressPhase::Filtering => "filtering",
    }
}

fn current_thread_id() -> std::thread::ThreadId {
    std::thread::current().id()
}
