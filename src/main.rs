mod config;
mod scanner;

use std::io::{self, BufWriter, Write};

use config::{Config, InputMode};
use scanner::{scan_mmap, scan_streaming};

fn main() {
	if let Err(err) = run() {
		eprintln!("Error: {err}");
		std::process::exit(1);
	}
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
	let config = Config::from_args()?;
	let stdout = io::stdout();
	let mut writer = BufWriter::with_capacity(1 << 20, stdout.lock());

	let total_emitted = match config.mode {
		InputMode::Mmap => scan_mmap(&config.input_path, config.min_tetrads, &mut writer)?,
		InputMode::Stream => scan_streaming(&config.input_path, config.min_tetrads, &mut writer)?,
	};

	writer.flush()?;
	eprintln!(
		"Detected {total_emitted} substrings with >= {} consecutive Gs",
		config.min_tetrads
	);
	Ok(())
}

