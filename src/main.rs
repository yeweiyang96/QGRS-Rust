use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use memchr::memchr;
use memmap2::MmapOptions;

fn main() {
	if let Err(err) = run() {
		eprintln!("Error: {err}");
		std::process::exit(1);
	}
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
	let config = Config::from_args()?;
	let file = File::open(&config.input_path)?;
	let mmap = unsafe { MmapOptions::new().map(&file)? };
	let data = &mmap[..];

	let stdout = io::stdout();
	let mut writer = BufWriter::with_capacity(1 << 20, stdout.lock());
	let mut total_emitted: u64 = 0;

	for (start, run_len) in GRunScanner::new(data, config.min_tetrads) {
		total_emitted += emit_substrings(&mut writer, start, run_len, config.min_tetrads)?;
	}

	writer.flush()?;
	eprintln!(
		"Detected {total_emitted} substrings with >= {} consecutive Gs",
		config.min_tetrads
	);
	Ok(())
}

#[derive(Debug)]
struct Config {
	min_tetrads: usize,
	input_path: PathBuf,
}

impl Config {
	fn from_args() -> Result<Self, String> {
		let mut args = env::args().skip(1);
		let min_tetrads = args
			.next()
			.ok_or_else(|| usage("missing <min_tetrads>"))?
			.parse::<usize>()
			.map_err(|_| usage("<min_tetrads> must be a positive integer"))?;

		if min_tetrads == 0 {
			return Err(usage("<min_tetrads> must be > 0"));
		}

		let input_path = args
			.next()
			.ok_or_else(|| usage("missing <path-to-sequence>"))?
			.into();

		if args.next().is_some() {
			return Err(usage("too many arguments"));
		}

		Ok(Self {
			min_tetrads,
			input_path,
		})
	}
}

fn usage(reason: &str) -> String {
	format!(
		"{reason}\nUsage: find-repeat-G <min_tetrads> <path-to-sequence>"
	)
}

struct GRunScanner<'a> {
	data: &'a [u8],
	cursor: usize,
	min_tetrads: usize,
}

impl<'a> GRunScanner<'a> {
	fn new(data: &'a [u8], min_tetrads: usize) -> Self {
		Self {
			data,
			cursor: 0,
			min_tetrads,
		}
	}
}

impl<'a> Iterator for GRunScanner<'a> {
	type Item = (usize, usize);

	fn next(&mut self) -> Option<Self::Item> {
		let len = self.data.len();

		while self.cursor < len {
			let search_slice = &self.data[self.cursor..];
			let relative = memchr(b'G', search_slice)?;
			let run_start = self.cursor + relative;
			let mut run_end = run_start;

			// SAFETY: run_end is bounded by the check in the loop condition.
			while run_end < len && unsafe { *self.data.get_unchecked(run_end) } == b'G' {
				run_end += 1;
			}

			self.cursor = if run_end < len { run_end + 1 } else { len };
			let run_len = run_end - run_start;

			if run_len >= self.min_tetrads {
				return Some((run_start, run_len));
			}
		}

		None
	}
}

fn emit_substrings<W: Write>(
	writer: &mut W,
	run_start: usize,
	run_len: usize,
	min_tetrads: usize,
) -> io::Result<u64> {
	let mut emitted = 0u64;
	let mut length = min_tetrads;

	while length <= run_len {
		let limit = run_len - length + 1;
		for offset in 0..limit {
			writeln!(writer, "{}\t{}", run_start + offset, length)?;
			emitted += 1;
		}
		length += 1;
	}

	Ok(emitted)
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn finds_runs_with_threshold() {
		let data = b"TTGGGGATGGGN";
		let runs: Vec<_> = GRunScanner::new(data, 3).collect();
		assert_eq!(runs, vec![(2, 4), (8, 3)]);
	}

	#[test]
	fn emits_all_substrings() {
		let mut buffer = Vec::new();
		let emitted = emit_substrings(&mut buffer, 5, 4, 3).unwrap();
		assert_eq!(emitted, 3);
		let buffer_str = String::from_utf8(buffer).unwrap();
		let lines: Vec<_> = buffer_str
			.trim()
			.split('\n')
			.map(|line| line.split('\t').collect::<Vec<_>>())
			.collect();
		assert_eq!(lines, vec![vec!["5", "3"], vec!["6", "3"], vec!["5", "4"]]);
	}
}
