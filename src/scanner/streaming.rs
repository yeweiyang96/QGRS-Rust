use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use super::{flush_run, is_g};

pub fn scan_streaming(
	path: &Path,
	min_tetrads: usize,
	writer: &mut dyn Write,
) -> Result<u64, Box<dyn std::error::Error>> {
	let file = File::open(path)?;
	let mut reader = BufReader::with_capacity(1 << 20, file);
	let mut total = 0u64;
	let mut global_index = 0usize;
	let mut run_start = 0usize;
	let mut run_len = 0usize;
	let mut line = String::new();

	loop {
		line.clear();
		if reader.read_line(&mut line)? == 0 {
			break;
		}
		if line.starts_with('>') {
			total += flush_run(writer, run_start, run_len, min_tetrads)?;
			run_len = 0;
			continue;
		}
		for byte in line.bytes() {
			if byte.is_ascii_whitespace() {
				continue;
			}
			if is_g(byte) {
				if run_len == 0 {
					run_start = global_index;
				}
				run_len += 1;
			} else {
				total += flush_run(writer, run_start, run_len, min_tetrads)?;
				run_len = 0;
			}
			global_index += 1;
		}
	}

	total += flush_run(writer, run_start, run_len, min_tetrads)?;
	Ok(total)
}
