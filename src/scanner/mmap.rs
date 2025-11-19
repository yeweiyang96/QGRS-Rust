use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use memmap2::MmapOptions;

use super::{emit_substrings, flush_run, is_g, GRunScanner};

pub fn scan_mmap(
	path: &Path,
	min_tetrads: usize,
	writer: &mut dyn Write,
) -> Result<u64, Box<dyn std::error::Error>> {
	let file = File::open(path)?;
	let mmap = unsafe { MmapOptions::new().map(&file)? };
	let data = &mmap[..];

	let total = if needs_fasta_handling(data) {
		scan_mmap_fasta(data, min_tetrads, writer)?
	} else {
		scan_mmap_plain(data, min_tetrads, writer)?
	};

	Ok(total)
}

fn needs_fasta_handling(data: &[u8]) -> bool {
	data.iter().any(|&byte| byte == b'>' || byte.is_ascii_whitespace())
}

fn scan_mmap_plain(
	data: &[u8],
	min_tetrads: usize,
	writer: &mut dyn Write,
) -> io::Result<u64> {
	let mut total = 0u64;
	for (start, run_len) in GRunScanner::new(data, min_tetrads) {
		total += emit_substrings(writer, start, run_len, min_tetrads)?;
	}
	Ok(total)
}

fn scan_mmap_fasta(
	data: &[u8],
	min_tetrads: usize,
	writer: &mut dyn Write,
) -> io::Result<u64> {
	let mut total = 0u64;
	let mut global_index = 0usize;
	let mut run_start = 0usize;
	let mut run_len = 0usize;
	let mut at_line_start = true;
	let mut skipping_header = false;

	for &byte in data {
		match byte {
			b'\n' | b'\r' => {
				at_line_start = true;
				skipping_header = false;
				continue;
			}
			_ => {}
		}

		if at_line_start {
			if byte == b'>' {
				skipping_header = true;
				at_line_start = false;
				continue;
			}
			at_line_start = false;
		}

		if skipping_header {
			continue;
		}

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

	total += flush_run(writer, run_start, run_len, min_tetrads)?;
	Ok(total)
}
