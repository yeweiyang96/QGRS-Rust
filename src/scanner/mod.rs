mod mmap;
mod streaming;

pub use mmap::scan_mmap;
pub use streaming::scan_streaming;

use std::io::{self, Write};

pub(crate) struct GRunScanner<'a> {
	data: &'a [u8],
	cursor: usize,
	min_tetrads: usize,
}

impl<'a> GRunScanner<'a> {
	pub(crate) fn new(data: &'a [u8], min_tetrads: usize) -> Self {
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
		use memchr::memchr2;

		let len = self.data.len();

		while self.cursor < len {
			let search_slice = &self.data[self.cursor..];
			let relative = memchr2(b'G', b'g', search_slice)?;
			let run_start = self.cursor + relative;
			let mut run_end = run_start;

			while run_end < len && is_g(unsafe { *self.data.get_unchecked(run_end) }) {
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

pub(crate) fn emit_substrings(
	writer: &mut dyn Write,
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

pub(crate) fn flush_run(
	writer: &mut dyn Write,
	run_start: usize,
	run_len: usize,
	min_tetrads: usize,
) -> io::Result<u64> {
	if run_len >= min_tetrads {
		emit_substrings(writer, run_start, run_len, min_tetrads)
	} else {
		Ok(0)
	}
}

#[inline(always)]
pub(crate) fn is_g(byte: u8) -> bool {
	byte == b'G' || byte == b'g'
}

#[cfg(test)]
mod tests {
	use super::*;
	use std::fs;
	use std::time::{SystemTime, UNIX_EPOCH};

	#[test]
	fn finds_runs_with_threshold() {
		let data = b"TTGGGGATGGGN";
		let runs: Vec<_> = GRunScanner::new(data, 3).collect();
		assert_eq!(runs, vec![(2, 4), (8, 3)]);
	}

	#[test]
	fn finds_lowercase_runs() {
		let data = b"ttggggatgggn";
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

	#[test]
	fn streaming_mode_reads_fasta() {
		let filename = format!(
			"find_repeat_g_test_{}.fa",
			SystemTime::now()
				.duration_since(UNIX_EPOCH)
				.unwrap()
				.as_nanos()
		);
		let path = std::env::temp_dir().join(&filename);
		fs::write(&path, b">chr1\nACGGGGTT\n>chr2\nGGGG\n").unwrap();
		let mut buffer = Vec::new();
		let emitted = super::scan_streaming(&path, 4, &mut buffer).unwrap();
		fs::remove_file(&path).unwrap();
		assert_eq!(emitted, 2);
		let rendered = String::from_utf8(buffer).unwrap();
		let lines: Vec<_> = rendered.trim().split('\n').collect();
		assert_eq!(lines, vec!["2\t4", "8\t4"]);
	}

	#[test]
	fn mmap_mode_handles_fasta_boundaries() {
		let filename = format!(
			"find_repeat_g_mmap_{}.fa",
			SystemTime::now()
				.duration_since(UNIX_EPOCH)
				.unwrap()
				.as_nanos()
		);
		let path = std::env::temp_dir().join(filename);
		let contents = b">chr1\nGGGG\nGG\n";
		fs::write(&path, contents).unwrap();
		let mut mmap_buffer = Vec::new();
		let mmap_total = super::scan_mmap(&path, 4, &mut mmap_buffer).unwrap();
		let mut stream_buffer = Vec::new();
		let stream_total = super::scan_streaming(&path, 4, &mut stream_buffer).unwrap();
		fs::remove_file(&path).unwrap();
		assert_eq!(stream_total, mmap_total);
		assert_eq!(stream_buffer, mmap_buffer);
	}
}
