use std::env;
use std::path::PathBuf;

#[derive(Debug)]
pub struct Config {
	pub min_tetrads: usize,
	pub input_path: PathBuf,
	pub mode: InputMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputMode {
	Mmap,
	Stream,
}

impl Config {
	pub fn from_args() -> Result<Self, String> {
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

		let mode = match args.next().as_deref() {
			Some("--mmap") => InputMode::Mmap,
			Some("--stream") => InputMode::Stream,
			Some(flag) => {
				return Err(usage(&format!(
					"unknown mode flag '{flag}'. expected --mmap or --stream"
				)));
			}
			None => return Err(usage("missing mode flag (--mmap | --stream)")),
		};

		if args.next().is_some() {
			return Err(usage("too many arguments"));
		}

		Ok(Self {
			min_tetrads,
			input_path,
			mode,
		})
	}
}

fn usage(reason: &str) -> String {
	format!(
		"{reason}\nUsage: find-repeat-G <min_tetrads> <path-to-sequence> (--mmap | --stream)"
	)
}
