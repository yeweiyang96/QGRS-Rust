pub mod stream;

mod chunks;
mod consolidation;
mod data;
mod export;
mod loaders;
mod search;
#[cfg(test)]
mod tests;

pub use chunks::{find_owned_bytes, find_owned_bytes_with_limits};
pub use consolidation::consolidate_g4s;
pub use data::{ChromSequence, DEFAULT_MAX_G_RUN, DEFAULT_MAX_G4_LENGTH, InputMode, ScanLimits};
pub use export::{
    ExportError, render_csv_results, render_family_ranges_csv, write_parquet_results,
};
pub use loaders::load_sequences_from_path;
pub use search::G4;

#[cfg(test)]
pub(crate) use chunks::find_with_sequence;
pub(crate) use chunks::{chunk_size_for_limits, compute_chunk_overlap, shift_g4};
pub(crate) use loaders::parse_chrom_name;
pub(crate) use search::find_raw_bytes_no_chunking;
