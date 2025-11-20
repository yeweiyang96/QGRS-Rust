use std::sync::Arc;

pub mod qgrs;
pub use qgrs::{DEFAULT_MAX_G_RUN, DEFAULT_MAX_G4_LENGTH, InputMode, ScanLimits};

pub type ProgressReporterHandle = Arc<dyn ProgressReporter + Send + Sync>;

pub trait ProgressReporter: Send + Sync {
    fn start_chromosome(&self, _name: &str, _length: usize) {}
    fn seed_progress(&self, _name: &str, _processed: usize, _total: usize) {}
    fn finish_chromosome(&self, _name: &str) {}
}

#[derive(Default)]
pub struct NoopProgressReporter;

impl ProgressReporter for NoopProgressReporter {}

pub fn noop_progress_handle() -> ProgressReporterHandle {
    Arc::new(NoopProgressReporter)
}
