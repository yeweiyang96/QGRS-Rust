use std::hash::{Hash, Hasher};
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InputMode {
    Mmap,
    Stream,
}

#[derive(Clone, Debug)]
pub struct ChromSequence {
    pub(crate) name: String,
    pub(crate) sequence: Arc<Vec<u8>>,
}

impl ChromSequence {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn sequence(&self) -> Arc<Vec<u8>> {
        Arc::clone(&self.sequence)
    }

    pub fn into_parts(self) -> (String, Arc<Vec<u8>>) {
        (self.name, self.sequence)
    }

    pub fn as_uppercase_string(&self) -> String {
        let mut seq = unsafe { String::from_utf8_unchecked(self.sequence.as_ref().clone()) };
        seq.make_ascii_uppercase();
        seq
    }
}

pub const DEFAULT_MAX_G4_LENGTH: usize = 45;
pub const DEFAULT_MAX_G_RUN: usize = 10;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ScanLimits {
    pub max_g4_length: usize,
    pub max_g_run: usize,
}

impl ScanLimits {
    pub const fn new(max_g4_length: usize, max_g_run: usize) -> Self {
        Self {
            max_g4_length,
            max_g_run,
        }
    }
}

impl Default for ScanLimits {
    fn default() -> Self {
        ScanLimits::new(DEFAULT_MAX_G4_LENGTH, DEFAULT_MAX_G_RUN)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SequenceData {
    pub(crate) normalized: Arc<Vec<u8>>,
}

impl SequenceData {
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn new(sequence: &str) -> Self {
        let normalized = Arc::new(sequence.to_ascii_lowercase().into_bytes());
        Self { normalized }
    }

    pub(crate) fn from_bytes(normalized: Arc<Vec<u8>>) -> Self {
        Self { normalized }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SequenceSlice {
    normalized: Arc<Vec<u8>>,
    start: usize,
    length: usize,
}

impl SequenceSlice {
    pub(crate) fn new(normalized: Arc<Vec<u8>>, start: usize, length: usize) -> Self {
        Self {
            normalized,
            start,
            length,
        }
    }

    pub(crate) fn bytes(&self) -> &[u8] {
        let end = self.start + self.length;
        &self.normalized[self.start..end]
    }

    pub(crate) fn to_uppercase_string(&self) -> String {
        let mut sequence = unsafe { String::from_utf8_unchecked(self.bytes().to_vec()) };
        sequence.make_ascii_uppercase();
        sequence
    }
}

impl PartialEq for SequenceSlice {
    fn eq(&self, other: &Self) -> bool {
        self.length == other.length && self.bytes() == other.bytes()
    }
}

impl Eq for SequenceSlice {}

impl Hash for SequenceSlice {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.length.hash(state);
        state.write(self.bytes());
    }
}

#[inline(always)]
pub(crate) fn is_g(byte: u8) -> bool {
    byte == b'G' || byte == b'g'
}
