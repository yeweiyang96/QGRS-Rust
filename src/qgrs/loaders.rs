use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;
use std::sync::Arc;

use memmap2::MmapOptions;

use crate::qgrs::data::{ChromSequence, InputMode};

use super::input::{is_gzip_path, open_input_reader};

pub fn load_sequences_from_path(path: &Path, mode: InputMode) -> io::Result<Vec<ChromSequence>> {
    match mode {
        InputMode::Mmap => load_sequences_mmap(path),
        InputMode::Stream => load_sequences_stream(path),
    }
}

fn load_sequences_stream(path: &Path) -> io::Result<Vec<ChromSequence>> {
    let mut reader = open_input_reader(path)?;
    parse_sequences_from_reader(reader.as_mut())
}

fn parse_sequences_from_reader(reader: &mut dyn BufRead) -> io::Result<Vec<ChromSequence>> {
    let mut sequences = Vec::new();
    let mut current_name: Option<String> = None;
    let mut sequence: Vec<u8> = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
            current_name = Some(parse_chrom_name(&line, sequences.len() + 1));
            continue;
        }
        for byte in line.bytes() {
            if byte.is_ascii_whitespace() {
                continue;
            }
            sequence.push(byte.to_ascii_lowercase());
        }
    }
    finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
    if !sequence.is_empty() {
        sequences.push(ChromSequence {
            name: format!("chromosome_{}", sequences.len() + 1),
            sequence: Arc::new(std::mem::take(&mut sequence)),
        });
    }
    Ok(sequences)
}

fn load_sequences_mmap(path: &Path) -> io::Result<Vec<ChromSequence>> {
    if is_gzip_path(path)? {
        let mut reader = open_input_reader(path)?;
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed)?;
        return Ok(parse_sequences_from_bytes(&decompressed));
    }
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    Ok(parse_sequences_from_bytes(&mmap))
}

fn parse_sequences_from_bytes(bytes: &[u8]) -> Vec<ChromSequence> {
    let mut sequences = Vec::new();
    let mut sequence = Vec::with_capacity(bytes.len());
    let mut current_name: Option<String> = None;
    let mut at_line_start = true;
    let mut i = 0;
    while i < bytes.len() {
        let byte = bytes[i];
        if byte == b'\n' || byte == b'\r' {
            at_line_start = true;
            i += 1;
            continue;
        }
        if at_line_start && byte == b'>' {
            finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
            i += 1;
            let header_start = i;
            while i < bytes.len() && bytes[i] != b'\n' && bytes[i] != b'\r' {
                i += 1;
            }
            let header = &bytes[header_start..i];
            current_name = Some(parse_chrom_name_bytes(header, sequences.len() + 1));
            at_line_start = true;
            continue;
        }
        at_line_start = false;
        if byte.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        sequence.push(byte.to_ascii_lowercase());
        i += 1;
    }
    finalize_sequence(&mut current_name, &mut sequence, &mut sequences);
    if !sequence.is_empty() {
        let fallback =
            current_name.unwrap_or_else(|| format!("chromosome_{}", sequences.len() + 1));
        sequences.push(ChromSequence {
            name: fallback,
            sequence: Arc::new(std::mem::take(&mut sequence)),
        });
    }
    sequences
}

fn finalize_sequence(
    current_name: &mut Option<String>,
    sequence: &mut Vec<u8>,
    sequences: &mut Vec<ChromSequence>,
) {
    if let Some(name) = current_name.take()
        && !sequence.is_empty()
    {
        sequences.push(ChromSequence {
            name,
            sequence: Arc::new(std::mem::take(sequence)),
        });
    }
}

pub(crate) fn parse_chrom_name(line: &str, index: usize) -> String {
    let header = line.trim_start_matches('>');
    header
        .split_whitespace()
        .next()
        .map(|s| s.to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| format!("chromosome_{index}"))
}

pub(crate) fn parse_chrom_name_bytes(header: &[u8], index: usize) -> String {
    let header_str = std::str::from_utf8(header).unwrap_or("");
    parse_chrom_name(header_str, index)
}
