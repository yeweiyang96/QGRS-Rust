use std::fmt;
use std::io::Write;
use std::sync::Arc;

use arrow_array::{ArrayRef, Int32Array, RecordBatch, StringArray, UInt64Array};
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::errors::ParquetError;

use crate::qgrs::data::SequenceTopology;
use crate::qgrs::search::G4;

pub fn render_family_ranges_csv(ranges: &[(usize, usize)]) -> String {
    render_family_ranges_csv_with_projection(ranges, SequenceTopology::Linear, 0)
}

pub fn render_family_ranges_csv_with_projection(
    ranges: &[(usize, usize)],
    topology: SequenceTopology,
    sequence_len: usize,
) -> String {
    let mut out = String::from("family_index,start,end\n");
    for (index, (start, end)) in ranges.iter().enumerate() {
        out.push_str(&format!(
            "{},{},{}\n",
            index + 1,
            start,
            projected_end(*end, topology, sequence_len)
        ));
    }
    out
}

pub fn render_csv_results(g4s: &[G4]) -> String {
    render_csv_results_with_projection(g4s, SequenceTopology::Linear, 0)
}

pub fn render_csv_results_with_projection(
    g4s: &[G4],
    topology: SequenceTopology,
    sequence_len: usize,
) -> String {
    let mut out = String::from("start,end,length,tetrads,y1,y2,y3,gscore,sequence\n");
    for g4 in g4s {
        let sequence_field = escape_csv_field(g4.sequence());
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{}\n",
            g4.start,
            projected_end(g4.end, topology, sequence_len),
            g4.length,
            g4.tetrads,
            g4.y1,
            g4.y2,
            g4.y3,
            g4.gscore,
            sequence_field
        ));
    }
    out
}

fn escape_csv_field(value: &str) -> String {
    if value.is_empty() {
        return String::new();
    }
    let needs_quotes = value.contains([',', '"', '\n']);
    if !needs_quotes {
        return value.to_string();
    }
    let mut escaped = String::from('"');
    for ch in value.chars() {
        if ch == '"' {
            escaped.push_str("\"\"");
        } else {
            escaped.push(ch);
        }
    }
    escaped.push('"');
    escaped
}

#[derive(Debug)]
pub enum ExportError {
    Arrow(arrow_schema::ArrowError),
    Parquet(ParquetError),
}

impl From<arrow_schema::ArrowError> for ExportError {
    fn from(value: arrow_schema::ArrowError) -> Self {
        ExportError::Arrow(value)
    }
}

impl From<ParquetError> for ExportError {
    fn from(value: ParquetError) -> Self {
        ExportError::Parquet(value)
    }
}

impl fmt::Display for ExportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExportError::Arrow(err) => write!(f, "arrow error: {err}"),
            ExportError::Parquet(err) => write!(f, "parquet error: {err}"),
        }
    }
}

impl std::error::Error for ExportError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ExportError::Arrow(err) => Some(err),
            ExportError::Parquet(err) => Some(err),
        }
    }
}

pub fn write_parquet_results<W: Write + Send + 'static>(
    g4s: &[G4],
    writer: W,
) -> Result<(), ExportError> {
    write_parquet_results_with_projection(g4s, writer, SequenceTopology::Linear, 0)
}

pub fn write_parquet_results_with_projection<W: Write + Send + 'static>(
    g4s: &[G4],
    writer: W,
    topology: SequenceTopology,
    sequence_len: usize,
) -> Result<(), ExportError> {
    write_parquet_from_results(g4s, writer, topology, sequence_len)
}

fn write_parquet_from_results<W: Write + Send + 'static>(
    g4s: &[G4],
    writer: W,
    topology: SequenceTopology,
    sequence_len: usize,
) -> Result<(), ExportError> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("start", DataType::UInt64, false),
        Field::new("end", DataType::UInt64, false),
        Field::new("length", DataType::UInt64, false),
        Field::new("tetrads", DataType::UInt64, false),
        Field::new("y1", DataType::Int32, false),
        Field::new("y2", DataType::Int32, false),
        Field::new("y3", DataType::Int32, false),
        Field::new("gscore", DataType::Int32, false),
        Field::new("sequence", DataType::Utf8, false),
    ]));

    let starts: Vec<u64> = g4s.iter().map(|g| g.start as u64).collect();
    let ends: Vec<u64> = g4s
        .iter()
        .map(|g| projected_end(g.end, topology, sequence_len) as u64)
        .collect();
    let lengths: Vec<u64> = g4s.iter().map(|g| g.length as u64).collect();
    let tetrads: Vec<u64> = g4s.iter().map(|g| g.tetrads as u64).collect();
    let y1s: Vec<i32> = g4s.iter().map(|g| g.y1).collect();
    let y2s: Vec<i32> = g4s.iter().map(|g| g.y2).collect();
    let y3s: Vec<i32> = g4s.iter().map(|g| g.y3).collect();
    let gscores: Vec<i32> = g4s.iter().map(|g| g.gscore).collect();
    let sequences: Vec<String> = g4s.iter().map(|g| g.sequence().to_string()).collect();

    let columns: Vec<ArrayRef> = vec![
        Arc::new(UInt64Array::from(starts)),
        Arc::new(UInt64Array::from(ends)),
        Arc::new(UInt64Array::from(lengths)),
        Arc::new(UInt64Array::from(tetrads)),
        Arc::new(Int32Array::from(y1s)),
        Arc::new(Int32Array::from(y2s)),
        Arc::new(Int32Array::from(y3s)),
        Arc::new(Int32Array::from(gscores)),
        Arc::new(StringArray::from(sequences)),
    ];

    let batch = RecordBatch::try_new(schema.clone(), columns)?;
    let mut arrow_writer = ArrowWriter::try_new(writer, schema, None)?;
    arrow_writer.write(&batch)?;
    arrow_writer.close()?;
    Ok(())
}

fn projected_end(end: usize, topology: SequenceTopology, sequence_len: usize) -> usize {
    if !topology.is_circular() || sequence_len == 0 || end == 0 {
        return end;
    }
    ((end - 1) % sequence_len) + 1
}
