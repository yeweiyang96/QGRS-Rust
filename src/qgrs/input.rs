use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use flate2::read::MultiGzDecoder;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
pub(crate) const INPUT_BUFFER_CAPACITY: usize = 1 << 20;

pub(crate) fn is_gzip_path(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let bytes_read = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    Ok(bytes_read == 2 && magic == GZIP_MAGIC)
}

pub(crate) fn open_input_reader(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    open_reader_from_file(file)
}

fn open_reader_from_file(mut file: File) -> io::Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let bytes_read = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if bytes_read == 2 && magic == GZIP_MAGIC {
        let reader = BufReader::with_capacity(INPUT_BUFFER_CAPACITY, MultiGzDecoder::new(file));
        Ok(Box::new(reader))
    } else {
        let reader = BufReader::with_capacity(INPUT_BUFFER_CAPACITY, file);
        Ok(Box::new(reader))
    }
}
