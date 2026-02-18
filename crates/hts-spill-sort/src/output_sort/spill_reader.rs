use crate::Result;
use std::io::Read;

use super::spill_format::{
    SPILL_FRAME_HEADER_LEN, SpillFormatEntry, decode_entry_from_parts, decode_frame_payload_len,
};

pub struct SpillReader<R: Read> {
    source: R,
    frame_header: [u8; SPILL_FRAME_HEADER_LEN],
    payload: Vec<u8>,
}

impl<R: Read> SpillReader<R> {
    pub fn new(source: R) -> Result<Self> {
        Ok(Self {
            source,
            frame_header: [0_u8; SPILL_FRAME_HEADER_LEN],
            payload: Vec::new(),
        })
    }

    pub fn read_next(&mut self) -> Result<Option<SpillFormatEntry>> {
        if !read_header_or_eof(&mut self.source, &mut self.frame_header)? {
            return Ok(None);
        }

        let payload_len = decode_frame_payload_len(self.frame_header.as_slice())?;
        self.payload.resize(payload_len, 0);
        self.source
            .read_exact(self.payload.as_mut_slice())
            .map_err(|error| {
                crate::SpillSortError::message(format!(
                    "failed reading encoded spill bytes from source: {}",
                    error
                ))
            })?;

        let entry = decode_entry_from_parts(self.frame_header.as_slice(), self.payload.as_slice())?;
        Ok(Some(entry))
    }
}

fn read_header_or_eof<R: Read>(source: &mut R, header: &mut [u8]) -> Result<bool> {
    let mut bytes_read = 0usize;
    while bytes_read < header.len() {
        let count = source.read(&mut header[bytes_read..]).map_err(|error| {
            crate::SpillSortError::message(format!(
                "failed reading encoded spill bytes from source: {}",
                error
            ))
        })?;
        if count == 0 {
            if bytes_read == 0 {
                return Ok(false);
            }
            return Err(crate::SpillSortError::message(
                "failed reading encoded spill bytes from source: truncated spill frame header"
                    .to_string(),
            ));
        }
        bytes_read = bytes_read.saturating_add(count);
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output_sort::spill_format::{SpillFormatEntry, SpillSourceOrder, encode_entry};
    use std::io;
    use std::sync::Arc;

    fn sample_entry(rid: u32, pos: i64, blob_ordinal: u64, run_id: u64) -> SpillFormatEntry {
        SpillFormatEntry {
            rid,
            pos,
            packed_alleles: Arc::from(b"\x00\x00\x00\x01g\x00\x00\x00\x01t".as_slice()),
            source: SpillSourceOrder {
                blob_ordinal,
                record_ordinal: 0,
            },
            run_id,
            record_payload: format!("payload-{rid}-{pos}-{run_id}").into(),
        }
    }

    struct ErrorSource;

    impl Read for ErrorSource {
        fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("reader failure"))
        }
    }

    #[test]
    fn reader_roundtrip_decodes_written_entry() {
        let entry = sample_entry(2, 99, 7, 31);
        let encoded =
            encode_entry(&entry).expect("spill format should encode for reader roundtrip");

        let mut reader =
            SpillReader::new(std::io::Cursor::new(encoded)).expect("reader should initialize");
        let observed = reader
            .read_next()
            .expect("reader should decode first entry")
            .expect("reader should return one decoded entry");

        assert_eq!(observed, entry);
        assert!(
            reader
                .read_next()
                .expect("reader should report stream exhaustion")
                .is_none()
        );
    }

    #[test]
    fn reader_rejects_corrupt_entry_bytes() {
        let mut reader = SpillReader::new(std::io::Cursor::new(vec![0xde, 0xad, 0xbe, 0xef]))
            .expect("reader should initialize");
        let error = reader
            .read_next()
            .expect_err("reader should fail on corrupt spill entry payload");
        assert!(!error.to_string().is_empty());
    }

    #[test]
    fn reader_surfaces_source_io_errors() {
        let mut reader = SpillReader::new(ErrorSource).expect("reader should initialize");
        let error = reader
            .read_next()
            .expect_err("reader should fail when source read fails");
        assert!(
            error
                .to_string()
                .contains("failed reading encoded spill bytes")
        );
    }
}
