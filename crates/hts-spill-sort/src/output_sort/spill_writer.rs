use crate::Result;
use std::io::Write;

use super::spill_format::{SpillFormatEntry, encode_entry_with_run_id_into};

pub struct SpillWriter<W: Write> {
    sink: W,
    encode_buffer: Vec<u8>,
}

impl<W: Write> SpillWriter<W> {
    pub fn new(sink: W) -> Result<Self> {
        Ok(Self {
            sink,
            encode_buffer: Vec::new(),
        })
    }

    pub fn write_entry(&mut self, entry: &SpillFormatEntry) -> Result<()> {
        self.write_entry_with_run_id(entry, entry.run_id)
    }

    pub fn write_entry_with_run_id(&mut self, entry: &SpillFormatEntry, run_id: u64) -> Result<()> {
        encode_entry_with_run_id_into(entry, run_id, &mut self.encode_buffer)?;
        self.sink
            .write_all(self.encode_buffer.as_slice())
            .map_err(|error| {
                crate::SpillSortError::message(format!(
                    "failed writing encoded spill entry to sink: {}",
                    error
                ))
            })?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<W> {
        self.sink.flush().map_err(|error| {
            crate::SpillSortError::message(format!("failed flushing spill sink: {}", error))
        })?;
        Ok(self.sink)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output_sort::spill_format::SpillSourceOrder;
    use crate::output_sort::spill_reader::SpillReader;
    use std::io;
    use std::sync::Arc;

    fn sample_entry(blob_ordinal: u64, record_ordinal: usize, run_id: u64) -> SpillFormatEntry {
        SpillFormatEntry {
            rid: 0,
            pos: 45,
            packed_alleles: Arc::from(b"\x00\x00\x00\x01a\x00\x00\x00\x01c".as_slice()),
            source: SpillSourceOrder {
                blob_ordinal,
                record_ordinal,
            },
            run_id,
            record_payload: format!("record-{blob_ordinal}-{record_ordinal}-{run_id}").into(),
        }
    }

    #[derive(Default)]
    struct AlwaysFailSink;

    impl Write for AlwaysFailSink {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("test sink write failure"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn writer_reader_roundtrip_preserves_source_order_and_run_id() {
        let input = vec![sample_entry(4, 0, 12), sample_entry(4, 1, 12)];
        let mut writer = SpillWriter::new(Vec::new()).expect("spill writer should initialize");
        for entry in &input {
            writer
                .write_entry(entry)
                .expect("spill writer should encode each entry");
        }

        let bytes = writer
            .finish()
            .expect("spill writer should finish and return sink");
        let mut reader = SpillReader::new(std::io::Cursor::new(bytes))
            .expect("spill reader should initialize from spill bytes");
        let mut observed = Vec::new();
        loop {
            let next = reader
                .read_next()
                .expect("spill reader should decode next spill entry");
            if let Some(entry) = next {
                observed.push(entry);
                continue;
            }
            break;
        }

        assert_eq!(observed, input);
    }

    #[test]
    fn writer_surfaces_sink_write_failures() {
        let mut writer = SpillWriter::new(AlwaysFailSink).expect("spill writer should initialize");
        let error = writer
            .write_entry(&sample_entry(9, 3, 44))
            .expect_err("write should fail when sink write fails");
        assert!(!error.to_string().is_empty());
    }

    #[test]
    fn writer_can_override_run_id_without_cloning_payload() {
        let entry = sample_entry(4, 0, 12);
        let mut overridden = entry.clone();
        overridden.run_id = 33;
        let mut writer = SpillWriter::new(Vec::new()).expect("spill writer should initialize");
        writer
            .write_entry(&overridden)
            .expect("spill writer should encode entry with overridden run_id");

        let bytes = writer
            .finish()
            .expect("spill writer should finish and return sink");
        let mut reader = SpillReader::new(std::io::Cursor::new(bytes))
            .expect("spill reader should initialize from spill bytes");
        let observed = reader
            .read_next()
            .expect("spill reader should decode next spill entry")
            .expect("spill entry should be present");

        assert_eq!(observed.run_id, 33);
        assert_eq!(observed.record_payload, entry.record_payload);
    }
}
