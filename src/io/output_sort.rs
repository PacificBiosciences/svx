use crate::{io::vcf_writer::VcfWriter, utils::util::Result};
use rust_htslib::bcf;

pub use hts_spill_sort::{SortConfig, ensure_records_sorted_for_output};

#[derive(Debug)]
pub struct OutputSorter {
    inner: hts_spill_sort::OutputSorter,
}

impl OutputSorter {
    pub fn new(header: bcf::Header, config: SortConfig) -> Result<Self> {
        let inner = hts_spill_sort::OutputSorter::new(header, config)?;
        Ok(Self { inner })
    }

    pub fn push_blob_sorted_run(
        &mut self,
        blob_ordinal: u64,
        records: Vec<bcf::Record>,
    ) -> Result<()> {
        self.inner.push_blob_sorted_run(blob_ordinal, records)?;
        Ok(())
    }

    pub fn finish_into(&mut self, writer: &mut VcfWriter) -> Result<()> {
        self.inner.finish_with(|record| {
            writer
                .writer
                .write(record)
                .map_err(hts_spill_sort::SpillSortError::from)
        })?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spill_sort_crate_api_bridge_smoke() {
        let config = SortConfig::new(1024, None, 8, 4)
            .expect("external spill-sort config should initialize");
        let description = format!("{config:?}");
        assert!(description.contains("max_mem: 1024"));
        assert!(description.contains("max_open_files: 8"));
        assert!(description.contains("merge_fan_in: 4"));
    }
}
