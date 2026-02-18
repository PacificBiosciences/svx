use super::{
    config::SortConfig,
    keys::{PackedSortKey, estimate_record_bytes},
    pending::{PendingBlobRun, SortEntry, SourceOrder, drain_pending_runs_into_sink},
    spill::{
        SpillRun, create_spill_writer, finish_spill_writer, merge_spill_runs_into_sink,
        remove_run_file,
    },
    spill_format::{SpillFormatEntry, SpillSourceOrder},
    spill_payload::{decode_record_payload_into, duplicate_header_view, encode_record_payload},
};
use crate::{Result, SpillSortError};
use rust_htslib::bcf;
use std::path::PathBuf;

const LOG_PREFIX: &str = "hts-spill-sort";

#[cfg(feature = "logging")]
macro_rules! sorter_debug {
    ($($arg:tt)*) => {
        log::debug!($($arg)*)
    };
}

#[cfg(not(feature = "logging"))]
macro_rules! sorter_debug {
    ($($arg:tt)*) => {{
        if false {
            let _ = format_args!($($arg)*);
        }
    }};
}

#[cfg(feature = "logging")]
macro_rules! sorter_warn {
    ($($arg:tt)*) => {
        log::warn!($($arg)*)
    };
}

#[cfg(not(feature = "logging"))]
macro_rules! sorter_warn {
    ($($arg:tt)*) => {{
        if false {
            let _ = format_args!($($arg)*);
        }
    }};
}

#[derive(Debug)]
pub struct OutputSorter {
    header: bcf::Header,
    config: SortConfig,
    pub temp_dir: tempfile::TempDir,
    pub pending_runs: Vec<PendingBlobRun>,
    pub pending_bytes: usize,
    pub spill_runs: Vec<SpillRun>,
    next_run_id: u64,
}

impl OutputSorter {
    pub fn new(header: bcf::Header, config: SortConfig) -> Result<Self> {
        let config = SortConfig::new(
            config.max_mem,
            config.tmp_dir.clone(),
            config.max_open_files,
            config.merge_fan_in,
        )?;

        let temp_dir = match config.tmp_dir.as_ref() {
            Some(root) => tempfile::Builder::new()
                .prefix("spill-sort-")
                .tempdir_in(root)
                .map_err(|error| {
                    SpillSortError::message(format!(
                        "Failed to create sort temp directory under {}: {}",
                        root.display(),
                        error
                    ))
                })?,
            None => tempfile::Builder::new()
                .prefix("spill-sort-")
                .tempdir()
                .map_err(|error| {
                    SpillSortError::message(format!(
                        "Failed to create sort temp directory: {}",
                        error
                    ))
                })?,
        };

        let configured_tmp_root = config.tmp_dir.as_ref().map_or_else(
            || "<system-temp>".to_string(),
            |path| path.display().to_string(),
        );
        sorter_debug!(
            "{LOG_PREFIX}: initialized sorter max_mem={} max_open_files={} merge_fan_in={} tmp_root={} temp_dir={}",
            config.max_mem,
            config.max_open_files,
            config.merge_fan_in,
            configured_tmp_root,
            temp_dir.path().display()
        );

        Ok(Self {
            header,
            config,
            temp_dir,
            pending_runs: Vec::new(),
            pending_bytes: 0,
            spill_runs: Vec::new(),
            next_run_id: 0,
        })
    }

    pub fn push_blob_sorted_run(
        &mut self,
        blob_ordinal: u64,
        records: Vec<bcf::Record>,
    ) -> Result<()> {
        if records.is_empty() {
            return Ok(());
        }

        let record_count = records.len();
        let mut entries = Vec::with_capacity(records.len());
        let mut previous_key: Option<PackedSortKey> = None;
        let mut blob_bytes = 0usize;

        for (record_ordinal, record) in records.into_iter().enumerate() {
            let key = PackedSortKey::from_record(&record)?;
            if matches!(previous_key.as_ref(), Some(previous) if key < *previous) {
                return Err(SpillSortError::message(format!(
                    "Blob {blob_ordinal} is not sorted by output comparator"
                )));
            }
            previous_key = Some(key.clone());

            let source = SourceOrder {
                blob_ordinal,
                record_ordinal,
            };
            let approx_bytes = estimate_record_bytes(&key, &record);
            blob_bytes = blob_bytes.saturating_add(approx_bytes);
            entries.push(SortEntry {
                key,
                source,
                record,
            });
        }

        self.pending_bytes = self.pending_bytes.saturating_add(blob_bytes);
        self.pending_runs
            .push(PendingBlobRun { entries, cursor: 0 });
        sorter_debug!(
            "{LOG_PREFIX}: buffered blob={} records={} approx_bytes={} pending_runs={} pending_bytes={}/{}",
            blob_ordinal,
            record_count,
            blob_bytes,
            self.pending_runs.len(),
            self.pending_bytes,
            self.config.max_mem
        );
        if self.pending_bytes > self.config.max_mem {
            sorter_debug!(
                "{LOG_PREFIX}: memory cap exceeded after blob={} (pending_bytes={}/{}), spilling pending runs",
                blob_ordinal,
                self.pending_bytes,
                self.config.max_mem
            );
            self.spill_pending_runs()?;
        }
        Ok(())
    }

    pub fn finish_with<F>(&mut self, mut emit_record: F) -> Result<()>
    where
        F: FnMut(&bcf::Record) -> Result<()>,
    {
        if self.spill_runs.is_empty() {
            sorter_debug!(
                "{LOG_PREFIX}: no spill runs detected; draining buffered records directly"
            );
            drain_pending_runs_into_sink(
                &mut self.pending_runs,
                &mut self.pending_bytes,
                |entry| emit_record(&entry.record),
            )?;
            return Ok(());
        }

        if !self.pending_runs.is_empty() {
            sorter_debug!(
                "{LOG_PREFIX}: spill runs already exist and pending buffers remain; spilling pending runs before final merge"
            );
            self.spill_pending_runs()?;
        }

        self.compact_spill_runs()?;
        sorter_debug!(
            "{LOG_PREFIX}: merging {} spill runs into final output",
            self.spill_runs.len()
        );
        self.merge_spills_with(emit_record)
    }

    pub fn spill_pending_runs(&mut self) -> Result<()> {
        if self.pending_runs.is_empty() {
            return Ok(());
        }

        sorter_debug!(
            "{LOG_PREFIX}: spilling {} pending runs (pending_bytes={})",
            self.pending_runs.len(),
            self.pending_bytes
        );
        let run = self.create_spill_run_from_pending()?;
        sorter_debug!(
            "{LOG_PREFIX}: created spill run id={} path={}",
            run.run_id,
            run.path.display()
        );
        self.spill_runs.push(run);
        self.spill_runs.sort_by_key(|candidate| candidate.run_id);
        sorter_debug!(
            "{LOG_PREFIX}: spill run count now {}",
            self.spill_runs.len()
        );
        Ok(())
    }

    fn create_spill_run_from_pending(&mut self) -> Result<SpillRun> {
        let run_id = self.next_run_id;
        self.next_run_id = self.next_run_id.saturating_add(1);
        let path = self.run_path(run_id);
        sorter_debug!(
            "{LOG_PREFIX}: materializing pending runs into spill run id={} path={}",
            run_id,
            path.display()
        );

        let mut run_writer = create_spill_writer(
            path.as_path(),
            "Failed to create spill run",
            "Failed to initialize spill writer",
        )?;

        drain_pending_runs_into_sink(&mut self.pending_runs, &mut self.pending_bytes, |entry| {
            let spill_entry = SpillFormatEntry {
                rid: entry.key.rid(),
                pos: entry.key.pos(),
                packed_alleles: entry.key.packed_alleles(),
                source: SpillSourceOrder {
                    blob_ordinal: entry.source.blob_ordinal,
                    record_ordinal: entry.source.record_ordinal,
                },
                run_id,
                record_payload: encode_record_payload(&entry.record)?,
            };
            run_writer.write_entry(&spill_entry).map_err(|error| {
                SpillSortError::message(format!(
                    "Failed to write spill run {}: {}",
                    path.display(),
                    error
                ))
            })
        })?;
        finish_spill_writer(run_writer, path.as_path(), "Failed to finalize spill run")?;

        Ok(SpillRun { path, run_id })
    }

    fn run_path(&self, run_id: u64) -> PathBuf {
        self.temp_dir.path().join(format!("run_{run_id:012}.spill"))
    }

    fn compact_spill_runs(&mut self) -> Result<()> {
        let fan_in = self.config.merge_fan_in;
        while self.spill_runs.len() > fan_in {
            let chunk_size = fan_in.min(self.spill_runs.len());
            let chunk = self.spill_runs.drain(0..chunk_size).collect::<Vec<_>>();
            let merged = self.merge_runs_to_new_spill(chunk)?;
            self.spill_runs.push(merged);
            self.spill_runs.sort_by_key(|candidate| candidate.run_id);
        }
        Ok(())
    }

    fn merge_runs_to_new_spill(&mut self, runs: Vec<SpillRun>) -> Result<SpillRun> {
        let run_id = self.next_run_id;
        self.next_run_id = self.next_run_id.saturating_add(1);
        let path = self.run_path(run_id);
        let input_run_ids = runs.iter().map(|run| run.run_id).collect::<Vec<_>>();
        sorter_debug!(
            "{LOG_PREFIX}: k-way merging runs {:?} into spill run id={} path={}",
            input_run_ids,
            run_id,
            path.display()
        );

        let mut writer = create_spill_writer(
            path.as_path(),
            "Failed to create merged spill run",
            "Failed to initialize merged spill writer",
        )?;

        merge_spill_runs_into_sink(runs.as_slice(), |entry| {
            writer.write_entry(entry).map_err(|error| {
                SpillSortError::message(format!(
                    "Failed to write merged spill run {}: {}",
                    path.display(),
                    error
                ))
            })
        })?;
        finish_spill_writer(
            writer,
            path.as_path(),
            "Failed to finalize merged spill run",
        )?;

        for run in runs {
            remove_run_file(run.path.as_path())?;
        }

        Ok(SpillRun { path, run_id })
    }

    fn merge_spills_with<F>(&mut self, mut emit_record: F) -> Result<()>
    where
        F: FnMut(&bcf::Record) -> Result<()>,
    {
        let runs = std::mem::take(&mut self.spill_runs);
        let decode_header = match duplicate_header_view(&self.header) {
            Ok(header) => header,
            Err(error) => {
                self.cleanup_spill_runs(runs.as_slice())?;
                return Err(error);
            }
        };
        let mut decoded_record = decode_header.empty_record();

        let merge_result = merge_spill_runs_into_sink(runs.as_slice(), |entry| {
            decode_record_payload_into(
                entry.record_payload.as_slice(),
                entry.rid,
                entry.pos,
                &mut decoded_record,
            )?;
            emit_record(&decoded_record)
        });
        let cleanup_result = self.cleanup_spill_runs(runs.as_slice());

        if let Err(error) = merge_result {
            if let Err(cleanup_error) = cleanup_result {
                sorter_warn!(
                    "{LOG_PREFIX}: cleanup after failed final merge also failed: {}",
                    cleanup_error
                );
            }
            return Err(error);
        }

        cleanup_result
    }

    fn cleanup_spill_runs(&self, runs: &[SpillRun]) -> Result<()> {
        let mut first_error: Option<SpillSortError> = None;
        for run in runs {
            match remove_run_file(run.path.as_path()) {
                Err(error) if first_error.is_none() => first_error = Some(error),
                _ => {}
            }
        }
        match first_error {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::LOG_PREFIX;

    #[test]
    fn log_prefix_matches_crate_identity() {
        assert_eq!(LOG_PREFIX, "hts-spill-sort");
    }

    #[test]
    fn logging_feature_enabled_by_default() {
        assert!(cfg!(feature = "logging"));
    }
}
