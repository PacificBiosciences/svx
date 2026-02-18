use super::{
    keys::{PackedSortKey, packed_sort_key_from_spill_entry},
    spill_format::SpillFormatEntry,
    spill_reader::SpillReader,
    spill_writer::SpillWriter,
};
use crate::{Result, SpillSortError};
use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    fs,
    io::{BufReader, BufWriter},
    path::{Path, PathBuf},
};

#[derive(Clone, Debug)]
pub struct SpillRun {
    pub path: PathBuf,
    pub run_id: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SpillHeapItem {
    key: PackedSortKey,
    tie_run_id: u64,
    state_index: usize,
}

impl Ord for SpillHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .cmp(&self.key)
            .then_with(|| other.tie_run_id.cmp(&self.tie_run_id))
    }
}

impl PartialOrd for SpillHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct RunReadState {
    state_run_id: u64,
    reader: SpillReader<BufReader<fs::File>>,
    current_entry: SpillFormatEntry,
    current_key: PackedSortKey,
    current_tie_run_id: u64,
}

impl RunReadState {
    fn from_run(run: &SpillRun) -> Result<Option<Self>> {
        let spill_file = fs::File::open(&run.path).map_err(|error| {
            SpillSortError::message(format!(
                "Failed to open spill run {}: {}",
                run.path.display(),
                error
            ))
        })?;
        let mut reader = SpillReader::new(BufReader::new(spill_file))?;
        match reader.read_next() {
            Ok(Some(entry)) => {
                let key = packed_sort_key_from_spill_entry(&entry);
                Ok(Some(Self {
                    state_run_id: run.run_id,
                    reader,
                    current_tie_run_id: entry.run_id,
                    current_entry: entry,
                    current_key: key,
                }))
            }
            Ok(None) => Ok(None),
            Err(error) => Err(SpillSortError::message(format!(
                "Failed to read spill run {}: {}",
                run.path.display(),
                error
            ))),
        }
    }

    fn advance(&mut self) -> Result<bool> {
        match self.reader.read_next() {
            Ok(Some(entry)) => {
                self.current_key = packed_sort_key_from_spill_entry(&entry);
                self.current_tie_run_id = entry.run_id;
                self.current_entry = entry;
                Ok(true)
            }
            Ok(None) => Ok(false),
            Err(error) => Err(SpillSortError::message(format!(
                "Failed to read spill run record run_id={}: {}",
                self.state_run_id, error
            ))),
        }
    }
}

pub fn merge_spill_runs_into_sink<F>(runs: &[SpillRun], mut emit: F) -> Result<()>
where
    F: FnMut(&SpillFormatEntry) -> Result<()>,
{
    let mut states = Vec::new();
    for run in runs {
        if let Some(state) = RunReadState::from_run(run)? {
            states.push(state);
        }
    }

    let mut heap = BinaryHeap::new();
    for (state_index, state) in states.iter().enumerate() {
        heap.push(SpillHeapItem {
            key: state.current_key.clone(),
            tie_run_id: state.current_tie_run_id,
            state_index,
        });
    }

    while let Some(item) = heap.pop() {
        let state = &mut states[item.state_index];
        emit(&state.current_entry)?;
        if state.advance()? {
            heap.push(SpillHeapItem {
                key: state.current_key.clone(),
                tie_run_id: state.current_tie_run_id,
                state_index: item.state_index,
            });
        }
    }

    Ok(())
}

pub fn remove_run_file(path: &Path) -> Result<()> {
    fs::remove_file(path).map_err(|error| {
        SpillSortError::message(format!(
            "Failed to remove spill run {}: {}",
            path.display(),
            error
        ))
    })
}

pub fn create_spill_writer(
    path: &Path,
    create_error_message: &str,
    initialize_error_message: &str,
) -> Result<SpillWriter<BufWriter<fs::File>>> {
    let spill_file = fs::File::create(path).map_err(|error| {
        SpillSortError::message(format!(
            "{create_error_message} {}: {}",
            path.display(),
            error
        ))
    })?;
    SpillWriter::new(BufWriter::new(spill_file)).map_err(|error| {
        SpillSortError::message(format!(
            "{initialize_error_message} {}: {}",
            path.display(),
            error
        ))
    })
}

pub fn finish_spill_writer(
    writer: SpillWriter<BufWriter<fs::File>>,
    path: &Path,
    finalize_error_message: &str,
) -> Result<()> {
    writer.finish().map_err(|error| {
        SpillSortError::message(format!(
            "{finalize_error_message} {}: {}",
            path.display(),
            error
        ))
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::remove_run_file;

    #[test]
    fn remove_run_file_error_includes_cleanup_path_context() {
        let temp_dir =
            tempfile::TempDir::new().expect("temporary directory for cleanup test should exist");
        let missing_path = temp_dir.path().join("missing.spill");
        let error = remove_run_file(missing_path.as_path())
            .expect_err("removing a missing spill path should fail");
        let message = error.to_string();
        assert!(
            message.contains("Failed to remove spill run"),
            "cleanup error should include contextual prefix: {message}"
        );
        assert!(
            message.contains(missing_path.to_string_lossy().as_ref()),
            "cleanup error should include failing path: {message}"
        );
    }
}
