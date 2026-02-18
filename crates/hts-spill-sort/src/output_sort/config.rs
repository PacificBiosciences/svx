use crate::{Result, SpillSortError};
use std::path::PathBuf;

#[derive(Clone, Debug)]
pub struct SortConfig {
    pub max_mem: usize,
    pub tmp_dir: Option<PathBuf>,
    pub max_open_files: usize,
    pub merge_fan_in: usize,
}

impl SortConfig {
    pub fn new(
        max_mem: usize,
        tmp_dir: Option<PathBuf>,
        max_open_files: usize,
        merge_fan_in: usize,
    ) -> Result<Self> {
        if max_mem == 0 {
            return Err(SpillSortError::message(
                "sort-max-mem must be >= 1".to_string(),
            ));
        }
        if max_open_files == 0 {
            return Err(SpillSortError::message(
                "sort-max-open-files must be >= 1".to_string(),
            ));
        }
        if merge_fan_in == 0 {
            return Err(SpillSortError::message(
                "sort-merge-fan-in must be >= 1".to_string(),
            ));
        }
        if merge_fan_in > max_open_files {
            return Err(SpillSortError::message(format!(
                "sort-merge-fan-in ({merge_fan_in}) must be <= sort-max-open-files ({max_open_files})"
            )));
        }

        Ok(Self {
            max_mem,
            tmp_dir,
            max_open_files,
            merge_fan_in,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::SortConfig;

    #[test]
    fn sort_config_new_rejects_zero_limits() {
        let max_mem_error =
            SortConfig::new(0, None, 8, 4).expect_err("zero max_mem should be rejected");
        assert!(max_mem_error.to_string().contains("sort-max-mem"));

        let max_open_files_error =
            SortConfig::new(1, None, 0, 1).expect_err("zero max_open_files should be rejected");
        assert!(
            max_open_files_error
                .to_string()
                .contains("sort-max-open-files")
        );

        let fan_in_error =
            SortConfig::new(1, None, 8, 0).expect_err("zero merge_fan_in should be rejected");
        assert!(fan_in_error.to_string().contains("sort-merge-fan-in"));
    }

    #[test]
    fn sort_config_new_rejects_fan_in_above_max_open_files() {
        let error = SortConfig::new(1, None, 4, 5)
            .expect_err("fan-in above max open files should be rejected");
        assert!(error.to_string().contains("sort-merge-fan-in"));
        assert!(error.to_string().contains("sort-max-open-files"));
    }

    #[test]
    fn sort_config_new_accepts_valid_values() {
        let config =
            SortConfig::new(1024, None, 8, 4).expect("valid sort config should be constructed");
        assert_eq!(config.max_mem, 1024);
        assert_eq!(config.max_open_files, 8);
        assert_eq!(config.merge_fan_in, 4);
        assert!(config.tmp_dir.is_none());
    }
}
