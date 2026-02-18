pub mod error;
pub mod output_sort;

pub use error::{Result, SpillSortError};
pub use output_sort::{OutputSorter, SortConfig, ensure_records_sorted_for_output};
