# hts-spill-sort

This is a crate that provides deterministic global ordering for `rust-htslib` BCF/VCF records with bounded memory via spill-to-disk and k-way merge.

## Public API

- `SortConfig` for memory/open-file/fan-in limits.
- `OutputSorter` for staged ingestion and final emission.
- `ensure_records_sorted_for_output` for in-memory batch normalization.
- `SpillSortError` as the crate-local error type.

## Example usage

Typical usage is: buffer many already-sorted “runs” (`Vec<bcf::Record>`), then call `finish_with` to emit a single globally-sorted stream (spilling to disk once `max_mem` is exceeded).

```rust
use hts_spill_sort::{OutputSorter, Result, SortConfig, ensure_records_sorted_for_output};
use rust_htslib::bcf;

fn write_sorted_runs(
    mut out: bcf::Writer,
    runs: impl IntoIterator<Item = (u64, Vec<bcf::Record>)>,
) -> Result<()> {
    let config = SortConfig::new(
        256 * 1024 * 1024, // max_mem (bytes)
        None,              // tmp_dir (None => system temp)
        128,               // max_open_files
        64,                // merge_fan_in (must be <= max_open_files)
    )?;
    let mut sorter = OutputSorter::new(bcf::Header::from_template(out.header()), config)?;

    for (blob_ordinal, mut records) in runs {
        // Each run must be sorted by (RID, POS, alleles). This helper sorts stably if needed.
        ensure_records_sorted_for_output(&mut records)?;
        sorter.push_blob_sorted_run(blob_ordinal, records)?;
    }

    sorter.finish_with(|record| {
        out.write(record)?;
        Ok(())
    })
}
```

Notes:

- `blob_ordinal` is part of the stable tie-break behavior for equal keys; assign it deterministically (not “arrival order” from parallel workers).
- Records must have `rid`, `pos`, and alleles set (the sorter hard-errors if it can’t derive the sort key).

## Guarantees

- Deterministic ordering by `(RID, POS, case-sensitive allele lexical key)`.
- Deterministic tie behavior:
  - pending-run merge ties by source order (`blob_ordinal`, `record_ordinal`).
  - spill-run merge ties by spill entry `run_id`.
- Strict hard-error behavior on spill open/read/decode/write/remove failures.
- No silent fallback to unsorted output.

## Notes

- This crate is currently BCF/VCF focused (`rust-htslib::bcf::Record`).
- Logging can be compiled out with `--no-default-features` (default feature `logging` controls `log` emission).
