# Changelog

## 0.6.0

### Highlights

- Added multi-sample VCF support for merge inputs
  - Enables incremental merging workflows and merging joint-called multi-sample sawfish output
  - Multi-sample records are treated as atomic records during merging and are not split into per-sample records
  - For better consistency, prefer using one-shot merging over incremental merging to reduce precision drift and preserve full variant context during merging
- Added optional global output sorting with `--sort-output`
  - Added automatic indexing for sorted compressed outputs (`.tbi` for VCF, `.csi` for BCF)
  - Sorting failures are hard errors (no fallback to unsorted output), and `--sort-tmp-dir` must exist (if unspecified will use system defaults)
- Added output support controls `--min-supp` and `--keep-monomorphic`
- Improved performance and internals across sharded merging, parse/write caching, lazy initialization, and shard/sorter refactors

### Important behavior changes

- CNVs are now always included (part of `--svtype ALL`)
  
### Output updates

- Legacy aggregate INFO stats (`START_AVG`, `START_VARIANCE`, `SVLEN_AVG`, `SVLEN_VARIANCE`, `SUPP_VEC`) are no longer written
  - Output uncertainty fields now use `CIPOS` / `CIEND`

## 0.5.0

**Major** update and refactor of svx with improved clustering and better parametrization, performance and reliability optimizations, added BND and CNV support, improved TR logic, and expanded documentation.

### Highlights

- Added BND and CNV (opt-in required) support in merge and output paths
- Reworked TR containment and matching
- Added configurable merge constraints
- Added size similarity as a new merging predicate
- Added TOML config support with the `--config` flag
- Added queue controls and shard-aware merging parallelizing within contigs for better scalability
- Added progress reporting and queue telemetry
- Added shared HTS I/O threadpool support for readers and writer

### Important behavior changes

- Input VCFs must contain exactly one sample each
- Sample names must be unique across all input VCFs
- CNV support is opt-in currently, `--svtype ALL` excludes `CNV` (use `--svtype ALL,CNV` or `--svtype CNV`).
- `--target-positions` now strictly rejects non-positive coordinates
- TR containment is now INS/DEL-only and `TR_CONTAINED` is emitted only when a TR ID is present

### Documentation

- Added and expanded documentation for usage, tuning/defaults, output interpretation, TR containment, and merge constraints

## 0.1.0

- Initial release of svx.
