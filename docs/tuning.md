# Tuning and defaults

This page describes the merge thresholds that most strongly affect candidate generation, merge decisions, and runtime.

Related documentation:

- [Usage guide](guide.md)
- [Merge constraints](merge_constraints.md)
- [TR containment](tr_containment.md)
- [Full CLI reference](cli.md)

## Current defaults

The current defaults are listed in the table below.

| Setting | Default | Meaning |
| --- | ---: | --- |
| `--svtype` | `ALL` | Processes `INS,DEL,INV,DUP,BND,CNV` by default. |
| `-t, --threads` | `1` | Number of merge worker threads. |
| `--io-threads` | `2` | Shared HTS I/O threads for readers/writer. |
| output destination | stdout | If `-o/--output` is omitted, output goes to stdout. |
| output format | uncompressed VCF | Default stream format when writing to stdout. |
| `--sort-output` | `false` | Keep input/runtime emission order unless global output sorting is explicitly enabled. |
| `--sort-max-mem` | `768M` | Memory cap for `--sort-output` buffering before spill-to-disk. |
| `--sort-tmp-dir` | unset | Spill files for `--sort-output` are placed under the system temp root unless set. |
| `--sort-max-open-files` (hidden) | `64` | Advanced spill-sort tuning: maximum number of spill runs kept open before partial merges are forced. |
| `--sort-merge-fan-in` (hidden) | `32` | Advanced spill-sort tuning for k-way merge fan-in; must be `<= --sort-max-open-files`. |
| `--min-supp` | `1` | Minimum supporting-sample count (`SUPP`) required to write a merged record. |
| `--keep-monomorphic` | `false` | Do not write `SUPP=0` merged records unless explicitly enabled. |
| `--blob-queue-capacity` | derived | Default is `clamp(2 * threads, 4, 64)` unless overridden. |
| `--result-queue-capacity` | derived | Default is `clamp(4 * threads, 8, 128)` unless overridden. |
| `--max-dist-linear` | `0.5` | Linear factor used to compute per-variant merge radius from `abs(SVLEN)`. |
| `--min-dist` | `100` | Floor for linear threshold (`-1` disables floor). |
| linear threshold mode | enabled | Effective per-variant distance is `max(0.5 * abs(SVLEN), 100)`. |
| `--max-dist` | `1000` | Fixed distance if linear thresholding is disabled. |
| `--min-sequence-similarity` | `0.8` | Minimum insertion sequence identity for merging. |
| `--min-size-similarity` | `0.0` | Minimum `min(abs(SVLEN_i),abs(SVLEN_j))/max(abs(SVLEN_i),abs(SVLEN_j))` required for merging. |
| `--min-reciprocal-overlap` | `0.0` | No reciprocal-overlap requirement by default for DEL/INV/DUP. |
| `--knn-search-k` | `4` | Number of nearest neighbors used for candidate generation (valid range: `1..=1024`). |
| `--no-shard` | `false` | Keeps shard-based within-block parallelism enabled by default. |
| `--min-shard-size` | `0` | Coalesce adjacent shards until size threshold is reached (`0` disables coalescing). |
| `--allow-intrasample` | `false` | Prevents same-sample merges by default. |
| mutual distance mode | enabled | Uses `min(v1.max_dist, v2.max_dist)` for pair thresholds. |
| `--merge-constraint` | `none` | Default single-linkage-style component growth. |
| `--filter-tr-contained` | `false` | Keep TR-contained calls by default. |
| `--tr-span-query-slop` | `0` | Expand DEL containment query on each side by this many bp. |
| `--tr-min-span-containment` | `800000` | DEL containment fraction threshold (scaled by 1,000,000; `800000 = 0.8`). |
| `--tr-min-span-overlap-bp` | `1` | Minimum DEL overlap in bp before containment can pass. |
| `--tr-ins-max-dist` | `0` | Insertion distance threshold for TR containment annotation. |
| `--tr-max-dist` | `4000` | TR-specific max distance when both calls share the same TR ID. |
| `--tr-min-sequence-similarity` | `0.6` | TR-specific insertion sequence identity threshold. |
| `--tr-min-recip-overlap` | `0.0` | TR-specific reciprocal-overlap threshold for DEL. |

## Distance threshold details

SVX first computes a per-variant merge radius (`max_dist`), then a per-pair threshold.

### 1. Per-variant radius (default linear mode)

With defaults:

`max_dist = max(0.5 * abs(SVLEN), 100)`

- `0.5` comes from `--max-dist-linear`
- `100` comes from `--min-dist`

Examples:

- `SVLEN=60` -> `0.5 * 60 = 30` -> floor applies -> `max_dist=100`
- `SVLEN=800` -> `0.5 * 800 = 400` -> above floor -> `max_dist=400`

If `--disable-linear-threshold` is set, SVX uses fixed `--max-dist` instead.

### 2. Per-pair threshold (mutual distance)

Default behavior:

`D(i,j) = min(max_dist_i, max_dist_j)`

With `--no-mutual-distance`:

`D(i,j) = max(max_dist_i, max_dist_j)`

Example:

- Variant A: `max_dist=100`
- Variant B: `max_dist=1000`
- Pair distance: `180`

Then:

- Mutual distance (default): `D=100` -> no merge
- `--no-mutual-distance`: `D=1000` -> merge can pass distance gate

This is why mutual distance generally reduces bridge/hub over-merging.

## Tuning patterns

- General-purpose collapse-like behavior: start with defaults, then adjust one parameter at a time.
- Strict deduplication (avoid over-merging): e.g. `--disable-linear-threshold --max-dist 100 --min-sequence-similarity 0.95 --min-size-similarity 0.7`.
- Loose harmonization across callsets: e.g. `--disable-linear-threshold --max-dist 1000 --min-sequence-similarity 0.7`; optionally add `--no-mutual-distance`.
- CNV / imprecise interval behavior: e.g. `--disable-linear-threshold --max-dist 2000 --min-sequence-similarity 0.0 --min-reciprocal-overlap 0.5`.
- Repeat-heavy regions: provide `--trs <bed>` and tune TR containment + TR merge flags separately from global thresholds.

If bridging/chaining is a concern, evaluate `--merge-constraint` modes in [Merge constraints](merge_constraints.md).

## Practical tuning

1. Restrict scope first (`--contig` or `--target-positions`) for faster iterations.
2. Fix one objective per run (recall, precision, or runtime).
3. Change one or two flags at a time.
4. Compare support fields (`SUPP`, `SUPP_CALLS`, `IDLIST`) and representative records in output.
