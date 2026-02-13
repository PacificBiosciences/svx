# Output interpretation

This page explains how to read merged VCF records from `svx merge`.

Related docs:

- [Usage guide](guide.md)
- [Tuning and defaults](tuning.md)
- [TR containment](tr_containment.md)

## Common support fields

Merged records are written with these INFO fields:

- `SUPP`: number of samples supporting the merged record
- `SUPP_CALLS`: number of input calls merged into the record
- `SUPP_VEC`: variable-length integer support vector (`1` = supported, `0` = not supported) in merged sample order
- `IDLIST`: IDs of input records merged into the output record

`SUPP_VEC` is aligned to the merged VCF sample columns: the *i*th integer corresponds to the *i*th sample name in the output VCF header (the columns after `FORMAT`). Sample order is constructed by concatenating the sample lists from each input VCF in input order (`--vcf` argument order or `--vcf-list` line order), preserving the sample order within each input header.

## Non-BND records (`INS/DEL/INV/DUP/CNV`)

For non-BND merged groups, svx picks one representative input record for core VCF fields:

- `CHROM`, `POS`, `ID`, `REF`, `ALT`
- base `FORMAT` payload

Note, output `POS` is not the cluster mean; it comes from the representative input record. Additional merge-level INFO summaries include:

- `END`: merged end coordinate (`POS` + `abs(SVLEN)` for interval-like classes)
- `START_AVG`: mean breakpoint start (1-based report space, 2 decimals)
- `START_VARIANCE`: variance of breakpoint starts (2 decimals)
- `SVLEN_AVG`: mean `SVLEN` (2 decimals)
- `SVLEN_VARIANCE`: variance of `SVLEN` (2 decimals)
- `TR_CONTAINED`: TR ID when any merged variant in the group is TR-contained

For CNV output specifically:

- CNV class is derived from all merged members in the group:
  - DEL-only groups emit `SVTYPE=DEL` with `ALT=<DEL>`
  - DUP-only groups emit `SVTYPE=DUP` with `ALT=<DUP>`
  - mixed DEL+DUP groups (or groups containing explicit `SVTYPE=CNV`) emit `SVTYPE=CNV` with `ALT=<CNV>`
- Claim tags are aligned to that chosen output class:
  - when output class is DEL or DUP, `SVCLAIM` is emitted only if all merged records have the same claim; otherwise `SVCLAIM_SET` is emitted
  - when output class is CNV, svx emits `SVCLAIM_SET` (if any claims are present) and does not emit scalar `SVCLAIM`

## Breakends (BND) / translocations

Svx handles `SVTYPE=BND` as paired breakends representing one event. During merging, breakends are paired per input VCF using `INFO/MATEID`, merging happens at event level, not single-breakend level. Here orientation is part of event identity and missing mates are dropped (warning can be emitted depending on selection scope). In the output, each merged event is written as two `SVTYPE=BND` records with the same `INFO/EVENT` value, a reciprocal `INFO/MATEID`, IDs that are derived as `<event>_A` and `<event>_B` and the positions are rounded means of contributing event endpoints. Each BND output record also includes:

- `CHR2`: mate contig
- `STRANDS`: breakend orientation corresponding to ALT bracket form

Inversion-associated BND metadata is preserved conservatively. If all merged members of a BND group carry compatible inversion metadata, svx keeps inversion linkage (`INFO/EVENT`), outputting `INFO/EVENTTYPE=INV`, and marks both output breakends with `FILTER=InvBreakpoint`. If group members disagree or only a subset carries inversion metadata, svx falls back to standard BND output (synthesized `EVENT`, no `EVENTTYPE`, no `InvBreakpoint` filter).

## BND orientation from ALT

Svx writes orientation via VCF bracket notation in `ALT`. If you need `STRANDS`-like interpretation, infer as:

- `ALT` starts with `[` -> `STRANDS=+-`
- `ALT` starts with `]` -> `STRANDS=--`
- otherwise if `ALT` contains `[` -> `STRANDS=++`
- otherwise if `ALT` contains `]` -> `STRANDS=-+`
