# TR containment

TR containment is an annotation step that assigns an INS or DEL call to a tandem-repeat (TR) interval when the call is sufficiently consistent with that interval. The input TR set is provided via `--trs <bed>`. The expected BED schema follows the [TRGT repeat-catalog BED configurations](https://github.com/PacificBiosciences/trgt?tab=readme-ov-file#repeat-catalogs) with 4 columns (`chrom`, `start`, `end`, `info`), where the `info` field is semicolon-delimited key/value entries and must include at least `ID=<repeat_id>` and `MOTIFS=<motif1[,motif2,...]>`. Containment is evaluated against the BED intervals as genomic coordinates; the result is a best-scoring TR assignment (or no assignment) per variant. The annotation is currently defined for INS and DEL only.

Containment is evaluated differently for deletions and insertions because the underlying representations differ. Deletions are evaluated against TR intervals using span overlap, while insertions are evaluated using breakpoint proximity and boundary-aware interior checks. In both cases, candidates are ranked and a single TR is selected using deterministic tie-breaking to ensure stable results.

When two variants are TR-contained and share the same TR identifier, TR-specific merge thresholds are applied, which allows repeat-local merging behavior to be tuned independently from the global thresholds.

Related documentation:

- [Usage guide](guide.md)
- [Tuning and defaults](tuning.md)
- [Full CLI reference](cli.md)

## Filtering TR-contained calls

If you want to exclude repeat-contained calls from your merge entirely, pass `--filter-tr-contained` together with `--trs`. Any variant that is assigned a TR identifier (currently `INS` and `DEL`) is dropped before clustering, so it cannot contribute to merges and does not appear in the output.

## TR unrolling for TR-matched insertions

Within tandem repeats, equivalent insertion alleles can be reported at slightly shifted breakpoints. That shift often shows up as a cyclic rotation (phase shift) of the inserted sequence rather than a biological difference. For insertion pairs that already share the same TR identifier, SVX therefore evaluates sequence similarity in two stages:

1. Direct alignment of the inserted sequences against `--tr-min-sequence-similarity`.
2. If direct alignment fails, TR-aware cyclic unrolled comparisons are attempted using the observed breakpoint offset and the repeat motif period, testing both signs to accommodate upstream and downstream anchoring drift.

This fallback is restricted to `INS` pairs that already share the same TR identifier. If no unrolled comparison meets `--tr-min-sequence-similarity`, the pair does not satisfy the TR-specific similarity requirement.

Example:

```text
Call A: POS=100, ALT=ACGTACGT
Call B: POS=101, ALT=CGTACGTA

Direct compare:
  ACGTACGT
  CGTACGTA

Observed breakpoint offset = -1
Motif period = 4
Tested phase shifts include ±1 and ±3

Rotate Call B by +1:
  CGTACGTA -> ACGTACGT

Re-compare:
  ACGTACGT
  ACGTACGT
  -> passes sequence similarity
```
