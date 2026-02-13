# Command-line interface of SVX

This file is generated from Clap help output.
Run `scripts/update_cli_docs.sh` after changing CLI flags/help text.

```text
Structural variant merger

Usage: svx merge [OPTIONS] <--vcf <VCF>...|--vcf-list <VCF_LIST>>

Options:
  -h, --help
          Print help (see a summary with '-h')

Inputs:
      --vcf <VCF>...
          VCF files to merge

      --vcf-list <VCF_LIST>
          File containing paths of VCF files to merge (one per line)

      --config <CONFIG>
          TOML config file with merge parameters (excluding input/output paths)

Output:
  -o, --output <FILE>
          Write output to a file [default: standard output]

  -O, --output-type <OUTPUT_TYPE>
          Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF

      --print-header
          Print only the merged header and exit

      --no-version
          Do not append version and command line to the header

Tandem repeats:
      --trs <TRS>
          BED file with tandem repeat definitions

      --filter-tr-contained
          Filter out all SVs that are TR contained

      --tr-span-query-slop <TR_SPAN_QUERY_SLOP>
          Expand DEL containment query interval by this many base pairs on each side
          
          [default: 0]

      --tr-min-span-containment <TR_MIN_SPAN_CONTAINMENT>
          Minimum DEL overlap fraction required for TR containment, scaled by 1,000,000 (e.g. 800000 = 0.8)
          
          [default: 800000]

      --tr-min-span-overlap-bp <TR_MIN_SPAN_OVERLAP_BP>
          Minimum overlap in base pairs required before DEL containment can pass
          
          [default: 1]

      --tr-ins-max-dist <TR_INS_MAX_DIST>
          Maximum distance (bp) for insertion TR containment proximity scoring
          
          [default: 0]

      --tr-max-dist <TR_MAX_DIST>
          Maximum distance for merging TR-contained variants (if they share the same TR ID)
          
          [default: 4000]

      --tr-min-sequence-similarity <TR_MIN_SEQUENCE_SIMILARITY>
          Minimum sequence similarity for merging TR-contained variants (if they share the same TR ID)
          
          [default: 0.6]

      --tr-min-recip-overlap <TR_MIN_RECIP_OVERLAP>
          Minimum reciprocal overlap for merging TR-contained DEL SVs (if they share the same TR ID)
          
          [default: 0]

Performance:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 1]

      --io-threads <IO_THREADS>
          Number of shared HTS I/O threads; applied to all input readers when >= 2
          
          [default: 2]

      --blob-queue-capacity <BLOB_QUEUE_CAPACITY>
          Capacity of reader -> worker queue (1..=65536, default derived from --threads)

      --result-queue-capacity <RESULT_QUEUE_CAPACITY>
          Capacity of worker -> writer queue (1..=65536, default derived from --threads)

      --no-shard
          Disable within-contig independent-shard parallelization

      --min-shard-size <MIN_SHARD_SIZE>
          Combine adjacent independent shards until each coalesced shard has at least this many variants (0 disables coalescing)
          
          [default: 0]

Diagnostics:
      --progress
          Force progress output (still requires interactive stderr and non-debug logging)

      --no-progress
          Disable progress output

      --dump-path <DUMP_PATH>
          Write merge diagnostics (KD-tree coordinates and group stats) to a TSV file

Selection:
      --force-single
          Run even if there is only one file on input

      --contig <CONTIG>
          Process only the specified contigs (comma-separated list), e.g., (chr1,chr2,chrX)

      --target-positions <TARGET-POSITIONS>
          Specific positions to merge in format (contig:start[-end]) e.g., (chr1:12345),(chr2:67890-67900)

      --svtype <SVTYPE>
          Restrict processing to specific SV types (comma-separated): INS, DEL, INV, DUP, BND, CNV, ALL (ALL excludes CNV)
          
          [default: ALL]
          [possible values: INS, DEL, INV, DUP, BND, CNV, ALL]

Distance thresholds:
      --max-dist <MAX_DIST>
          Maximum distance for merging variants
          
          [default: 1000]

      --min-dist <MIN_DIST>
          The minimum distance threshold a variant can have when using max_dist_linear (-1 for no minimum)
          
          [default: 100]

      --min-sequence-similarity <MIN_SEQUENCE_SIMILARITY>
          The minimum sequence identity for two insertions to be merged
          
          [default: 0.8]

      --min-size-similarity <MIN_SIZE_SIMILARITY>
          The minimum size similarity for two variants to be merged
          
          [default: 0]

      --max-dist-linear <MAX_DIST_LINEAR>
          The proportion of the length of each variant to set distance threshold to
          
          [default: 0.5]

      --disable-linear-threshold
          Disable distance threshold depending on variant length and use max_dist instead

      --min-reciprocal-overlap <MIN_RECIP_OVERLAP>
          The minimum reciprocal overlap for DEL/INV/DUP SVs
          
          [default: 0]

Clustering:
      --knn-search-k <KNN_SEARCH_K>
          Number of nearest neighbors to search in KD-tree
          
          [default: 4]

      --allow-intrasample
          Allow merging variants from the same sample

      --no-mutual-distance
          Disable mutual distance for merging

      --merge-constraint <MERGE_CONSTRAINT>
          Additional constraints to apply when forming merged clusters

          Possible values:
          - none:          Default SVX behavior (single-linkage-style connected components)
          - bbox_diameter: Prevent “bridging” merges by bounding overall cluster span in (start,end) space
          - clique:        Require every pair of variants in a cluster to be distance-mergeable
          - centroid:      Require every variant in a cluster to be within its own max distance of the cluster centroid
          
          [default: none]

Global:
      --color <COLOR>
          Enable or disable color output in logging
          
          [default: auto]
          [possible values: always, auto, never]

  -v, --verbose...
          Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
```
