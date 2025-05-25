# Command-line interface of SVX

Commands:

- `merge`

Options:
- `-v, --verbose` Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
- `-h, --help` Print help
- `-V, --version` Print version


## Merge command-line

Options:

`--vcf <VCF>...`         VCF files to merge

`--vcf-list <VCF_LIST>`  File containing paths of VCF files to merge (one per line)

`-o`, `--output <FILE>`       Write output to a file [default: standard output]

`-@ <THREADS>`     Number of threads to use [default: 1]

Advanced:

`-O`, `--output-type <OUTPUT_TYPE>`
    Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF

`--print-header`
    Print only the merged header and exit

`--force-single`
    Run even if there is only one file on input

`--no-version`
    Do not append version and command line to the header

`--contig <CONTIG>`
    Process only the specified contigs (comma-separated list), e.g., (chr1,chr2,chrX)

`--target-positions <TARGET-POSITIONS>`
    Specific positions or ranges to merge in format (contig:start[-end]) e.g., "(chr1:12345),(chr2:67890-67900)"

`--trs <TRS>`
    BED file with tandem repeat definitions

`--dump`
    Dump coordinates of variants

`--max-dist <MAX_DIST>`
    Maximum distance for merging variants [default: 1000]

`--min-dist <MIN_DIST>`
    The minimum distance threshold a variant can have when using max_dist_linear (-1 for no minimum) [default: 100]

`--min-sequence-similarity <MIN_SEQUENCE_SIMILARITY>`
    The minimum sequence identity for two insertions to be merged [default: 0.8]

`--max-dist-linear <MAX_DIST_LINEAR>`
    The proportion of the length of each variant to set distance threshold to [default: 0.5]

`--disable-linear-threshold`
    Disable distance threshold depending on variant length and use max_dist instead

`--kd-tree-norm <KD_TREE_NORM>`
    Norm to use for KD-tree distance calculations [default: 2]

`--knn-search-k <KNN_SEARCH_K>`
    Number of nearest neighbors to search in KD-tree [default: 4]

`--min-reciprocal-overlap <MIN_RECIP_OVERLAP>`
    The minimum reciprocal overlap for DEL/INV/DUP SVs [default: 0]

`--allow-intrasample`
    Allow merging variants from the same sample

`--require-mutual-distance`
    Require mutual distance for merging

`--filter-tr-contained`
    Filter out all SVs that are TR contained

`--tr-span-ratio-threshold <TR_SPAN_RATIO_THRESHOLD>`
    Maximum ratio of SV span to TR span for an SV to be considered TR-contained (SV_span / TR_span <= 1.0 + value) [default: 30]

`--tr-max-dist <TR_MAX_DIST>`
    Maximum distance for merging TR-contained variants (if they share the same TR ID) [default: 4000]

`--tr-min-sequence-similarity <TR_MIN_SEQUENCE_SIMILARITY>`
    Minimum sequence similarity for merging TR-contained variants (if they share the same TR ID) [default: 0.6]

`--tr-min-recip-overlap <TR_MIN_RECIP_OVERLAP>`
    Minimum reciprocal overlap for merging TR-contained DEL/INV/DUP SVs (if they share the same TR ID) [default: 0]

`--refinement-iterations <REFINEMENT_ITERATIONS>`
    Number of iterations for cluster refinement [default: 6]