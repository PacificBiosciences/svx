use crate::{
    constants::*,
    core::svtype::SvType,
    io::{
        positions_reader::{PositionTuple, parse_position_tuple},
        vcf_writer::OutputType,
    },
};
use anyhow::{Context, Result, anyhow};
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum};
use log::{Level, LevelFilter};
use owo_colors::{
    OwoColorize, Stream, Style,
    colors::{Blue, Green, Magenta, Red, Yellow},
};
use serde::Deserialize;
use std::{
    collections::HashSet,
    ffi::OsString,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

#[cfg(has_git_describe)]
pub const FULL_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "-", env!("VERGEN_GIT_DESCRIBE"));

#[cfg(not(has_git_describe))]
pub const FULL_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser, Debug)]
#[command(name="svx",
          author="Tom Mokveld <tmokveld@pacificbiosciences.com>", 
          version=FULL_VERSION,
          about="Structural variation merging tool",
          long_about = None,
          after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
          This program comes with ABSOLUTELY NO WARRANTY; it is intended for
          Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Enable or disable color output in logging
    #[arg(long, value_enum, default_value_t = Color::Auto, global = true, help_heading = "Global")]
    color: Color,

    /// Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
    #[arg(
        short = 'v',
        long = "verbose",
        action = ArgAction::Count,
        global = true,
        help_heading = "Global"
    )]
    pub verbosity: u8,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    #[clap(about = "Structural variant merger")]
    Merge(MergeArgs),
}

impl Command {
    pub fn name(&self) -> &'static str {
        match self {
            Command::Merge(_) => "merge",
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(group(
    ArgGroup::new("input")
        .required(true)
        .args(["vcfs", "vcf_list"]),
))]
#[command(arg_required_else_help(true))]
pub struct MergeArgs {
    /// VCF files to merge
    #[arg(
        long = "vcf",
        value_name = "VCF",
        num_args = 1..,
        value_parser = check_file_exists,
        help_heading = "Inputs"
    )]
    pub vcfs: Option<Vec<PathBuf>>,

    /// File containing paths of VCF files to merge (one per line)
    #[arg(
        long = "vcf-list",
        value_name = "VCF_LIST",
        value_parser = check_file_exists,
        help_heading = "Inputs"
    )]
    pub vcf_list: Option<PathBuf>,

    /// TOML config file with merge parameters (excluding input/output paths)
    #[arg(
        long = "config",
        value_name = "CONFIG",
        value_parser = check_file_exists,
        help_heading = "Inputs"
    )]
    pub config: Option<PathBuf>,

    /// Write output to a file [default: standard output]
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        value_parser = check_prefix_path,
        help_heading = "Output"
    )]
    pub output: Option<String>,

    /// BED file with tandem repeat definitions
    #[arg(
        long = "trs",
        value_name = "TRS",
        value_parser = check_file_exists,
        help_heading = "Tandem repeats",
    )]
    pub tr_bed_path: Option<PathBuf>,

    /// Number of threads to use
    #[arg(
        short = 't',
        long = "threads",
        value_name = "THREADS",
        default_value = "1",
        value_parser = threads_in_range,
        help_heading = "Performance"
    )]
    pub num_threads: usize,

    /// Number of shared HTS I/O threads; applied to all input readers when >= 2
    #[arg(
        long = "io-threads",
        value_name = "IO_THREADS",
        default_value = "2",
        value_parser = threads_in_range,
        help_heading = "Performance"
    )]
    pub num_io_threads: usize,

    /// Capacity of reader -> worker queue (1..=65536, default derived from --threads)
    #[arg(
        long = "blob-queue-capacity",
        value_name = "BLOB_QUEUE_CAPACITY",
        value_parser = queue_capacity_in_range,
        help_heading = "Performance"
    )]
    pub blob_queue_capacity: Option<usize>,

    /// Capacity of worker -> writer queue (1..=65536, default derived from --threads)
    #[arg(
        long = "result-queue-capacity",
        value_name = "RESULT_QUEUE_CAPACITY",
        value_parser = queue_capacity_in_range,
        help_heading = "Performance"
    )]
    pub result_queue_capacity: Option<usize>,

    /// Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "OUTPUT_TYPE",
        value_parser = merge_validate_output_type,
        help_heading = "Output"
    )]
    pub output_type: Option<OutputType>,

    /// Print only the merged header and exit
    #[arg(long = "print-header", help_heading = "Output")]
    pub print_header: bool,

    /// Sort output records globally by coordinate and allele key before writing
    #[arg(long = "sort-output", help_heading = "Output")]
    pub sort_output: bool,

    /// Maximum memory for sorted output buffering (supports k/m/g suffixes, decimal units)
    #[arg(
        long = "sort-max-mem",
        value_name = "SORT_MAX_MEM",
        default_value = "768M",
        value_parser = parse_sort_max_mem,
        help_heading = "Output"
    )]
    pub sort_max_mem: usize,

    /// Temporary directory root for sorted output spill runs (must already exist)
    #[arg(
        long = "sort-tmp-dir",
        value_name = "SORT_TMP_DIR",
        value_parser = check_existing_dir_path,
        help_heading = "Output"
    )]
    pub sort_tmp_dir: Option<PathBuf>,

    /// Maximum number of spill runs that may remain before partial merges are forced
    #[arg(
        long = "sort-max-open-files",
        value_name = "SORT_MAX_OPEN_FILES",
        default_value_t = DEFAULT_SORT_MAX_OPEN_FILES,
        value_parser = positive_usize,
        hide = true
    )]
    pub sort_max_open_files: usize,

    /// Fan-in used for partial/final k-way merge over spill runs
    #[arg(
        long = "sort-merge-fan-in",
        value_name = "SORT_MERGE_FAN_IN",
        default_value_t = DEFAULT_SORT_MERGE_FAN_IN,
        value_parser = positive_usize,
        hide = true
    )]
    pub sort_merge_fan_in: usize,

    /// Force progress output (still requires interactive stderr and non-debug logging)
    #[arg(
        long = "progress",
        help_heading = "Diagnostics",
        conflicts_with = "no_progress"
    )]
    pub progress: bool,

    /// Disable progress output
    #[arg(
        long = "no-progress",
        help_heading = "Diagnostics",
        conflicts_with = "progress"
    )]
    pub no_progress: bool,

    /// Run even if there is only one file on input
    #[arg(long = "force-single", help_heading = "Selection")]
    pub force_single: bool,

    /// Do not append version and command line to the header
    #[arg(long = "no-version", help_heading = "Output")]
    pub no_version: bool,

    /// Process only the specified contigs (comma-separated list), e.g., (chr1,chr2,chrX)
    #[arg(
        long = "contig",
        value_name = "CONTIG",
        value_delimiter = ',',
        help_heading = "Selection"
    )]
    pub contigs: Option<Vec<String>>,

    /// Specific positions to merge in format (contig:start[-end]) e.g., (chr1:12345),(chr2:67890-67900)
    #[arg(
        long = "target-positions",
        value_name = "TARGET-POSITIONS",
        value_delimiter = ',',
        value_parser = parse_position_tuple,
        help_heading = "Selection"
    )]
    pub target_positions: Option<Vec<PositionTuple>>,

    /// Restrict processing to specific SV types (comma-separated): INS, DEL, INV, DUP, BND, CNV, ALL
    #[arg(
        long = "svtype",
        value_name = "SVTYPE",
        value_enum,
        value_delimiter = ',',
        ignore_case = true,
        default_value = "ALL",
        help_heading = "Selection"
    )]
    pub svtypes: Vec<MergeSvType>,

    #[command(flatten)]
    pub merge_args: MergeArgsInner,
}

#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
pub enum MergeConstraint {
    /// Default SVX behavior (single-linkage-style connected components)
    #[value(name = "none")]
    None,

    /// Prevent “bridging” merges by bounding overall cluster span in (start,end) space
    #[value(name = "bbox_diameter", alias = "bbox-diameter")]
    BboxDiameter,

    /// Require every pair of variants in a cluster to be distance-mergeable
    #[value(name = "clique")]
    Clique,

    /// Require every variant in a cluster to be within its own max distance of the cluster centroid
    #[value(name = "centroid")]
    Centroid,
}

impl Default for MergeConstraint {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq, Hash)]
pub enum MergeSvType {
    #[value(name = "INS")]
    Ins,
    #[value(name = "DEL")]
    Del,
    #[value(name = "INV")]
    Inv,
    #[value(name = "DUP")]
    Dup,
    #[value(name = "BND")]
    Bnd,
    #[value(name = "CNV")]
    Cnv,
    #[value(name = "ALL")]
    All,
}

impl MergeSvType {
    const ORDERED_DEFAULT_ALL: [Self; 6] = [
        Self::Ins,
        Self::Del,
        Self::Inv,
        Self::Dup,
        Self::Bnd,
        Self::Cnv,
    ];
    const ORDERED_SELECTABLE_NON_ALL: [Self; 6] = [
        Self::Ins,
        Self::Del,
        Self::Inv,
        Self::Dup,
        Self::Bnd,
        Self::Cnv,
    ];

    fn as_svtype(self) -> Option<SvType> {
        match self {
            Self::Ins => Some(SvType::INSERTION),
            Self::Del => Some(SvType::DELETION),
            Self::Inv => Some(SvType::INVERSION),
            Self::Dup => Some(SvType::DUPLICATION),
            Self::Bnd => Some(SvType::BND),
            Self::Cnv => Some(SvType::CNV),
            Self::All => None,
        }
    }
}

#[derive(Parser, Debug, Clone)]
pub struct MergeArgsInner {
    /// Write merge diagnostics (KD-tree coordinates and group stats) to a TSV file
    #[arg(
        help_heading = "Diagnostics",
        long = "dump-path",
        value_name = "DUMP_PATH",
        value_parser = check_prefix_pathbuf
    )]
    pub dump_path: Option<PathBuf>,

    /// Maximum distance for merging variants
    #[arg(
        help_heading = "Distance thresholds",
        long,
        default_value_t = DEFAULT_MAX_DIST,
        value_parser = max_dist_in_range
    )]
    pub max_dist: i32,

    /// The minimum distance threshold a variant can have when using max_dist_linear (-1 for no minimum)
    #[arg(
        help_heading = "Distance thresholds",
        long,
        default_value_t = DEFAULT_MIN_DIST,
        value_parser = min_dist_in_range
    )]
    pub min_dist: i32,

    /// The minimum sequence identity for two insertions to be merged
    #[arg(
        help_heading = "Distance thresholds",
        long,
        default_value_t = DEFAULT_MIN_SEQUENCE_SIMILARITY
    )]
    pub min_sequence_similarity: f32,

    /// The minimum size similarity for two variants to be merged
    #[arg(
        help_heading = "Distance thresholds",
        long,
        default_value_t = DEFAULT_MIN_SIZE_SIMILARITY
    )]
    pub min_size_similarity: f64,

    /// The proportion of the length of each variant to set distance threshold to
    #[arg(
        help_heading = "Distance thresholds",
        long,
        default_value_t = DEFAULT_MAX_DIST_LINEAR
    )]
    pub max_dist_linear: f32,

    /// Disable distance threshold depending on variant length and use max_dist instead
    #[arg(
        help_heading = "Distance thresholds",
        long = "disable-linear-threshold",
        action = ArgAction::SetFalse,
        default_value_t = DEFAULT_USE_LINEAR_THRESHOLD
    )]
    pub use_linear_threshold: bool,

    /// Number of nearest neighbors to search in KD-tree
    #[arg(
        help_heading = "Clustering",
        long,
        default_value_t = DEFAULT_KNN_SEARCH_K,
        value_parser = knn_search_k_in_range
    )]
    pub knn_search_k: usize,

    /// Disable within-contig independent-shard parallelization
    #[arg(help_heading = "Performance", long = "no-shard")]
    pub no_shard: bool,

    /// Combine adjacent independent shards until each coalesced shard has at least this many variants (0 disables coalescing)
    #[arg(
        help_heading = "Performance",
        long = "min-shard-size",
        default_value_t = DEFAULT_MIN_SHARD_SIZE
    )]
    pub min_shard_size: usize,

    /// The minimum reciprocal overlap for DEL/INV/DUP SVs
    #[arg(
        help_heading = "Distance thresholds",
        long = "min-reciprocal-overlap",
        default_value_t = DEFAULT_MIN_RECIP_OVERLAP
    )]
    pub min_recip_overlap: f32,

    /// Allow merging variants from the same sample
    #[arg(
        help_heading = "Clustering",
        long,
        default_value_t = DEFAULT_ALLOW_INTRASAMPLE
    )]
    pub allow_intrasample: bool,

    /// Minimum number of supporting samples (SUPP) required to write a merged record
    #[arg(
        help_heading = "Output",
        long = "min-supp",
        value_name = "MIN_SUPP",
        default_value_t = DEFAULT_MIN_SUPP,
        value_parser = positive_usize
    )]
    pub min_supp: usize,

    /// Allow writing monomorphic merged records with SUPP=0
    #[arg(
        help_heading = "Output",
        long = "keep-monomorphic",
        default_value_t = DEFAULT_KEEP_MONOMORPHIC
    )]
    pub keep_monomorphic: bool,

    /// Disable mutual distance for merging
    #[arg(
        help_heading = "Clustering",
        long = "no-mutual-distance",
        action = ArgAction::SetFalse,
        default_value_t = DEFAULT_REQUIRE_MUTUAL_DISTANCE
    )]
    pub require_mutual_distance: bool,

    /// Additional constraints to apply when forming merged clusters
    #[arg(
        help_heading = "Clustering",
        long = "merge-constraint",
        value_enum,
        default_value_t = MergeConstraint::None
    )]
    pub merge_constraint: MergeConstraint,

    /// Filter out all SVs that are TR contained
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_FILTER_TR_CONTAINED
    )]
    pub filter_tr_contained: bool,

    /// Expand DEL containment query interval by this many base pairs on each side
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_SPAN_QUERY_SLOP
    )]
    pub tr_span_query_slop: u32,

    /// Minimum DEL overlap fraction required for TR containment, scaled by 1,000,000 (e.g. 800000 = 0.8)
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_MIN_SPAN_CONTAINMENT
    )]
    pub tr_min_span_containment: u32,

    /// Minimum overlap in base pairs required before DEL containment can pass
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_MIN_SPAN_OVERLAP_BP
    )]
    pub tr_min_span_overlap_bp: u32,

    /// Maximum distance (bp) for insertion TR containment proximity scoring
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_INS_MAX_DIST
    )]
    pub tr_ins_max_dist: u32,

    /// Maximum distance for merging TR-contained variants (if they share the same TR ID).
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_MAX_DIST,
        value_parser = tr_max_dist_in_range
    )]
    pub tr_max_dist: i32,

    /// Minimum sequence similarity for merging TR-contained variants (if they share the same TR ID).
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_MIN_SEQUENCE_SIMILARITY
    )]
    pub tr_min_sequence_similarity: f32,

    /// Minimum reciprocal overlap for merging TR-contained DEL SVs (if they share the same TR ID).
    #[arg(
        help_heading = "Tandem repeats",
        long,
        default_value_t = DEFAULT_TR_MIN_RECIP_OVERLAP
    )]
    pub tr_min_recip_overlap: f32,
}

impl MergeArgsInner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with(mut self, updates: impl FnOnce(&mut Self)) -> Self {
        updates(&mut self);
        self
    }
}

impl Default for MergeArgsInner {
    fn default() -> Self {
        Self {
            dump_path: None,
            max_dist: DEFAULT_MAX_DIST,
            min_dist: DEFAULT_MIN_DIST,
            max_dist_linear: DEFAULT_MAX_DIST_LINEAR,
            use_linear_threshold: DEFAULT_USE_LINEAR_THRESHOLD,
            min_sequence_similarity: DEFAULT_MIN_SEQUENCE_SIMILARITY,
            min_size_similarity: DEFAULT_MIN_SIZE_SIMILARITY,
            knn_search_k: DEFAULT_KNN_SEARCH_K,
            no_shard: false,
            min_shard_size: DEFAULT_MIN_SHARD_SIZE,
            allow_intrasample: DEFAULT_ALLOW_INTRASAMPLE,
            min_supp: DEFAULT_MIN_SUPP,
            keep_monomorphic: DEFAULT_KEEP_MONOMORPHIC,
            require_mutual_distance: DEFAULT_REQUIRE_MUTUAL_DISTANCE,
            merge_constraint: MergeConstraint::None,
            min_recip_overlap: DEFAULT_MIN_RECIP_OVERLAP,
            filter_tr_contained: DEFAULT_FILTER_TR_CONTAINED,
            tr_span_query_slop: DEFAULT_TR_SPAN_QUERY_SLOP,
            tr_min_span_containment: DEFAULT_TR_MIN_SPAN_CONTAINMENT,
            tr_min_span_overlap_bp: DEFAULT_TR_MIN_SPAN_OVERLAP_BP,
            tr_ins_max_dist: DEFAULT_TR_INS_MAX_DIST,
            tr_max_dist: DEFAULT_TR_MAX_DIST,
            tr_min_sequence_similarity: DEFAULT_TR_MIN_SEQUENCE_SIMILARITY,
            tr_min_recip_overlap: DEFAULT_TR_MIN_RECIP_OVERLAP,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct MergeTomlConfig {
    trs: Option<PathBuf>,
    threads: Option<usize>,
    io_threads: Option<usize>,
    blob_queue_capacity: Option<usize>,
    result_queue_capacity: Option<usize>,
    print_header: Option<bool>,
    sort_output: Option<bool>,
    sort_max_mem: Option<SortMaxMemConfigValue>,
    sort_tmp_dir: Option<PathBuf>,
    sort_max_open_files: Option<usize>,
    sort_merge_fan_in: Option<usize>,
    progress: Option<bool>,
    no_progress: Option<bool>,
    force_single: Option<bool>,
    no_version: Option<bool>,
    contig: Option<Vec<String>>,
    target_positions: Option<Vec<String>>,
    svtype: Option<Vec<String>>,
    dump_path: Option<PathBuf>,
    max_dist: Option<i32>,
    min_dist: Option<i32>,
    min_sequence_similarity: Option<f32>,
    min_size_similarity: Option<f64>,
    max_dist_linear: Option<f32>,
    disable_linear_threshold: Option<bool>,
    knn_search_k: Option<usize>,
    no_shard: Option<bool>,
    min_shard_size: Option<usize>,
    min_reciprocal_overlap: Option<f32>,
    allow_intrasample: Option<bool>,
    min_supp: Option<usize>,
    keep_monomorphic: Option<bool>,
    no_mutual_distance: Option<bool>,
    merge_constraint: Option<String>,
    filter_tr_contained: Option<bool>,
    tr_span_query_slop: Option<u32>,
    tr_min_span_containment: Option<u32>,
    tr_min_span_overlap_bp: Option<u32>,
    tr_ins_max_dist: Option<u32>,
    tr_max_dist: Option<i32>,
    tr_min_sequence_similarity: Option<f32>,
    tr_min_recip_overlap: Option<f32>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum SortMaxMemConfigValue {
    Integer(usize),
    Float(f64),
    Text(String),
}

impl SortMaxMemConfigValue {
    fn parse_bytes(&self) -> Result<usize> {
        match self {
            SortMaxMemConfigValue::Integer(value) => parse_sort_max_mem(&value.to_string()),
            SortMaxMemConfigValue::Float(value) => parse_sort_max_mem(&value.to_string()),
            SortMaxMemConfigValue::Text(value) => parse_sort_max_mem(value),
        }
    }
}

impl MergeTomlConfig {
    fn from_path(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;
        toml::from_str(&contents)
            .with_context(|| format!("Failed to parse TOML config file: {}", path.display()))
    }
}

fn resolve_config_relative_path(config_path: &Path, value: PathBuf) -> PathBuf {
    if value.is_absolute() {
        value
    } else {
        config_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(value)
    }
}

fn is_cli_value_set(merge_matches: &clap::ArgMatches, arg_id: &str) -> bool {
    merge_matches
        .value_source(arg_id)
        .is_some_and(|source| source == clap::parser::ValueSource::CommandLine)
}

fn parse_config_svtypes(values: &[String]) -> Result<Vec<MergeSvType>> {
    if values.is_empty() {
        return Err(anyhow!(
            "`svtype` in config must contain at least one value"
        ));
    }

    values
        .iter()
        .map(|value| {
            MergeSvType::from_str(value, true).map_err(|error| {
                anyhow!(
                    "Invalid `svtype` value `{}` in config: {}",
                    value,
                    error.replace('\n', " ")
                )
            })
        })
        .collect()
}

fn parse_config_merge_constraint(value: &str) -> Result<MergeConstraint> {
    MergeConstraint::from_str(value, true).map_err(|error| {
        anyhow!(
            "Invalid `merge_constraint` value `{}` in config: {}",
            value,
            error.replace('\n', " ")
        )
    })
}

fn parse_config_target_positions(values: &[String]) -> Result<Vec<PositionTuple>> {
    values
        .iter()
        .map(|value| {
            parse_position_tuple(value)
                .map_err(|error| anyhow!("Invalid target_positions entry `{value}`: {error}"))
        })
        .collect()
}

fn apply_merge_toml_config(
    merge_args: &mut MergeArgs,
    merge_matches: &clap::ArgMatches,
) -> Result<()> {
    let Some(config_path) = merge_args.config.as_deref() else {
        return Ok(());
    };

    let config = MergeTomlConfig::from_path(config_path)?;

    if let (false, Some(path)) = (is_cli_value_set(merge_matches, "tr_bed_path"), config.trs) {
        let resolved = resolve_config_relative_path(config_path, path);
        if !resolved.exists() {
            return Err(anyhow!(
                "File does not exist for `trs` in config: {}",
                resolved.display()
            ));
        }
        merge_args.tr_bed_path = Some(resolved);
    }

    if let (false, Some(threads)) = (
        is_cli_value_set(merge_matches, "num_threads"),
        config.threads,
    ) {
        merge_args.num_threads = threads_in_range(&threads.to_string())?;
    }

    if let (false, Some(io_threads)) = (
        is_cli_value_set(merge_matches, "num_io_threads"),
        config.io_threads,
    ) {
        merge_args.num_io_threads = threads_in_range(&io_threads.to_string())?;
    }

    if let (false, Some(capacity)) = (
        is_cli_value_set(merge_matches, "blob_queue_capacity"),
        config.blob_queue_capacity,
    ) {
        merge_args.blob_queue_capacity = Some(queue_capacity_in_range(&capacity.to_string())?);
    }

    if let (false, Some(capacity)) = (
        is_cli_value_set(merge_matches, "result_queue_capacity"),
        config.result_queue_capacity,
    ) {
        merge_args.result_queue_capacity = Some(queue_capacity_in_range(&capacity.to_string())?);
    }

    if let (false, Some(print_header)) = (
        is_cli_value_set(merge_matches, "print_header"),
        config.print_header,
    ) {
        merge_args.print_header = print_header;
    }

    if let (false, Some(sort_output)) = (
        is_cli_value_set(merge_matches, "sort_output"),
        config.sort_output,
    ) {
        merge_args.sort_output = sort_output;
    }

    if let (false, Some(sort_max_mem)) = (
        is_cli_value_set(merge_matches, "sort_max_mem"),
        config.sort_max_mem,
    ) {
        merge_args.sort_max_mem = sort_max_mem.parse_bytes()?;
    }

    if let (false, Some(sort_tmp_dir)) = (
        is_cli_value_set(merge_matches, "sort_tmp_dir"),
        config.sort_tmp_dir,
    ) {
        let resolved = resolve_config_relative_path(config_path, sort_tmp_dir);
        validate_existing_dir_path(&resolved)?;
        merge_args.sort_tmp_dir = Some(resolved);
    }

    if let (false, Some(sort_max_open_files)) = (
        is_cli_value_set(merge_matches, "sort_max_open_files"),
        config.sort_max_open_files,
    ) {
        merge_args.sort_max_open_files = positive_usize(&sort_max_open_files.to_string())?;
    }

    if let (false, Some(sort_merge_fan_in)) = (
        is_cli_value_set(merge_matches, "sort_merge_fan_in"),
        config.sort_merge_fan_in,
    ) {
        merge_args.sort_merge_fan_in = positive_usize(&sort_merge_fan_in.to_string())?;
    }

    if let (false, Some(progress)) = (is_cli_value_set(merge_matches, "progress"), config.progress)
    {
        merge_args.progress = progress;
    }

    if let (false, Some(no_progress)) = (
        is_cli_value_set(merge_matches, "no_progress"),
        config.no_progress,
    ) {
        merge_args.no_progress = no_progress;
    }

    if merge_args.progress && merge_args.no_progress {
        return Err(anyhow!(
            "Config conflict: `progress` and `no_progress` cannot both be true"
        ));
    }

    if let (false, Some(force_single)) = (
        is_cli_value_set(merge_matches, "force_single"),
        config.force_single,
    ) {
        merge_args.force_single = force_single;
    }

    if let (false, Some(no_version)) = (
        is_cli_value_set(merge_matches, "no_version"),
        config.no_version,
    ) {
        merge_args.no_version = no_version;
    }

    if let (false, Some(contigs)) = (is_cli_value_set(merge_matches, "contigs"), config.contig) {
        merge_args.contigs = Some(contigs);
    }

    if let (false, Some(target_positions)) = (
        is_cli_value_set(merge_matches, "target_positions"),
        config.target_positions,
    ) {
        merge_args.target_positions = Some(parse_config_target_positions(&target_positions)?);
    }

    if let (false, Some(svtypes)) = (is_cli_value_set(merge_matches, "svtypes"), config.svtype) {
        merge_args.svtypes = parse_config_svtypes(&svtypes)?;
    }

    if let (false, Some(dump_path)) = (
        is_cli_value_set(merge_matches, "dump_path"),
        config.dump_path,
    ) {
        let resolved = resolve_config_relative_path(config_path, dump_path);
        validate_parent_dir_exists(&resolved)?;
        merge_args.merge_args.dump_path = Some(resolved);
    }

    if let (false, Some(max_dist)) = (is_cli_value_set(merge_matches, "max_dist"), config.max_dist)
    {
        merge_args.merge_args.max_dist = max_dist_in_range(&max_dist.to_string())?;
    }

    if let (false, Some(min_dist)) = (is_cli_value_set(merge_matches, "min_dist"), config.min_dist)
    {
        merge_args.merge_args.min_dist = min_dist_in_range(&min_dist.to_string())?;
    }

    if let (false, Some(min_sequence_similarity)) = (
        is_cli_value_set(merge_matches, "min_sequence_similarity"),
        config.min_sequence_similarity,
    ) {
        merge_args.merge_args.min_sequence_similarity = min_sequence_similarity;
    }

    if let (false, Some(min_size_similarity)) = (
        is_cli_value_set(merge_matches, "min_size_similarity"),
        config.min_size_similarity,
    ) {
        merge_args.merge_args.min_size_similarity = min_size_similarity;
    }

    if let (false, Some(max_dist_linear)) = (
        is_cli_value_set(merge_matches, "max_dist_linear"),
        config.max_dist_linear,
    ) {
        merge_args.merge_args.max_dist_linear = max_dist_linear;
    }

    if let (false, Some(disable_linear_threshold)) = (
        is_cli_value_set(merge_matches, "use_linear_threshold"),
        config.disable_linear_threshold,
    ) {
        merge_args.merge_args.use_linear_threshold = !disable_linear_threshold;
    }

    if let (false, Some(knn_search_k)) = (
        is_cli_value_set(merge_matches, "knn_search_k"),
        config.knn_search_k,
    ) {
        merge_args.merge_args.knn_search_k = knn_search_k_in_range(&knn_search_k.to_string())?;
    }

    if let (false, Some(no_shard)) = (is_cli_value_set(merge_matches, "no_shard"), config.no_shard)
    {
        merge_args.merge_args.no_shard = no_shard;
    }

    if let (false, Some(min_shard_size)) = (
        is_cli_value_set(merge_matches, "min_shard_size"),
        config.min_shard_size,
    ) {
        merge_args.merge_args.min_shard_size = min_shard_size;
    }

    if let (false, Some(min_reciprocal_overlap)) = (
        is_cli_value_set(merge_matches, "min_recip_overlap"),
        config.min_reciprocal_overlap,
    ) {
        merge_args.merge_args.min_recip_overlap = min_reciprocal_overlap;
    }

    if let (false, Some(allow_intrasample)) = (
        is_cli_value_set(merge_matches, "allow_intrasample"),
        config.allow_intrasample,
    ) {
        merge_args.merge_args.allow_intrasample = allow_intrasample;
    }

    if let (false, Some(min_supp)) = (is_cli_value_set(merge_matches, "min_supp"), config.min_supp)
    {
        merge_args.merge_args.min_supp = positive_usize(&min_supp.to_string())?;
    }

    if let (false, Some(keep_monomorphic)) = (
        is_cli_value_set(merge_matches, "keep_monomorphic"),
        config.keep_monomorphic,
    ) {
        merge_args.merge_args.keep_monomorphic = keep_monomorphic;
    }

    if let (false, Some(no_mutual_distance)) = (
        is_cli_value_set(merge_matches, "require_mutual_distance"),
        config.no_mutual_distance,
    ) {
        merge_args.merge_args.require_mutual_distance = !no_mutual_distance;
    }

    if let (false, Some(merge_constraint)) = (
        is_cli_value_set(merge_matches, "merge_constraint"),
        config.merge_constraint,
    ) {
        merge_args.merge_args.merge_constraint = parse_config_merge_constraint(&merge_constraint)?;
    }

    if let (false, Some(filter_tr_contained)) = (
        is_cli_value_set(merge_matches, "filter_tr_contained"),
        config.filter_tr_contained,
    ) {
        merge_args.merge_args.filter_tr_contained = filter_tr_contained;
    }

    if let (false, Some(tr_span_query_slop)) = (
        is_cli_value_set(merge_matches, "tr_span_query_slop"),
        config.tr_span_query_slop,
    ) {
        merge_args.merge_args.tr_span_query_slop = tr_span_query_slop;
    }

    if let (false, Some(tr_min_span_containment)) = (
        is_cli_value_set(merge_matches, "tr_min_span_containment"),
        config.tr_min_span_containment,
    ) {
        merge_args.merge_args.tr_min_span_containment = tr_min_span_containment;
    }

    if let (false, Some(tr_min_span_overlap_bp)) = (
        is_cli_value_set(merge_matches, "tr_min_span_overlap_bp"),
        config.tr_min_span_overlap_bp,
    ) {
        merge_args.merge_args.tr_min_span_overlap_bp = tr_min_span_overlap_bp;
    }

    if let (false, Some(tr_ins_max_dist)) = (
        is_cli_value_set(merge_matches, "tr_ins_max_dist"),
        config.tr_ins_max_dist,
    ) {
        merge_args.merge_args.tr_ins_max_dist = tr_ins_max_dist;
    }

    if let (false, Some(tr_max_dist)) = (
        is_cli_value_set(merge_matches, "tr_max_dist"),
        config.tr_max_dist,
    ) {
        merge_args.merge_args.tr_max_dist = tr_max_dist_in_range(&tr_max_dist.to_string())?;
    }

    if let (false, Some(tr_min_sequence_similarity)) = (
        is_cli_value_set(merge_matches, "tr_min_sequence_similarity"),
        config.tr_min_sequence_similarity,
    ) {
        merge_args.merge_args.tr_min_sequence_similarity = tr_min_sequence_similarity;
    }

    if let (false, Some(tr_min_recip_overlap)) = (
        is_cli_value_set(merge_matches, "tr_min_recip_overlap"),
        config.tr_min_recip_overlap,
    ) {
        merge_args.merge_args.tr_min_recip_overlap = tr_min_recip_overlap;
    }

    Ok(())
}

fn try_parse_cli_with_config_from<I, T>(args: I) -> std::result::Result<Cli, clap::Error>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let args: Vec<OsString> = args.into_iter().map(Into::into).collect();
    let matches = Cli::command().try_get_matches_from(args)?;
    let mut cli = Cli::from_arg_matches(&matches)?;

    let Command::Merge(merge_args) = &mut cli.command;
    if let Some(merge_matches) = matches.subcommand_matches("merge") {
        apply_merge_toml_config(merge_args, merge_matches).map_err(|error| {
            clap::Error::raw(
                clap::error::ErrorKind::ValueValidation,
                format!("{error:#}"),
            )
        })?;
        validate_sort_tuning(merge_args).map_err(|error| {
            clap::Error::raw(
                clap::error::ErrorKind::ValueValidation,
                format!("{error:#}"),
            )
        })?;
        distinct_output_and_dump_paths(merge_args).map_err(|error| {
            clap::Error::raw(
                clap::error::ErrorKind::ValueValidation,
                format!("{error:#}"),
            )
        })?;
    }

    Ok(cli)
}

pub fn parse_cli_with_config() -> Cli {
    match try_parse_cli_with_config_from(std::env::args_os()) {
        Ok(cli) => cli,
        Err(error) => error.exit(),
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Color {
    Always,
    Auto,
    Never,
}

impl Color {
    fn apply(self) {
        match self {
            Color::Always => owo_colors::set_override(true),
            Color::Auto => {}
            Color::Never => owo_colors::set_override(false),
        }
    }
}

pub fn init_verbose(args: &Cli) {
    args.color.apply();

    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::Builder::from_default_env()
        .format(format_log)
        .filter_level(filter_level)
        .init();
}

#[inline(always)]
fn level_style(level: Level) -> (&'static str, Style) {
    match level {
        Level::Error => ("ERROR", Style::new().fg::<Red>().bold()),
        Level::Warn => ("WARN", Style::new().fg::<Yellow>()),
        Level::Info => ("INFO", Style::new().fg::<Green>()),
        Level::Debug => ("DEBUG", Style::new().fg::<Blue>()),
        Level::Trace => ("TRACE", Style::new().fg::<Magenta>()),
    }
}

fn format_log(buf: &mut env_logger::fmt::Formatter, record: &log::Record) -> std::io::Result<()> {
    let (label, style) = level_style(record.level());
    let ts = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
    let painted_label = label.if_supports_color(Stream::Stderr, |t| style.style(t));
    writeln!(buf, "{ts} [{}] - {}", painted_label, record.args())
}

fn threads_in_range(s: &str) -> Result<usize> {
    let thread: usize = s
        .parse::<usize>()
        .map_err(|_| anyhow!("`{}` is not a valid thread number", s))?;
    if thread == 0 {
        return Err(anyhow!("Number of threads must be >= 1"));
    }
    Ok(thread)
}

fn queue_capacity_in_range(s: &str) -> Result<usize> {
    let capacity = s
        .parse::<usize>()
        .map_err(|_| anyhow!("`{}` is not a valid queue capacity", s))?;
    if capacity == 0 {
        return Err(anyhow!("Queue capacity must be >= 1"));
    }
    if capacity > MAX_PIPELINE_QUEUE_CAPACITY {
        return Err(anyhow!(
            "Queue capacity must be <= {}",
            MAX_PIPELINE_QUEUE_CAPACITY
        ));
    }
    Ok(capacity)
}

fn check_file_exists(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if !path.exists() {
        return Err(anyhow!("File does not exist: {}", path.display()));
    }
    Ok(path.to_path_buf())
}

fn check_existing_dir_path(s: &str) -> Result<PathBuf> {
    let path = PathBuf::from(s);
    validate_existing_dir_path(&path)?;
    Ok(path)
}

fn check_prefix_path(s: &str) -> Result<String> {
    let path = Path::new(s);
    validate_parent_dir_exists(path)?;
    Ok(s.to_string())
}

fn check_prefix_pathbuf(s: &str) -> Result<PathBuf> {
    let path = PathBuf::from(s);
    validate_parent_dir_exists(&path)?;
    Ok(path)
}

fn validate_parent_dir_exists(path: &Path) -> Result<()> {
    match path.parent() {
        Some(parent_dir) if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() => {
            Err(anyhow!("Path does not exist: {}", parent_dir.display()))
        }
        _ => Ok(()),
    }
}

fn validate_existing_dir_path(path: &Path) -> Result<()> {
    if !path.exists() {
        return Err(anyhow!("Path does not exist: {}", path.display()));
    }
    if !path.is_dir() {
        return Err(anyhow!("Path is not a directory: {}", path.display()));
    }
    Ok(())
}

fn positive_usize(s: &str) -> Result<usize> {
    let value = s
        .parse::<usize>()
        .map_err(|_| anyhow!("`{s}` is not a valid positive integer"))?;
    if value == 0 {
        return Err(anyhow!("Value must be >= 1"));
    }
    Ok(value)
}

fn knn_search_k_in_range(s: &str) -> Result<usize> {
    let value = positive_usize(s)?;
    if value > MAX_KNN_SEARCH_K {
        return Err(anyhow!("knn_search_k must be <= {MAX_KNN_SEARCH_K}"));
    }
    Ok(value)
}

fn parse_sort_max_mem(value: &str) -> Result<usize> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Err(anyhow!("sort-max-mem cannot be empty"));
    }

    let (numeric, multiplier) = match trimmed.chars().last() {
        Some(last) if last.is_ascii_alphabetic() => {
            let unit = last.to_ascii_lowercase();
            let multiplier = match unit {
                'k' => 1_000_f64,
                'm' => 1_000_000_f64,
                'g' => 1_000_000_000_f64,
                _ => {
                    return Err(anyhow!(
                        "Invalid sort-max-mem suffix `{last}`. Supported suffixes: k, m, g"
                    ));
                }
            };
            (&trimmed[..trimmed.len() - last.len_utf8()], multiplier)
        }
        Some(_) => (trimmed, 1_f64),
        None => return Err(anyhow!("sort-max-mem cannot be empty")),
    };

    let numeric = numeric.trim();
    if numeric.is_empty() {
        return Err(anyhow!("sort-max-mem must include a numeric value"));
    }

    let parsed = numeric
        .parse::<f64>()
        .map_err(|_| anyhow!("`{value}` is not a valid sort-max-mem value"))?;
    if !parsed.is_finite() || parsed <= 0.0 {
        return Err(anyhow!("sort-max-mem must be > 0"));
    }

    let bytes = parsed * multiplier;
    if !bytes.is_finite() || bytes < 1.0 || bytes > usize::MAX as f64 {
        return Err(anyhow!(
            "sort-max-mem `{value}` is outside the supported range"
        ));
    }

    Ok(bytes as usize)
}

fn validate_sort_tuning(args: &MergeArgs) -> Result<()> {
    if args.sort_merge_fan_in > args.sort_max_open_files {
        return Err(anyhow!(
            "sort-merge-fan-in ({}) must be <= sort-max-open-files ({})",
            args.sort_merge_fan_in,
            args.sort_max_open_files
        ));
    }
    Ok(())
}

fn distinct_output_and_dump_paths(args: &MergeArgs) -> Result<()> {
    match (args.output.as_deref(), args.merge_args.dump_path.as_deref()) {
        (Some(output_path), Some(dump_path)) if Path::new(output_path) == dump_path => {
            Err(anyhow!(
                "The dump path and output path must be different: {}",
                dump_path.display()
            ))
        }
        _ => Ok(()),
    }
}

fn min_dist_in_range(s: &str) -> Result<i32> {
    let dist: i32 = s
        .parse::<i32>()
        .map_err(|_| anyhow!("`{}` is not a valid distance", s))?;
    if dist < -1 {
        return Err(anyhow!(
            "min_dist must be >= 0 or exactly -1 (for no minimum)"
        ));
    }
    Ok(dist)
}

fn parse_non_negative_distance(s: &str, arg_name: &str) -> Result<i32> {
    let dist = s
        .parse::<i32>()
        .map_err(|_| anyhow!("`{}` is not a valid distance", s))?;
    if dist < 0 {
        return Err(anyhow!("{arg_name} must be >= 0"));
    }
    Ok(dist)
}

fn max_dist_in_range(s: &str) -> Result<i32> {
    parse_non_negative_distance(s, "max_dist")
}

fn tr_max_dist_in_range(s: &str) -> Result<i32> {
    parse_non_negative_distance(s, "tr_max_dist")
}

fn merge_validate_output_type(s: &str) -> Result<OutputType> {
    let valid_prefixes = ["u", "b", "v", "z"];
    if valid_prefixes.contains(&s) {
        return match s {
            "u" => Ok(OutputType::Bcf {
                is_uncompressed: true,
                level: None,
            }),
            "v" => Ok(OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }),
            "b" => Ok(OutputType::Bcf {
                is_uncompressed: false,
                level: None,
            }),
            "z" => Ok(OutputType::Vcf {
                is_uncompressed: false,
                level: None,
            }),
            _ => unreachable!(),
        };
    }

    // NOTE: Can't actually set compression level in rust/htslib at the moment
    // if s.len() == 2 {
    //     let (prefix, suffix) = s.split_at(1);
    //     if (prefix == "b" || prefix == "z") && suffix.chars().all(|c| c.is_digit(10)) {
    //         return match prefix {
    //             "b" => Ok(OutputType::Bcf {
    //                 is_uncompressed: false,
    //                 level: Some(suffix.parse().unwrap()),
    //             }),
    //             "z" => Ok(OutputType::Vcf {
    //                 is_uncompressed: false,
    //                 level: Some(suffix.parse().unwrap()),
    //             }),
    //             _ => unreachable!(),
    //         };
    //     } else if (prefix == "u" || prefix == "v") && suffix.chars().all(|c| c.is_digit(10)) {
    //         return Err(format!(
    //             "Error: compression level ({}) cannot be set on uncompressed streams ({})",
    //             suffix, prefix
    //         ));
    //     }
    // }

    Err(anyhow!(
        "Invalid output type: {}. Must be one of u, b, v, z.",
        s
    ))
}

impl MergeArgs {
    pub fn selected_svtypes(&self) -> Vec<SvType> {
        let requested: HashSet<MergeSvType> = self.svtypes.iter().copied().collect();
        let include_default_all = requested.contains(&MergeSvType::All);

        MergeSvType::ORDERED_SELECTABLE_NON_ALL
            .into_iter()
            .filter(|svtype| {
                requested.contains(svtype)
                    || (include_default_all && MergeSvType::ORDERED_DEFAULT_ALL.contains(svtype))
            })
            .filter_map(MergeSvType::as_svtype)
            .collect()
    }

    pub fn process_vcf_paths(&self) -> Result<Vec<PathBuf>> {
        match (&self.vcfs, &self.vcf_list) {
            (Some(vcfs), None) => Ok(vcfs.clone()),
            (None, Some(list_path)) => Self::read_vcf_paths_from_file(list_path),
            _ => unreachable!("Either --vcf or --vcf-list is provided, never both"),
        }
    }

    fn read_vcf_paths_from_file(path: &Path) -> Result<Vec<PathBuf>> {
        let file = File::open(path)
            .map_err(|e| anyhow!("Failed to open VCF list file {}: {}", path.display(), e))?;
        let reader = BufReader::new(file);

        let mut paths = Vec::new();
        for (line_num, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| anyhow!("Error reading line {}: {}", line_num + 1, e))?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let path = PathBuf::from(trimmed);
            if !path.exists() {
                Err(anyhow!("VCF file does not exist: {}", path.display()))?;
            }
            paths.push(path);
        }

        if paths.is_empty() {
            Err(anyhow!("No VCF paths found in the input file".to_string()))?;
        }

        Ok(paths)
    }
}
