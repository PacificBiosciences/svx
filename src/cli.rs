use crate::{
    constants::*,
    io::{
        positions_reader::{parse_position_tuples, PositionTuple},
        vcf_writer::OutputType,
    },
    utils::util::Result,
};
use anyhow::anyhow;
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, Parser, Subcommand};
use env_logger::fmt::Color;
use log::{Level, LevelFilter};
use once_cell::sync::Lazy;
use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

/// Full version string including the crate version and git description.
///
/// This version string is used in the command-line interface to provide detailed version information.
/// It includes the crate version from Cargo.toml and additional build information such as the git commit hash.
/// # Examples
/// * `0.1.0-1ba958a-dirty` - while on a dirty branch
/// * `0.1.0-1ba958a` - with a fresh commit
pub static FULL_VERSION: Lazy<String> = Lazy::new(|| {
    let git_describe = env!("VERGEN_GIT_DESCRIBE");
    if git_describe.is_empty() {
        env!("CARGO_PKG_VERSION").to_string()
    } else {
        format!("{}-{}", env!("CARGO_PKG_VERSION"), git_describe)
    }
});

#[derive(Parser, Debug)]
#[command(name="svx",
          author="Tom Mokveld <tmokveld@pacificbiosciences.com>", 
          version=&**FULL_VERSION,
          about="Structural variant merger",
          long_about = None,
          after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
          This program comes with ABSOLUTELY NO WARRANTY; it is intended for
          Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
    #[arg(
        short = 'v',
        long = "verbose",
        action = ArgAction::Count,
        global = true
    )]
    pub verbosity: u8,
}

#[derive(Subcommand, Debug)]
pub enum Command {
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
        value_parser = check_file_exists
    )]
    pub vcfs: Option<Vec<PathBuf>>,

    /// File containing paths of VCF files to merge (one per line)
    #[arg(
        long = "vcf-list",
        value_name = "VCF_LIST",
        value_parser = check_file_exists
    )]
    pub vcf_list: Option<PathBuf>,

    /// Write output to a file [default: standard output]
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        value_parser = check_prefix_path
    )]
    pub output: Option<String>,

    /// BED file with tandem repeat definitions
    #[arg(
        long = "trs",
        value_name = "TRS",
        value_parser = check_file_exists,
    )]
    pub tr_bed_path: Option<PathBuf>,

    /// Number of threads to use
    #[arg(
        short = '@',
        value_name = "THREADS",
        default_value = "1",
        value_parser = threads_in_range
    )]
    pub num_threads: usize,

    /// Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "OUTPUT_TYPE",
        value_parser = merge_validate_output_type,
        help_heading = "Advanced"
    )]
    pub output_type: Option<OutputType>,

    /// Print only the merged header and exit
    #[arg(long = "print-header", help_heading = "Advanced")]
    pub print_header: bool,

    /// Run even if there is only one file on input
    #[arg(long = "force-single", help_heading = "Advanced")]
    pub force_single: bool,

    /// Resolve duplicate sample names
    #[arg(long = "force-samples", help_heading = "Advanced", hide = true)]
    pub force_samples: bool,

    /// Do not append version and command line to the header
    #[arg(long = "no-version", help_heading = "Advanced")]
    pub no_version: bool,

    /// Process only the specified contigs (comma-separated list), e.g., (chr1,chr2,chrX)
    #[arg(
        long = "contig",
        value_name = "CONTIG",
        value_delimiter = ',',
        help_heading = "Advanced"
    )]
    pub contigs: Option<Vec<String>>,

    /// Specific positions to merge in format (contig:start[-end]) e.g., (chr1:12345),(chr2:67890-67900)
    #[arg(
        long = "target-positions",
        value_name = "TARGET-POSITIONS",
        value_parser = parse_position_tuples,
        help_heading = "Advanced"
    )]
    pub target_positions: Option<std::vec::Vec<PositionTuple>>,

    #[command(flatten)]
    pub merge_args: MergeArgsInner,
}

#[derive(Parser, Debug, Clone)]
pub struct MergeArgsInner {
    /// Dump coordinates of variants
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_DUMP
    )]
    pub dump: bool,

    /// Maximum distance for merging variants
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_MAX_DIST
    )]
    pub max_dist: i32,

    /// The minimum distance threshold a variant can have when using max_dist_linear (-1 for no minimum)
    #[arg(
        help_heading = "Advanced",
        long,
        default_value_t = DEFAULT_MIN_DIST,
        value_parser = min_dist_in_range
    )]
    pub min_dist: i32,

    /// The minimum sequence identity for two insertions to be merged
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_MIN_SEQUENCE_SIMILARITY
    )]
    pub min_sequence_similarity: f32,

    /// The proportion of the length of each variant to set distance threshold to
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_MAX_DIST_LINEAR
    )]
    pub max_dist_linear: f32,

    /// Disable distance threshold depending on variant length and use max_dist instead
    #[arg(
        help_heading("Advanced"),
        long = "disable-linear-threshold", 
        action = ArgAction::SetFalse,
        default_value_t = DEFAULT_USE_LINEAR_THRESHOLD
    )]
    pub use_linear_threshold: bool,

    /// Norm to use for KD-tree distance calculations
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_KD_TREE_NORM
    )]
    pub kd_tree_norm: i32,

    /// Number of nearest neighbors to search in KD-tree
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_KNN_SEARCH_K
    )]
    pub knn_search_k: usize,

    /// The minimum reciprocal overlap for DEL/INV/DUP SVs
    #[arg(
        help_heading("Advanced"),
        long = "min-reciprocal-overlap",
        default_value_t = DEFAULT_MIN_RECIP_OVERLAP
    )]
    pub min_recip_overlap: f32,

    /// Allow merging variants from the same sample
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_ALLOW_INTRASAMPLE
    )]
    pub allow_intrasample: bool,

    /// Require mutual distance for merging
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_REQUIRE_MUTUAL_DISTANCE
    )]
    pub require_mutual_distance: bool,

    /// Filter out all SVs that are TR contained
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_FILTER_TR_CONTAINED
    )]
    pub filter_tr_contained: bool,

    /// Maximum ratio of SV span to TR span for an SV to be considered TR-contained (SV_span / TR_span <= 1.0 + value)
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_TR_SPAN_RATIO_THRESHOLD
    )]
    pub tr_span_ratio_threshold: f32,

    /// Maximum distance for merging TR-contained variants (if they share the same TR ID).
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_TR_MAX_DIST
    )]
    pub tr_max_dist: i32,

    /// Minimum sequence similarity for merging TR-contained variants (if they share the same TR ID).
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_TR_MIN_SEQUENCE_SIMILARITY
    )]
    pub tr_min_sequence_similarity: f32,

    /// Minimum reciprocal overlap for merging TR-contained DEL/INV/DUP SVs (if they share the same TR ID).
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_TR_MIN_RECIP_OVERLAP
    )]
    pub tr_min_recip_overlap: f32,

    /// Number of iterations for cluster refinement
    #[arg(
        help_heading("Advanced"),
        long,
        default_value_t = DEFAULT_REFINEMENT_ITERATIONS
    )]
    pub refinement_iterations: usize,
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
            dump: DEFAULT_DUMP,
            max_dist: DEFAULT_MAX_DIST,
            min_dist: DEFAULT_MIN_DIST,
            max_dist_linear: DEFAULT_MAX_DIST_LINEAR,
            use_linear_threshold: DEFAULT_USE_LINEAR_THRESHOLD,
            min_sequence_similarity: DEFAULT_MIN_SEQUENCE_SIMILARITY,
            kd_tree_norm: DEFAULT_KD_TREE_NORM,
            knn_search_k: DEFAULT_KNN_SEARCH_K,
            allow_intrasample: DEFAULT_ALLOW_INTRASAMPLE,
            require_mutual_distance: DEFAULT_REQUIRE_MUTUAL_DISTANCE,
            min_recip_overlap: DEFAULT_MIN_RECIP_OVERLAP,
            filter_tr_contained: DEFAULT_FILTER_TR_CONTAINED,
            tr_span_ratio_threshold: DEFAULT_TR_SPAN_RATIO_THRESHOLD,
            tr_max_dist: DEFAULT_TR_MAX_DIST,
            tr_min_sequence_similarity: DEFAULT_TR_MIN_SEQUENCE_SIMILARITY,
            tr_min_recip_overlap: DEFAULT_TR_MIN_RECIP_OVERLAP,
            refinement_iterations: DEFAULT_REFINEMENT_ITERATIONS,
        }
    }
}

/// Initializes the verbosity level for logging based on the command-line arguments.
///
/// Sets up the logger with a specific verbosity level that is determined
/// by the number of occurrences of the `-v` or `--verbose` flag in the command-line arguments.
///
/// # Arguments
///
/// * `args` - A reference to the parsed command-line arguments.
pub fn init_verbose(args: &Cli) {
    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::Builder::from_default_env()
        .format(|buf, record| {
            let level = record.level();
            let mut style = buf.style();
            match record.level() {
                Level::Error => style.set_color(Color::Red),
                Level::Warn => style.set_color(Color::Yellow),
                Level::Info => style.set_color(Color::Green),
                Level::Debug => style.set_color(Color::Blue),
                Level::Trace => style.set_color(Color::Cyan),
            };

            writeln!(
                buf,
                "{} [{}] {} - {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                style.value(level),
                record.module_path().unwrap_or("unknown_module"),
                record.args()
            )
        })
        .filter_level(filter_level)
        .init();
}

/// Validates that the provided string represents a valid number of threads.
///
/// Checks if the string argument can be parsed into a non-zero positive integer
/// that represents the number of threads to use. It ensures that the number of threads is within
/// a valid range.
///
/// # Arguments
///
/// * `s` - A string slice representing the number of threads.
///
/// # Returns
///
/// Returns a `Result<usize>` which is Ok if the number is valid, or an Err with a descriptive message if not.
fn threads_in_range(s: &str) -> Result<usize> {
    let thread: usize = s
        .parse::<usize>()
        .map_err(|_| anyhow!("`{}` is not a valid thread number", s))?;
    if thread == 0 {
        return Err(anyhow!("Number of threads must be >= 1"));
    }
    Ok(thread)
}

/// Checks if the provided file path exists.
///
/// Validates that the file path provided as an argument exists in the file system.
/// It is used to ensure that the file paths provided for input files are valid before attempting to process them.
///
/// # Arguments
///
/// * `s` - A string slice representing the file path to check.
///
/// # Returns
///
/// Returns a `Result<PathBuf>` which is Ok if the file exists, or an Err with a descriptive message if not.
fn check_file_exists(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if !path.exists() {
        return Err(anyhow!("File does not exist: {}", path.display()));
    }
    Ok(path.to_path_buf())
}

fn check_prefix_path(s: &str) -> Result<String> {
    let path = Path::new(s);
    if let Some(parent_dir) = path.parent() {
        if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() {
            return Err(anyhow!("Path does not exist: {}", parent_dir.display()));
        }
    }
    Ok(s.to_string())
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
            // Skip empty or comment lines
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
