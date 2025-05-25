use once_cell::sync::Lazy;
use std::sync::OnceLock;

pub const DEFAULT_MAX_DIST_LINEAR: f32 = 0.5;
pub const DEFAULT_USE_LINEAR_THRESHOLD: bool = true;
pub const DEFAULT_MIN_DIST: i32 = 100; // -1 means no minimum
pub const DEFAULT_MAX_DIST: i32 = 1000;
pub const DEFAULT_MIN_SEQUENCE_SIMILARITY: f32 = 0.8;
pub const DEFAULT_KD_TREE_NORM: i32 = 2;
pub const DEFAULT_KNN_SEARCH_K: usize = 4;
pub const DEFAULT_ALLOW_INTRASAMPLE: bool = false;
pub const DEFAULT_REQUIRE_MUTUAL_DISTANCE: bool = true;
pub const DEFAULT_DUMP: bool = false;
pub const DEFAULT_MIN_RECIP_OVERLAP: f32 = 0.0;
pub const DEFAULT_FILTER_TR_CONTAINED: bool = false;
pub const DEFAULT_TR_SPAN_RATIO_THRESHOLD: f32 = 30.0;
pub const DEFAULT_TR_MAX_DIST: i32 = 4000;
pub const DEFAULT_TR_MIN_SEQUENCE_SIMILARITY: f32 = 0.6;
pub const DEFAULT_TR_MIN_RECIP_OVERLAP: f32 = 0.0;
pub const DEFAULT_REFINEMENT_ITERATIONS: usize = 6;

pub const DISTANCE_OFFSET: f32 = 1e-9f32;

#[derive(Debug, Clone)]
pub struct SvxConfig {
    pub kd_tree_norm: i32,
    pub dump: bool,
}

static CONFIG_INITIALIZER: OnceLock<SvxConfig> = OnceLock::new();

pub static CONFIG: Lazy<SvxConfig> = Lazy::new(|| {
    CONFIG_INITIALIZER
        .get()
        .cloned()
        .expect("CONFIG_INITIALIZER must be initialized before accessing it")
});

pub fn init_config(config: SvxConfig) {
    CONFIG_INITIALIZER
        .set(config)
        .expect("Config has already been initialized");
}
