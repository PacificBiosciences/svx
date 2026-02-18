pub mod cli;
pub mod commands;
pub mod error;

pub mod core {
    pub mod aligner;
    pub mod svtype;
    pub mod variant;
    pub mod variant_block;
    pub mod containers {
        pub mod forest;
        pub mod interval_tree;
        pub mod kd_tree;
    }
}

pub mod io {
    pub mod bed_reader;
    pub mod merge_reader;
    pub mod merge_writer;
    pub mod output_sort;
    pub mod positions_reader;
    pub mod readers;
    pub mod tpool;
    pub mod vcf_reader;
    pub mod vcf_writer;
}

pub mod utils {
    pub mod util;
    pub mod util_intern;
    pub mod wfa2;
}

pub mod constants;

pub use constants::*;
