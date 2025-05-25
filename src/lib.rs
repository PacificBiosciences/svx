pub mod cli;
pub mod commands;

pub mod core {
    pub mod aligner;
    pub mod svtype;
    pub mod variant_block;
    pub mod variant_internal;
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
    pub mod positions_reader;
    pub mod readers;
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
