use std::{
    num::{ParseFloatError, ParseIntError, TryFromIntError},
    path::PathBuf,
    str::Utf8Error,
};
use thiserror::Error;

pub type SvxResult<T> = std::result::Result<T, SvxError>;

#[derive(Debug, Error)]
pub enum SvxError {
    #[error("{0}")]
    Message(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Htslib(#[from] rust_htslib::errors::Error),
    #[error(transparent)]
    Utf8(#[from] Utf8Error),
    #[error(transparent)]
    ParseInt(#[from] ParseIntError),
    #[error(transparent)]
    ParseFloat(#[from] ParseFloatError),
    #[error(transparent)]
    TryFromInt(#[from] TryFromIntError),
    #[error("SVTYPE missing in VCF record")]
    MissingSvtype,
    #[error("Invalid SVTYPE: {value}")]
    InvalidSvtype { value: String },
    #[error("Error reading SVTYPE INFO from VCF record: {message}")]
    SvtypeInfoRead { message: String },
    #[error(
        "Reference index file not found: {}. Create it using 'samtools faidx {}'",
        fai_path.display(),
        reference_path.display()
    )]
    MissingReferenceIndex {
        fai_path: PathBuf,
        reference_path: PathBuf,
    },
    #[error("Invalid gzip header: {}", path.display())]
    InvalidGzipHeader { path: PathBuf },
}

impl SvxError {
    pub fn message(message: impl Into<String>) -> Self {
        Self::Message(message.into())
    }
}

#[macro_export]
macro_rules! svx_error {
    ($($arg:tt)*) => {
        $crate::error::SvxError::message(format!($($arg)*))
    };
}
