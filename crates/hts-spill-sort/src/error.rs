use std::num::TryFromIntError;
use thiserror::Error;

pub type Result<T> = std::result::Result<T, SpillSortError>;

#[derive(Debug, Error)]
pub enum SpillSortError {
    #[error("{0}")]
    Message(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Htslib(#[from] rust_htslib::errors::Error),
    #[error(transparent)]
    TryFromInt(#[from] TryFromIntError),
}

impl SpillSortError {
    pub fn message(message: impl Into<String>) -> Self {
        Self::Message(message.into())
    }
}
