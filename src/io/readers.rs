use crate::{error::SvxError, utils::util::Result};
use flate2::read::MultiGzDecoder;
use rust_htslib::faidx;
use std::{
    fs::File,
    io::{BufReader, Read as ioRead},
    path::{Path, PathBuf},
};

pub fn open_catalog_reader(path: &Path) -> Result<BufReader<Box<dyn ioRead>>> {
    fn is_gzipped(path: &Path) -> bool {
        let path_str = path.to_string_lossy().to_lowercase();
        path_str.ends_with(".gz") || path_str.ends_with(".gzip")
    }
    let file = File::open(path)
        .map_err(|error| crate::svx_error!("Failed to open file {}: {error}", path.display()))?;
    if is_gzipped(path) {
        let gz_decoder = MultiGzDecoder::new(file);
        if gz_decoder.header().is_some() {
            Ok(BufReader::new(Box::new(gz_decoder)))
        } else {
            Err(SvxError::InvalidGzipHeader {
                path: path.to_path_buf(),
            })
        }
    } else {
        Ok(BufReader::new(Box::new(file)))
    }
}

pub fn open_genome_reader(path: &Path) -> Result<faidx::Reader> {
    let fai_path = {
        let mut fai_path = path.as_os_str().to_os_string();
        fai_path.push(".fai");
        PathBuf::from(fai_path)
    };
    if !fai_path.exists() {
        return Err(SvxError::MissingReferenceIndex {
            fai_path,
            reference_path: path.to_path_buf(),
        });
    }
    faidx::Reader::from_path(path).map_err(|e| e.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::SvxError;
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;
    use std::path::PathBuf;
    use tempfile::tempdir;

    #[test]
    fn open_genome_reader_returns_error_for_path_without_extension() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let genome_path = temp_dir.path().join("reference");
        let result = open_genome_reader(&genome_path);
        let err = result.unwrap_err();
        assert!(matches!(err, SvxError::MissingReferenceIndex { .. }));
    }

    #[cfg(unix)]
    #[test]
    fn open_genome_reader_returns_error_for_non_utf8_extension() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let genome_name = OsString::from_vec(vec![b'r', b'e', b'f', b'.', 0xff]);
        let genome_path = temp_dir.path().join(PathBuf::from(genome_name));
        let result = open_genome_reader(&genome_path);
        let err = result.unwrap_err();
        assert!(matches!(err, SvxError::MissingReferenceIndex { .. }));
    }
}
