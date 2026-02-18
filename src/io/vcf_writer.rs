use crate::{io::tpool::HtsThreadPool, utils::util::Result};
use rust_htslib::bcf;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub enum OutputType {
    Vcf {
        is_uncompressed: bool,
        level: Option<u8>,
    },
    Bcf {
        is_uncompressed: bool,
        level: Option<u8>,
    },
}

const CSI_MIN_SHIFT_DEFAULT: u32 = 14;

type OutputIndexSpec = (&'static str, bcf::index::Type);

impl OutputType {
    fn is_uncompressed(&self) -> bool {
        match self {
            Self::Vcf {
                is_uncompressed, ..
            }
            | Self::Bcf {
                is_uncompressed, ..
            } => *is_uncompressed,
        }
    }

    fn output_index_spec(&self) -> Option<OutputIndexSpec> {
        if self.is_uncompressed() {
            return None;
        }

        match self {
            Self::Vcf { .. } => Some(("tabix", bcf::index::Type::Tbx)),
            Self::Bcf { .. } => Some(("csi", bcf::index::Type::Csi(CSI_MIN_SHIFT_DEFAULT))),
        }
    }
}

pub struct VcfWriter {
    pub writer: bcf::Writer,
    pub dummy_record: bcf::Record,
    _io_tpool: Option<Arc<HtsThreadPool>>,
}

impl VcfWriter {
    fn resolve_output_type(
        output_type: &Option<OutputType>,
        output: Option<&str>,
    ) -> Result<OutputType> {
        let resolved = match (output_type, output) {
            (Some(output_type), _) => output_type.clone(),
            (None, Some(path)) => Self::infer_output_type_from_extension(path)?,
            (None, None) => OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            },
        };
        Ok(resolved)
    }

    fn output_index_spec_for_target(
        output_type: &Option<OutputType>,
        output: Option<&str>,
    ) -> Result<Option<OutputIndexSpec>> {
        if output.is_none() {
            return Ok(None);
        }

        let resolved_output_type = Self::resolve_output_type(output_type, output)?;
        Ok(resolved_output_type.output_index_spec())
    }

    pub fn build_sorted_output_index(
        output_type: &Option<OutputType>,
        output: Option<&str>,
    ) -> Result<()> {
        let Some(path) = output else {
            return Ok(());
        };
        let Some((index_name, index_type)) =
            Self::output_index_spec_for_target(output_type, output)?
        else {
            return Ok(());
        };

        log::debug!(
            "Writer: Building {} index for sorted output {}",
            index_name,
            path
        );
        bcf::index::build(path, None, 1, index_type).map_err(|error| {
            crate::svx_error!(
                "Failed to build {} index for sorted output {}: {}",
                index_name,
                path,
                error
            )
        })?;
        Ok(())
    }

    pub fn new(
        header: &bcf::Header,
        output_type: &Option<OutputType>,
        output: Option<&String>,
        io_tpool: Option<Arc<HtsThreadPool>>,
    ) -> Result<Self> {
        let output_type = Self::resolve_output_type(output_type, output.map(String::as_str))?;
        log::trace!("{:?}", &output_type);

        let mut writer = match output {
            Some(path) => {
                let (is_uncompressed, format) = match output_type {
                    OutputType::Vcf {
                        is_uncompressed, ..
                    } => (is_uncompressed, bcf::Format::Vcf),
                    OutputType::Bcf {
                        is_uncompressed, ..
                    } => (is_uncompressed, bcf::Format::Bcf),
                };
                bcf::Writer::from_path(path, header, is_uncompressed, format)
            }
            None => {
                let (is_uncompressed, format) = match output_type {
                    OutputType::Vcf {
                        is_uncompressed, ..
                    } => (is_uncompressed, bcf::Format::Vcf),
                    OutputType::Bcf {
                        is_uncompressed, ..
                    } => (is_uncompressed, bcf::Format::Bcf),
                };
                bcf::Writer::from_stdout(header, is_uncompressed, format)
            }
        }
        .map_err(|e| crate::svx_error!("Failed to create writer: {}", e))?;

        if let Some(ref pool) = io_tpool {
            unsafe {
                pool.attach_to_writer(&mut writer).map_err(|e| {
                    crate::svx_error!(
                        "Failed to attach shared HTS thread pool to output writer: {}",
                        e
                    )
                })?;
                log::trace!("Attached shared HTS thread pool to output writer");
            }
        }

        let dummy_record = writer.empty_record();
        Ok(VcfWriter {
            writer,
            dummy_record,
            _io_tpool: io_tpool,
        })
    }

    fn infer_output_type_from_extension(path: &str) -> Result<OutputType> {
        let path_lower = path.to_lowercase();
        match path_lower.as_str() {
            s if s.ends_with(".bcf.gz") => Ok(OutputType::Bcf {
                is_uncompressed: false,
                level: None,
            }),
            s if s.ends_with(".vcf.gz") || s.ends_with(".vcf.bgz") => Ok(OutputType::Vcf {
                is_uncompressed: false,
                level: None,
            }),
            s if s.ends_with(".bcf") => Ok(OutputType::Bcf {
                is_uncompressed: true,
                level: None,
            }),
            s if s.ends_with(".vcf") => Ok(OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }),
            _ => Ok(OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::tpool::HtsThreadPool;
    use std::{
        env::temp_dir,
        fs,
        sync::Arc,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn create_temp_path(ext: &str) -> String {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let mut path = temp_dir();
        path.push(format!("svx_test_vcf_writer_{nanos}.{ext}"));
        path.to_string_lossy().into_owned()
    }

    fn make_header() -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(br#"##fileformat=VCFv4.3"#);
        header.push_record(br#"##contig=<ID=chr1,length=1000>"#);
        header.push_sample(b"S1");
        header
    }

    #[test]
    fn vcf_writer_keeps_shared_io_pool_when_provided() {
        let header = make_header();
        let out_path = create_temp_path("vcf.gz");
        let io_tpool = Arc::new(HtsThreadPool::new(2).unwrap());

        let writer =
            VcfWriter::new(&header, &None, Some(&out_path), Some(io_tpool.clone())).unwrap();

        let writer_pool = writer._io_tpool.as_ref().unwrap();
        assert!(Arc::ptr_eq(writer_pool, &io_tpool));

        fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn vcf_writer_has_no_io_pool_when_not_provided() {
        let header = make_header();
        let out_path = create_temp_path("vcf.gz");

        let writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        assert!(writer._io_tpool.is_none());

        fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn output_index_spec_is_generated_for_compressed_vcf_file_output() {
        let out_path = create_temp_path("vcf.gz");
        let index_spec = VcfWriter::output_index_spec_for_target(&None, Some(out_path.as_str()))
            .expect("compressed output should create index spec");

        assert!(matches!(index_spec, Some(("tabix", bcf::index::Type::Tbx))));
    }

    #[test]
    fn output_index_spec_is_none_for_uncompressed_or_stdout_output() {
        let uncompressed_path = create_temp_path("vcf");
        let uncompressed_spec =
            VcfWriter::output_index_spec_for_target(&None, Some(uncompressed_path.as_str()))
                .expect("uncompressed output should not error");
        assert!(uncompressed_spec.is_none());

        let stdout_spec = VcfWriter::output_index_spec_for_target(&None, None)
            .expect("stdout output should not error");
        assert!(stdout_spec.is_none());
    }
}
