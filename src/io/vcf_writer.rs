use crate::utils::util::Result;
use anyhow::anyhow;
use rust_htslib::bcf;

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

pub struct VcfWriter {
    pub writer: bcf::Writer,
    pub dummy_record: bcf::Record,
}

impl VcfWriter {
    pub fn new(
        header: &bcf::Header,
        output_type: &Option<OutputType>,
        output: Option<&String>,
    ) -> Result<Self> {
        let output_type = match (output_type, output) {
            (Some(output_type), _) => output_type.clone(),
            (None, Some(path)) => Self::infer_output_type_from_extension(path)?,
            (None, None) => OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            },
        };
        log::trace!("{:?}", &output_type);

        let writer = match output {
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
        .map_err(|e| anyhow!("Failed to create writer: {}", e))?;
        let dummy_record = writer.empty_record();
        Ok(VcfWriter {
            writer,
            dummy_record,
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
