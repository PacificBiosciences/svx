use crate::{error::SvxError, utils::util::Result};
use rust_htslib::bcf;

#[derive(Debug, Clone, PartialEq, Copy)]
pub enum SvType {
    INSERTION,
    DELETION,
    INVERSION,
    DUPLICATION,
    CNV,
    BND,
}

impl SvType {
    pub const BUCKET_COUNT: usize = 6;
    pub const NON_BND_BUCKET_TYPES: [Self; 5] = [
        SvType::INSERTION,
        SvType::DELETION,
        SvType::INVERSION,
        SvType::DUPLICATION,
        SvType::CNV,
    ];

    pub fn from_u8(bytes: &[u8]) -> Result<Self> {
        match bytes {
            b"INS" => Ok(SvType::INSERTION),
            b"DEL" => Ok(SvType::DELETION),
            b"INV" => Ok(SvType::INVERSION),
            b"DUP" => Ok(SvType::DUPLICATION),
            b"CNV" => Ok(SvType::CNV),
            b"BND" => Ok(SvType::BND),
            _ => Err(SvxError::InvalidSvtype {
                value: String::from_utf8_lossy(bytes).into_owned(),
            }),
        }
    }

    pub fn classify_vcf_record(record: &bcf::Record) -> Result<Self> {
        let svtype = Self::from_vcf_record(record)?;
        if svtype == SvType::CNV || Self::id_is_cnv(record.id().as_slice()) {
            return Ok(SvType::CNV);
        }
        Ok(svtype)
    }

    pub fn id_is_cnv(id: &[u8]) -> bool {
        id.split(|byte| *byte == b':').nth(1) == Some(b"CNV")
    }

    pub fn from_vcf_record(record: &bcf::Record) -> Result<Self> {
        match record.info(b"SVTYPE").string() {
            Ok(Some(svtypes)) => {
                let svtype = svtypes.first().copied().ok_or(SvxError::MissingSvtype)?;
                Self::from_u8(svtype)
            }
            Ok(None) => Err(SvxError::MissingSvtype),
            Err(error) => Err(SvxError::SvtypeInfoRead {
                message: error.to_string(),
            }),
        }
    }

    pub const fn bucket_index(self) -> usize {
        match self {
            SvType::INSERTION => 0,
            SvType::DELETION => 1,
            SvType::INVERSION => 2,
            SvType::DUPLICATION => 3,
            SvType::CNV => 4,
            SvType::BND => 5,
        }
    }
}

impl std::str::FromStr for SvType {
    type Err = SvxError;
    fn from_str(s: &str) -> Result<Self> {
        Self::from_u8(s.as_bytes())
    }
}

impl std::fmt::Display for SvType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SvType::INSERTION => write!(f, "INS"),
            SvType::DELETION => write!(f, "DEL"),
            SvType::INVERSION => write!(f, "INV"),
            SvType::DUPLICATION => write!(f, "DUP"),
            SvType::CNV => write!(f, "CNV"),
            SvType::BND => write!(f, "BND"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::SvxError;
    use rust_htslib::bcf::{self, Read, record::GenotypeAllele};
    use std::time::SystemTime;

    fn make_temp_vcf_with_optional_svtype(svtype: Option<&[u8]>, id: &[u8]) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!("svx_test_svtype_{nanos}.vcf"));

        let mut header = bcf::Header::new();
        header.push_record(br#"##fileformat=VCFv4.2"#);
        header.push_record(br#"##contig=<ID=chr1,length=1000000>"#);
        header.push_record(br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
        header.push_sample(b"sample1");

        let mut writer = bcf::Writer::from_path(&path, &header, false, bcf::Format::Vcf).unwrap();
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(100);
        record.set_id(id).unwrap();
        record.set_alleles(&[b"N", b"<INS>"]).unwrap();
        if let Some(svtype) = svtype {
            record.push_info_string(b"SVTYPE", &[svtype]).unwrap();
        }
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        writer.write(&record).unwrap();
        drop(writer);
        bcf::index::build(&path, None, 1, bcf::index::Type::Tbx).unwrap();

        path
    }

    #[test]
    fn from_vcf_record_parses_svtype_info() {
        let path = make_temp_vcf_with_optional_svtype(Some(b"INS"), b"v1");
        let mut reader = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();
        let rid = reader.header().name2rid(b"chr1").unwrap();
        reader.fetch(rid, 0, None).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let svtype = SvType::from_vcf_record(&record).unwrap();
        assert_eq!(svtype, SvType::INSERTION);
    }

    #[test]
    fn from_vcf_record_errors_when_svtype_missing() {
        let path = make_temp_vcf_with_optional_svtype(None, b"v1");
        let mut reader = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();
        let rid = reader.header().name2rid(b"chr1").unwrap();
        reader.fetch(rid, 0, None).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let err = SvType::from_vcf_record(&record).unwrap_err();
        assert!(matches!(err, SvxError::MissingSvtype));
    }

    #[test]
    fn classify_vcf_record_maps_cnv_id_to_cnv_bucket() {
        let path = make_temp_vcf_with_optional_svtype(Some(b"DEL"), b"sawfish:CNV:0:0:200");
        let mut reader = rust_htslib::bcf::IndexedReader::from_path(&path).unwrap();
        let rid = reader.header().name2rid(b"chr1").unwrap();
        reader.fetch(rid, 0, None).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let classified = SvType::classify_vcf_record(&record).unwrap();
        assert_eq!(classified, SvType::CNV);
    }
}
