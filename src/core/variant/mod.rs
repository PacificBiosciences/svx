mod bnd;
mod similarity;
mod tr;
mod vcf;

#[cfg(test)]
pub mod test_utils;

use crate::core::svtype::SvType;
use crate::io::bed_reader::TrId;
use rust_htslib::bcf::{header::TagType, record::GenotypeAllele};
use std::{collections::HashMap, fmt};

pub trait VariantSource {
    fn as_variant(&self) -> &VariantInternal;
}

impl VariantSource for VariantInternal {
    fn as_variant(&self) -> &VariantInternal {
        self
    }
}

impl VariantSource for &VariantInternal {
    fn as_variant(&self) -> &VariantInternal {
        self
    }
}

#[derive(Debug, Clone)]
pub struct HeaderFormatCache {
    pub tag_by_id: HashMap<i32, Vec<u8>>,
    pub type_by_tag: HashMap<Vec<u8>, TagType>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BndData {
    pub contig: String,
    pub pos0: i64,
    pub mate_id: String,
    pub mate_contig: String,
    pub mate_pos: i64,
    pub strands: [u8; 2],
    pub homlen: Option<i64>,
    pub inv_event_id: Option<String>,
    pub has_inv_breakpoint_filter: bool,
}

#[derive(Debug, Clone)]
pub struct BndEventData {
    pub a_contig: String,
    pub a_pos0: i64,
    pub a_strands: [u8; 2],
    pub a_ref_base: u8,
    pub a_vcf: VcfWriteData,
    pub b_contig: String,
    pub b_pos0: i64,
    pub b_strands: [u8; 2],
    pub b_ref_base: u8,
    pub b_vcf: VcfWriteData,
    pub inv_event_id: Option<String>,
}

#[derive(Debug, Clone)]
pub struct VcfWriteData {
    pub rid: Option<u32>,
    pub pos: i64,
    pub alleles: Vec<Vec<u8>>,
    pub gt: Vec<GenotypeAllele>,
    pub sample_format: SampleFormatData,
    pub format_data: Option<RecordFormatData>,
}

#[derive(Debug, Clone)]
pub struct RecordFormatData {
    pub format_order: Vec<Vec<u8>>,
    pub sample_gts: Vec<Vec<GenotypeAllele>>,
    pub fields: Vec<FormatField>,
}

#[derive(Debug, Clone)]
pub struct FormatField {
    pub tag: Vec<u8>,
    pub values: FormatFieldValues,
}

#[derive(Debug, Clone)]
pub enum FormatFieldValues {
    Integer(Vec<Vec<i32>>),
    Float(Vec<Vec<f32>>),
    String(Vec<Vec<u8>>),
}

#[derive(Debug, Clone)]
pub enum SampleFormatData {
    Sv(SvSampleFormatData),
    Cnv(CnvSampleFormatData),
}

#[derive(Debug, Clone)]
pub struct SvSampleFormatData {
    pub gq: i32,
    pub pl: Vec<i32>,
    pub ad: Vec<i32>,
}

#[derive(Debug, Clone)]
pub struct CnvSampleFormatData {
    pub cn: f32,
    pub cnq: f32,
}

impl VcfWriteData {
    pub fn sv_format(&self) -> Option<&SvSampleFormatData> {
        if let SampleFormatData::Sv(format_data) = &self.sample_format {
            return Some(format_data);
        }
        None
    }

    pub fn cnv_format(&self) -> Option<&CnvSampleFormatData> {
        if let SampleFormatData::Cnv(format_data) = &self.sample_format {
            return Some(format_data);
        }
        None
    }

    pub fn sample_count(&self) -> usize {
        self.format_data
            .as_ref()
            .map_or(1, |format_data| format_data.sample_gts.len())
    }

    pub fn sample_gt(&self, sample_idx: usize) -> Option<&[GenotypeAllele]> {
        if let Some(format_data) = &self.format_data {
            return format_data
                .sample_gts
                .get(sample_idx)
                .map(|genotype| genotype.as_slice());
        }
        if sample_idx == 0 {
            return Some(self.gt.as_slice());
        }
        None
    }

    pub fn format_order(&self) -> Vec<Vec<u8>> {
        if let Some(format_data) = &self.format_data {
            return format_data.format_order.clone();
        }
        match &self.sample_format {
            SampleFormatData::Sv(_) => vec![
                b"GT".to_vec(),
                b"GQ".to_vec(),
                b"PL".to_vec(),
                b"AD".to_vec(),
            ],
            SampleFormatData::Cnv(_) => vec![b"GT".to_vec(), b"CN".to_vec(), b"CNQ".to_vec()],
        }
    }

    pub fn format_field(&self, tag: &[u8]) -> Option<FormatFieldValues> {
        if let Some(format_data) = &self.format_data {
            return format_data
                .fields
                .iter()
                .find(|field| field.tag == tag)
                .map(|field| field.values.clone());
        }

        match (tag, &self.sample_format) {
            (b"GQ", SampleFormatData::Sv(format_data)) => {
                Some(FormatFieldValues::Integer(vec![vec![format_data.gq]]))
            }
            (b"PL", SampleFormatData::Sv(format_data)) => {
                Some(FormatFieldValues::Integer(vec![format_data.pl.clone()]))
            }
            (b"AD", SampleFormatData::Sv(format_data)) => {
                Some(FormatFieldValues::Integer(vec![format_data.ad.clone()]))
            }
            (b"CN", SampleFormatData::Cnv(format_data)) => {
                Some(FormatFieldValues::Float(vec![vec![format_data.cn]]))
            }
            (b"CNQ", SampleFormatData::Cnv(format_data)) => {
                Some(FormatFieldValues::Float(vec![vec![format_data.cnq]]))
            }
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct VariantInternal {
    pub start: f64,
    pub end: f64, // abs(svlen)
    pub interval: Option<[f64; 2]>,
    pub svlen: f64,
    pub vcf_id: usize,
    pub svtype: SvType,
    pub index: usize,
    pub max_dist: f64,
    pub info_hash: i32,
    pub sequence: Option<Vec<u8>>,
    pub support_mask: Vec<u64>,
    pub support_calls: i64,
    pub sample_id: usize,
    pub start_mean: f64,
    pub start_variance: f64,
    pub svlen_mean: f64,
    pub svlen_variance: f64,
    pub id: String,
    pub id_list: Vec<String>,
    pub svclaim: Option<String>,
    pub svclaims: Vec<String>,
    pub bnd: Option<BndData>,
    pub bnd_event: Option<BndEventData>,
    pub trid: Option<TrId>,
    pub vcf: Option<VcfWriteData>,
}

impl fmt::Display for VariantInternal {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "(id: {}, sample: {}, start: {}, end: {})",
            self.id, self.sample_id, self.start, self.end
        )
    }
}
