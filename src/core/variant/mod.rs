mod bnd;
mod similarity;
mod tr;
mod vcf;

#[cfg(test)]
pub(crate) mod test_utils;

use crate::core::svtype::SvType;
use crate::io::bed_reader::TrId;
use rust_htslib::bcf::record::GenotypeAllele;
use std::fmt;

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
}

// TODO: Refactor
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
    pub sample_id: usize, // TODO: This is not used for now, currently assumes we have 1 sample per VCF
    pub id: String,
    pub svclaim: Option<String>,
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
