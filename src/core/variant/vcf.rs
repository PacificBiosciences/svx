use super::{
    BndData, CnvSampleFormatData, SampleFormatData, SvSampleFormatData, VariantInternal,
    VcfWriteData,
};
use crate::{
    cli::MergeArgsInner,
    core::variant::tr::{TrContainmentConfig, TrContainmentQuery, annotate_tr_containment},
    io::bed_reader::TrId,
    utils::util::{MISSING_INTEGER, Result, VECTOR_END_INTEGER, hash_bytes, stable_hash64},
};
use rust_htslib::bcf;

use crate::core::{containers::interval_tree::IntervalTree, svtype::SvType};

// TODO: Currently we assume a variant corresponds to a single sample
impl VariantInternal {
    pub fn from_vcf_record(
        record: &bcf::Record,
        vcf_id: usize,
        contig: &str,
        args: &MergeArgsInner,
        tr_it: &Option<&IntervalTree<u32, TrId>>,
    ) -> Result<Self> {
        let svtype = SvType::from_vcf_record(record)?;

        let sequence = if svtype == SvType::INSERTION {
            Self::insertion_sequence_from_record(record)
        } else {
            None
        };

        let svlen = {
            let svlen = Self::parse_info_i64(record, b"SVLEN")?;
            if svtype == SvType::BND {
                svlen.unwrap_or(0) as f64
            } else {
                svlen
                    .map(|value| value as f64)
                    .ok_or_else(|| crate::svx_error!("SVLEN missing in VCF record"))?
            }
        };

        let start_i64 = record.pos();
        let start = start_i64 as f64;
        let (start, end, bnd) = if svtype == SvType::BND {
            let alleles = record.alleles();
            let alt = alleles
                .get(1)
                .copied()
                .ok_or_else(|| crate::svx_error!("BND record missing ALT allele"))?;

            let (mate_contig, mate_pos0) = Self::parse_bnd_mate_from_record(record, alt)?;
            let strands = Self::parse_bnd_strands_from_record(record, alt)?;

            let mate_ids = record
                .info(b"MATEID")
                .string()
                .map_err(|e| crate::svx_error!("Error reading MATEID INFO from BND record: {e}"))?
                .ok_or_else(|| crate::svx_error!("MATEID missing in BND record"))?;
            let mate_id_raw = mate_ids
                .first()
                .copied()
                .ok_or_else(|| crate::svx_error!("MATEID is empty in BND record"))?;
            let mate_id_raw = std::str::from_utf8(mate_id_raw)
                .map_err(|e| crate::svx_error!("MATEID is not valid UTF-8 in BND record: {e}"))?;
            let mate_id = format!("{vcf_id}_{mate_id_raw}");
            let homlen = Self::parse_info_i64(record, b"HOMLEN")?;
            let event_type = Self::parse_info_string(record, b"EVENTTYPE")?;
            let inv_event_id = if event_type.as_deref() == Some("INV") {
                Self::parse_info_string(record, b"EVENT")?
            } else {
                None
            };
            let has_inv_breakpoint_filter = Self::record_has_filter(record, "InvBreakpoint")?;

            let (merge_start, merge_end) =
                Self::bnd_merge_coords(record, contig, &mate_contig, mate_pos0)?;

            (
                merge_start,
                merge_end,
                Some(BndData {
                    contig: contig.to_string(),
                    pos0: record.pos(),
                    mate_id,
                    mate_contig,
                    mate_pos: mate_pos0,
                    strands,
                    homlen,
                    inv_event_id,
                    has_inv_breakpoint_filter,
                }),
            )
        } else {
            (start, svlen.abs(), None)
        };

        let raw_id_bytes = record.id();
        let id = Self::internal_variant_id(vcf_id, raw_id_bytes.as_slice());

        let is_cn = svtype == SvType::CNV || SvType::id_is_cnv(raw_id_bytes.as_slice());
        let svclaim = Self::parse_info_string(record, b"SVCLAIM")?;

        let mut max_dist: f64 = f64::from(args.max_dist);

        // There may be a per-variant distance threshold (i.e., a VCF is given that comes from svx)
        // There may also be a per-sample distance threshold
        // max_dist = PER_SAMPLE_DISTS[sample]

        // Alternatively, there is a length-based threshold
        if args.use_linear_threshold && args.max_dist_linear > 0.0 {
            max_dist = f64::from(args.max_dist_linear) * svlen.abs();
            // TODO: Some code that will still use max_dist but ONLY if it was explicitly set
            if args.min_dist != -1 {
                max_dist = max_dist.max(f64::from(args.min_dist));
            }
        }

        let interval = match svtype {
            SvType::DELETION | SvType::INVERSION | SvType::DUPLICATION | SvType::CNV => {
                if args.min_recip_overlap > 0.0 {
                    let end_pos = if let Some(end_pos) = Self::parse_info_i64(record, b"END")? {
                        end_pos as f64
                    } else {
                        start + svlen.abs()
                    };
                    Some([start, end_pos])
                } else {
                    None
                }
            }
            _ => None,
        };

        let vcf = {
            let alleles: Vec<Vec<u8>> = record.alleles().iter().map(|s| s.to_vec()).collect();

            let gt_alleles = record
                .genotypes()
                .map_err(|e| crate::svx_error!("Error reading genotypes: {e}"))?
                .get(0)
                .iter()
                .copied()
                .collect::<Vec<_>>();

            let sample_format = if is_cn {
                let cn = Self::parse_required_format_scalar_f32(record, b"CN")?;
                let cnq = Self::parse_required_format_scalar_f32(record, b"CNQ")?;
                SampleFormatData::Cnv(CnvSampleFormatData { cn, cnq })
            } else {
                let gq: i32 = match record.format(b"GQ").integer() {
                    Ok(buffer_backed) => buffer_backed
                        .iter()
                        .next()
                        .and_then(|slice| slice.first())
                        .copied()
                        .unwrap_or(0),
                    Err(_) => 0,
                };

                let pl: Vec<i32> = match record.format(b"PL").integer() {
                    Ok(pl_values_per_sample) => match pl_values_per_sample.first() {
                        Some(pl_slice) if !pl_slice.is_empty() => pl_slice.to_vec(),
                        _ => vec![MISSING_INTEGER, MISSING_INTEGER, MISSING_INTEGER],
                    },
                    Err(_) => vec![MISSING_INTEGER, MISSING_INTEGER, MISSING_INTEGER],
                };

                let ad: Vec<i32> = match record.format(b"AD").integer() {
                    Ok(ad_values_per_sample) => match ad_values_per_sample.first() {
                        Some(ad_slice) if !ad_slice.is_empty() => {
                            if ad_slice.len() == 1 && ad_slice[0] == MISSING_INTEGER {
                                vec![MISSING_INTEGER, VECTOR_END_INTEGER]
                            } else {
                                ad_slice.to_vec()
                            }
                        }
                        _ => vec![MISSING_INTEGER, VECTOR_END_INTEGER],
                    },
                    Err(_) => vec![MISSING_INTEGER, VECTOR_END_INTEGER],
                };

                SampleFormatData::Sv(SvSampleFormatData { gq, pl, ad })
            };

            Some(VcfWriteData {
                rid: record.rid(),
                pos: record.pos(),
                alleles,
                gt: gt_alleles,
                sample_format,
            })
        };

        let trid = if tr_it.is_some() && matches!(svtype, SvType::INSERTION | SvType::DELETION) {
            let containment_query = {
                let to_u32 = |label: &str, value: i64| -> Result<u32> {
                    u32::try_from(value).map_err(|_| {
                        crate::svx_error!(
                            "VCF {} coordinate out of range for TR containment: {}",
                            label,
                            value
                        )
                    })
                };

                match svtype {
                    SvType::INSERTION => TrContainmentQuery::Insertion {
                        pos: to_u32("POS", start_i64.saturating_add(1))?,
                    },
                    SvType::DELETION => {
                        let end_i64 = if let Some(end_pos) = Self::parse_info_i64(record, b"END")? {
                            end_pos
                        } else {
                            start_i64.saturating_add(svlen.abs() as i64)
                        };
                        TrContainmentQuery::Deletion {
                            start: to_u32("POS", start_i64)?,
                            end: to_u32("END", end_i64)?,
                        }
                    }
                    _ => unreachable!(),
                }
            };
            let containment_cfg = TrContainmentConfig {
                span_query_slop: args.tr_span_query_slop,
                min_span_containment_scaled: args.tr_min_span_containment,
                min_span_overlap_bp: args.tr_min_span_overlap_bp,
                ins_max_dist: args.tr_ins_max_dist,
            };
            annotate_tr_containment(containment_query, tr_it, containment_cfg)
        } else {
            None
        };

        log::trace!(
            "Variant: {}, start: {}, end: {}, svlen: {}, svtype: {}, is_cn: {}",
            String::from_utf8_lossy(raw_id_bytes.as_slice()),
            start,
            end,
            svlen,
            svtype,
            is_cn
        );

        let info_hash = Self::info_hash_from_record(record);

        Ok(Self {
            start,
            end,
            interval,
            svlen,
            vcf_id,
            sample_id: vcf_id,
            svtype,
            id,
            svclaim,
            bnd,
            bnd_event: None,
            index: 0,
            max_dist,
            info_hash,
            trid,
            sequence,
            vcf,
        })
    }

    fn internal_variant_id(vcf_id: usize, raw_id_bytes: &[u8]) -> String {
        if let Ok(vcf_id_str) = std::str::from_utf8(raw_id_bytes) {
            return format!("{vcf_id}_{vcf_id_str}");
        }

        let stable_hash = stable_hash64(raw_id_bytes);
        let secondary_hash = hash_bytes(raw_id_bytes) as u32;
        format!(
            "{vcf_id}__svx_nonutf8_{stable_hash:016x}_{secondary_hash:08x}_{}",
            raw_id_bytes.len()
        )
    }

    fn format_record_vcf_line(record: &bcf::Record, context: &str) -> Result<String> {
        let formatted =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| record.to_vcf_string()));

        match formatted {
            Ok(Ok(vcf_line)) => Ok(vcf_line),
            Ok(Err(e)) => Err(crate::svx_error!(
                "Failed to format record as VCF while {context}: {e}"
            )),
            Err(_) => Err(crate::svx_error!(
                "Failed to format record as VCF while {context}: rust-htslib panicked"
            )),
        }
    }

    pub(super) fn info_hash_from_record(record: &bcf::Record) -> i32 {
        let vcf_line = match Self::format_record_vcf_line(record, "hashing the INFO field") {
            Ok(s) => s,
            Err(e) => {
                log::debug!("{e}");
                return 0;
            }
        };

        Self::info_hash_from_vcf_line(&vcf_line)
    }

    pub(super) fn info_hash_from_vcf_line(vcf_line: &str) -> i32 {
        let bytes = vcf_line.as_bytes();
        let mut info_start = None;
        let mut info_end = None;
        let mut tab_count = 0usize;
        for (i, &b) in bytes.iter().enumerate() {
            if b != b'\t' {
                continue;
            }
            tab_count += 1;
            if tab_count == 7 {
                info_start = Some(i + 1);
            } else if tab_count == 8 {
                info_end = Some(i);
                break;
            }
        }

        let Some(start) = info_start else {
            log::debug!("Failed to locate INFO field in formatted VCF line for hashing.");
            return 0;
        };

        let mut end = info_end.unwrap_or(bytes.len());
        while end > start && (bytes[end - 1] == b'\n' || bytes[end - 1] == b'\r') {
            end -= 1;
        }

        hash_bytes(&bytes[start..end])
    }

    pub(super) fn parse_info_string(record: &bcf::Record, key: &[u8]) -> Result<Option<String>> {
        let raw = match record.info(key).string() {
            Ok(v) => v,
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("undefined in BCF/VCF header") {
                    return Ok(None);
                }
                return Err(crate::svx_error!(
                    "Error reading {} INFO from record: {e}",
                    String::from_utf8_lossy(key)
                ));
            }
        };
        let Some(values) = raw else {
            return Ok(None);
        };
        let Some(v) = values.first().copied() else {
            return Ok(None);
        };
        if v.is_empty() {
            return Ok(None);
        }
        let s = std::str::from_utf8(v).map_err(|e| {
            crate::svx_error!(
                "{} INFO is not valid UTF-8: {e}",
                String::from_utf8_lossy(key)
            )
        })?;
        Ok(Some(s.to_string()))
    }

    pub(super) fn record_has_filter(record: &bcf::Record, filter_id: &str) -> Result<bool> {
        let vcf_line = Self::format_record_vcf_line(record, "parsing the FILTER field")?;
        let filter_field = vcf_line
            .split('\t')
            .nth(6)
            .ok_or_else(|| crate::svx_error!("Failed to locate FILTER field in VCF record"))?;
        let filter_field = filter_field.trim_end_matches(['\r', '\n']);
        if filter_field == "." {
            return Ok(false);
        }
        Ok(filter_field.split(';').any(|value| value == filter_id))
    }

    pub(super) fn parse_info_i64(record: &bcf::Record, key: &[u8]) -> Result<Option<i64>> {
        let raw = match record.info(key).integer() {
            Ok(v) => v,
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("undefined in BCF/VCF header") {
                    return Ok(None);
                }
                return Err(crate::svx_error!(
                    "Error reading {} INFO from record: {e}",
                    String::from_utf8_lossy(key)
                ));
            }
        };
        let Some(values) = raw else {
            return Ok(None);
        };
        let Some(v) = values.first().copied() else {
            return Ok(None);
        };
        if v == MISSING_INTEGER || v == VECTOR_END_INTEGER {
            return Ok(None);
        }
        Ok(Some(i64::from(v)))
    }

    pub(super) fn parse_info_f64(record: &bcf::Record, key: &[u8]) -> Result<Option<f64>> {
        let s = Self::parse_info_string(record, key)?;
        let Some(s) = s else {
            return Ok(None);
        };
        let v: f64 = s.parse().map_err(|e| {
            crate::svx_error!(
                "{} INFO value {s:?} is not a valid float: {e}",
                String::from_utf8_lossy(key)
            )
        })?;
        Ok(Some(v))
    }

    fn insertion_sequence_from_record(record: &bcf::Record) -> Option<Vec<u8>> {
        let alleles = record.alleles();
        let alt = alleles.get(1).copied()?;
        if Self::is_symbolic_or_nonliteral_alt(alt) {
            return None;
        }
        Some(alt.to_vec())
    }

    fn is_symbolic_or_nonliteral_alt(alt: &[u8]) -> bool {
        alt.is_empty()
            || alt == b"."
            || alt == b"*"
            || (alt.first() == Some(&b'<') && alt.last() == Some(&b'>'))
            || alt.contains(&b'[')
            || alt.contains(&b']')
    }

    fn parse_required_format_scalar_f32(record: &bcf::Record, key: &[u8]) -> Result<f32> {
        let value = Self::parse_format_scalar_f32(record, key)?;
        let Some(value) = value else {
            return Err(crate::svx_error!(
                "Missing {} FORMAT value in VCF record",
                String::from_utf8_lossy(key)
            ));
        };
        Ok(value)
    }

    fn parse_format_scalar_f32(record: &bcf::Record, key: &[u8]) -> Result<Option<f32>> {
        if let Ok(values_per_sample) = record.format(key).float() {
            if let Some(sample_values) = values_per_sample.first() {
                if let Some(value) = sample_values.first() {
                    if value.is_nan() {
                        return Ok(None);
                    }
                    return Ok(Some(*value));
                }
            }
            return Ok(None);
        }

        if let Ok(values_per_sample) = record.format(key).integer() {
            if let Some(sample_values) = values_per_sample.first() {
                if let Some(value) = sample_values.first() {
                    if *value == MISSING_INTEGER || *value == VECTOR_END_INTEGER {
                        return Ok(None);
                    }
                    return Ok(Some(*value as f32));
                }
            }
            return Ok(None);
        }

        Err(crate::svx_error!(
            "Error reading {} FORMAT from VCF record as float or integer",
            String::from_utf8_lossy(key)
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::containers::interval_tree::{Interval, IntervalTree};
    use crate::core::variant::test_utils::make_temp_vcf;
    use crate::io::bed_reader::TrId;
    use crate::utils::util::init_logger;
    use rust_htslib::bcf::{self, Read};

    #[test]
    fn test_info_hash_is_hash_of_formatted_vcf_info_field() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=X,Number=1,Type=Integer,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1\tv1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=10;X=1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let vcf_line = record.to_vcf_string().unwrap();
        let bytes = vcf_line.as_bytes();
        let mut info_start = None;
        let mut info_end = None;
        let mut tab_count = 0usize;
        for (i, &b) in bytes.iter().enumerate() {
            if b != b'\t' {
                continue;
            }
            tab_count += 1;
            if tab_count == 7 {
                info_start = Some(i + 1);
            } else if tab_count == 8 {
                info_end = Some(i);
                break;
            }
        }
        let start = info_start.expect("expected to find INFO field start");
        let mut end = info_end.unwrap_or(bytes.len());
        while end > start && (bytes[end - 1] == b'\n' || bytes[end - 1] == b'\r') {
            end -= 1;
        }

        let expected = hash_bytes(&bytes[start..end]);
        let got = VariantInternal::info_hash_from_record(&record);
        assert_eq!(got, expected);
    }

    #[test]
    fn cnv_record_parses_cn_and_cnq_format() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10\tsawfish:CNV:0:0:10\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:CN:CNQ\t0/1:1:8
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap();
        let format_data = variant
            .vcf
            .as_ref()
            .and_then(VcfWriteData::cnv_format)
            .expect("expected CNV format payload");
        assert_eq!(format_data.cn, 1.0);
        assert_eq!(format_data.cnq, 8.0);
    }

    #[test]
    fn cnv_record_parses_svclaim_info() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SVCLAIM,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10\tsawfish:CNV:0:0:10\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;SVCLAIM=D\tGT:CN:CNQ\t0/1:1:8
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap();
        assert_eq!(variant.svclaim.as_deref(), Some("D"));
    }

    #[test]
    fn non_cnv_del_record_parses_svclaim_info() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SVCLAIM,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;SVCLAIM=DJ\tGT:GQ:PL:AD\t0/1:20:0,20,200:10,8
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap();
        assert_eq!(variant.svclaim.as_deref(), Some("DJ"));
    }

    #[test]
    fn cnv_record_errors_when_cnq_missing() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10\tsawfish:CNV:0:0:10\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:CN:CNQ\t0/1:1:.
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let err = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap_err();
        assert!(
            err.to_string()
                .contains("Missing CNQ FORMAT value in VCF record"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn tr_annotation_is_limited_to_ins_and_del() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t101\tins1\tN\tNAAAA\t.\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1
chr1\t102\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10\tGT\t0/1
chr1\t103\tinv1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=10\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let args = MergeArgsInner::default();

        let tr_tree = IntervalTree::new(vec![Interval::new(
            90,
            200,
            TrId {
                id: "TR1".to_string(),
                motif_len: 2,
            },
        )]);
        let tr_it = Some(&tr_tree);

        let mut records = reader.records();
        let ins_record = records.next().unwrap().unwrap();
        let del_record = records.next().unwrap().unwrap();
        let inv_record = records.next().unwrap().unwrap();

        let ins_variant =
            VariantInternal::from_vcf_record(&ins_record, 0, "chr1", &args, &tr_it).unwrap();
        let del_variant =
            VariantInternal::from_vcf_record(&del_record, 0, "chr1", &args, &tr_it).unwrap();
        let inv_variant =
            VariantInternal::from_vcf_record(&inv_record, 0, "chr1", &args, &tr_it).unwrap();

        assert!(ins_variant.trid.is_some());
        assert!(del_variant.trid.is_some());
        assert!(inv_variant.trid.is_none());
    }

    #[test]
    fn ins_symbolic_alt_does_not_set_sequence() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t101\tins_symbolic\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap();
        assert!(variant.sequence.is_none());
    }

    #[test]
    fn ins_symbolic_subtype_alt_does_not_set_sequence() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t101\tins_symbolic\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap();
        assert!(variant.sequence.is_none());
    }

    #[test]
    fn ins_literal_alt_sets_sequence() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t101\tins_literal\tN\tNAAAA\t.\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap();
        assert_eq!(variant.sequence, Some(b"NAAAA".to_vec()));
    }

    #[test]
    fn svlen_type_error_is_propagated_for_non_bnd_records() {
        init_logger();

        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t101\tins_bad\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=oops\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let err = VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None).unwrap_err();
        let err_s = err.to_string();
        assert!(
            err_s.contains("Error reading SVLEN INFO from record"),
            "unexpected error: {err_s}"
        );
        assert!(
            !err_s.contains("SVLEN missing in VCF record"),
            "SVLEN type error was incorrectly treated as missing: {err_s}"
        );
    }
}
