use super::{
    BndData, CnvSampleFormatData, FormatField, FormatFieldValues, HeaderFormatCache,
    RecordFormatData, SampleFormatData, SvSampleFormatData, VariantInternal, VcfWriteData,
};
use crate::{
    cli::MergeArgsInner,
    core::variant::tr::{TrContainmentConfig, TrContainmentQuery, annotate_tr_containment},
    io::{bed_reader::TrId, vcf_reader::SampleMapping},
    utils::util::{MISSING_INTEGER, Result, VECTOR_END_INTEGER, stable_hash, stable_hash64},
};
use rust_htslib::{
    bcf::{
        self,
        header::{HeaderRecord, TagType},
    },
    htslib,
};
use std::{
    collections::{BTreeSet, HashMap},
    ffi::CString,
};

use crate::core::{containers::interval_tree::IntervalTree, svtype::SvType};

impl VariantInternal {
    pub fn from_vcf_record(
        record: &bcf::Record,
        vcf_id: usize,
        contig: &str,
        args: &MergeArgsInner,
        tr_it: &Option<&IntervalTree<u32, TrId>>,
        format_lookup: &HeaderFormatCache,
        local_to_merged: Option<&[usize]>,
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
            let mate_id =
                Self::internal_variant_id(vcf_id, mate_id_raw.as_bytes(), local_to_merged);
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
        let id = Self::internal_variant_id(vcf_id, raw_id_bytes.as_slice(), local_to_merged);
        let id_list = {
            let id_list_values = Self::parse_info_string_values(record, b"IDLIST")?;
            if id_list_values.is_empty() {
                vec![id.clone()]
            } else {
                id_list_values
                    .iter()
                    .map(|raw_id| {
                        Self::internal_variant_id(vcf_id, raw_id.as_slice(), local_to_merged)
                    })
                    .collect()
            }
        };
        let support_calls = {
            let support_calls = Self::parse_info_i64(record, b"SUPP_CALLS")?.unwrap_or(1);
            if support_calls < 1 {
                return Err(crate::svx_error!(
                    "SUPP_CALLS={} is invalid in record {}; expected a positive integer",
                    support_calls,
                    String::from_utf8_lossy(raw_id_bytes.as_slice())
                ));
            }
            support_calls
        };
        let parse_ci_offsets = |values: &[i64], key: &[u8]| -> Result<Option<(f64, f64)>> {
            if values.is_empty() {
                return Ok(None);
            }
            if values.len() != 2 {
                return Err(crate::svx_error!(
                    "{} must contain exactly two integer offsets in record {}",
                    String::from_utf8_lossy(key),
                    String::from_utf8_lossy(raw_id_bytes.as_slice())
                ));
            }
            let lower = values[0];
            let upper = values[1];
            if lower > upper {
                return Err(crate::svx_error!(
                    "{} has invalid bounds [{lower}, {upper}] in record {}; expected lower <= upper",
                    String::from_utf8_lossy(key),
                    String::from_utf8_lossy(raw_id_bytes.as_slice())
                ));
            }

            let center_offset = (lower as f64 + upper as f64) / 2.0;
            let half_width = (upper - lower) as f64 / 2.0;
            let sigma = half_width / 1.96;
            let variance = sigma * sigma;

            if !variance.is_finite() {
                return Err(crate::svx_error!(
                    "{} produced non-finite variance in record {}",
                    String::from_utf8_lossy(key),
                    String::from_utf8_lossy(raw_id_bytes.as_slice())
                ));
            }

            Ok(Some((center_offset, variance.max(0.0))))
        };
        let cipos = parse_ci_offsets(&Self::parse_info_i64_values(record, b"CIPOS")?, b"CIPOS")?;
        let ciend = parse_ci_offsets(&Self::parse_info_i64_values(record, b"CIEND")?, b"CIEND")?;

        let (start_mean, start_variance, svlen_mean, svlen_variance) = {
            let default_start_mean = start + 1.0;
            let (start_mean, start_variance) = if let Some((center_offset, variance)) = cipos {
                (default_start_mean + center_offset, variance)
            } else {
                (default_start_mean, 0.0)
            };

            if svtype == SvType::BND {
                (start_mean, start_variance, svlen, 0.0)
            } else {
                let interval_like = matches!(
                    svtype,
                    SvType::DELETION | SvType::INVERSION | SvType::DUPLICATION | SvType::CNV
                );
                if !interval_like {
                    if ciend.is_some() {
                        return Err(crate::svx_error!(
                            "CIEND is unsupported for {} in record {}",
                            svtype,
                            String::from_utf8_lossy(raw_id_bytes.as_slice())
                        ));
                    }
                    (start_mean, start_variance, svlen, 0.0)
                } else {
                    match (cipos, ciend) {
                        (None, None) => (start_mean, start_variance, svlen, 0.0),
                        (Some(_), Some((end_center_offset, end_variance)))
                        | (None, Some((end_center_offset, end_variance))) => {
                            let end_anchor_1based = if let Some(end_from_info) =
                                Self::parse_info_i64(record, b"END")?
                            {
                                end_from_info as f64
                            } else {
                                default_start_mean + svlen.abs()
                            };
                            let end_mean = end_anchor_1based + end_center_offset;
                            let abs_svlen_mean = (end_mean - start_mean).abs();
                            let svlen_sign = if svtype == SvType::DELETION {
                                -1.0
                            } else {
                                1.0
                            };
                            (
                                start_mean,
                                start_variance,
                                svlen_sign * abs_svlen_mean,
                                (start_variance + end_variance).max(0.0),
                            )
                        }
                        (Some(_), None) => {
                            return Err(crate::svx_error!(
                                "Record {} has incomplete CI fields; expected CIEND when CIPOS is present",
                                String::from_utf8_lossy(raw_id_bytes.as_slice())
                            ));
                        }
                    }
                }
            }
        };

        let is_cn = svtype == SvType::CNV || SvType::id_is_cnv(raw_id_bytes.as_slice());
        let svclaim = Self::parse_info_string(record, b"SVCLAIM")?;
        let svclaims = {
            let mut merged_claims = BTreeSet::new();
            if let Some(claim) = svclaim.as_ref() {
                merged_claims.insert(claim.clone());
            }
            for claim in Self::parse_info_utf8_values(record, b"SVCLAIM_SET")? {
                merged_claims.insert(claim);
            }
            merged_claims.into_iter().collect::<Vec<_>>()
        };

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

        let (vcf, support_mask) = {
            let alleles: Vec<Vec<u8>> = record.alleles().iter().map(|s| s.to_vec()).collect();
            let sample_count = record.sample_count() as usize;
            if sample_count == 0 {
                return Err(crate::svx_error!(
                    "Record {} has no samples; sites-only records are unsupported",
                    String::from_utf8_lossy(raw_id_bytes.as_slice())
                ));
            }

            let parsed_genotypes = record
                .genotypes()
                .map_err(|e| crate::svx_error!("Error reading genotypes: {e}"))?;
            let sample_gts = (0..sample_count)
                .map(|sample_idx| {
                    parsed_genotypes
                        .get(sample_idx)
                        .iter()
                        .copied()
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();

            let gt_alleles = sample_gts.first().cloned().ok_or_else(|| {
                crate::svx_error!(
                    "Record {} is missing genotype data for sample 0",
                    String::from_utf8_lossy(raw_id_bytes.as_slice())
                )
            })?;

            let format_order = Self::parse_format_order_from_record(record, format_lookup)?;
            let mut format_fields = Vec::new();
            for tag in &format_order {
                if tag.as_slice() == b"GT" {
                    continue;
                }

                let tag_type = format_lookup.type_by_tag.get(tag).copied().ok_or_else(|| {
                    crate::svx_error!(
                        "FORMAT tag {} is undefined or malformed in header for record {}",
                        String::from_utf8_lossy(tag),
                        String::from_utf8_lossy(raw_id_bytes.as_slice())
                    )
                })?;

                let values = match tag_type {
                    TagType::Integer => {
                        let values_per_sample =
                            record.format(tag.as_slice()).integer().map_err(|e| {
                                crate::svx_error!(
                                    "Error reading {} FORMAT: {e}",
                                    String::from_utf8_lossy(tag)
                                )
                            })?;
                        if values_per_sample.len() != sample_count {
                            return Err(crate::svx_error!(
                                "Malformed {} FORMAT in record {}: expected {} samples, got {}",
                                String::from_utf8_lossy(tag),
                                String::from_utf8_lossy(raw_id_bytes.as_slice()),
                                sample_count,
                                values_per_sample.len()
                            ));
                        }
                        FormatFieldValues::Integer(
                            values_per_sample
                                .iter()
                                .map(|values| values.to_vec())
                                .collect(),
                        )
                    }
                    TagType::Float => {
                        let values_per_sample =
                            record.format(tag.as_slice()).float().map_err(|e| {
                                crate::svx_error!(
                                    "Error reading {} FORMAT: {e}",
                                    String::from_utf8_lossy(tag)
                                )
                            })?;
                        if values_per_sample.len() != sample_count {
                            return Err(crate::svx_error!(
                                "Malformed {} FORMAT in record {}: expected {} samples, got {}",
                                String::from_utf8_lossy(tag),
                                String::from_utf8_lossy(raw_id_bytes.as_slice()),
                                sample_count,
                                values_per_sample.len()
                            ));
                        }
                        FormatFieldValues::Float(
                            values_per_sample
                                .iter()
                                .map(|values| values.to_vec())
                                .collect(),
                        )
                    }
                    TagType::String => {
                        let values_per_sample =
                            record.format(tag.as_slice()).string().map_err(|e| {
                                crate::svx_error!(
                                    "Error reading {} FORMAT: {e}",
                                    String::from_utf8_lossy(tag)
                                )
                            })?;
                        if values_per_sample.len() != sample_count {
                            return Err(crate::svx_error!(
                                "Malformed {} FORMAT in record {}: expected {} samples, got {}",
                                String::from_utf8_lossy(tag),
                                String::from_utf8_lossy(raw_id_bytes.as_slice()),
                                sample_count,
                                values_per_sample.len()
                            ));
                        }
                        FormatFieldValues::String(
                            values_per_sample
                                .iter()
                                .map(|values| values.to_vec())
                                .collect(),
                        )
                    }
                    TagType::Flag => {
                        return Err(crate::svx_error!(
                            "Unsupported FORMAT type Flag for tag {} in record {}",
                            String::from_utf8_lossy(tag),
                            String::from_utf8_lossy(raw_id_bytes.as_slice())
                        ));
                    }
                };

                format_fields.push(FormatField {
                    tag: tag.clone(),
                    values,
                });
            }

            let mut support_mask = vec![0u64; sample_count.div_ceil(u64::BITS as usize)];
            let mut first_supporting_sample_idx = None;
            for (sample_idx, genotype) in sample_gts.iter().enumerate() {
                if genotype
                    .iter()
                    .any(|allele| allele.index().is_some_and(|idx| idx > 0))
                {
                    let word = sample_idx / (u64::BITS as usize);
                    let bit = sample_idx % (u64::BITS as usize);
                    support_mask[word] |= 1u64 << bit;
                    if first_supporting_sample_idx.is_none() {
                        first_supporting_sample_idx = Some(sample_idx);
                    }
                }
            }

            let sample_format = if is_cn {
                let cn_sample_idx = first_supporting_sample_idx.unwrap_or(0);
                let cn = Self::parse_required_format_scalar_f32_for_sample(
                    &format_fields,
                    b"CN",
                    cn_sample_idx,
                )?;
                let cnq = Self::parse_required_format_scalar_f32_for_sample(
                    &format_fields,
                    b"CNQ",
                    cn_sample_idx,
                )?;
                SampleFormatData::Cnv(CnvSampleFormatData { cn, cnq })
            } else {
                let gq: i32 = format_fields
                    .iter()
                    .find(|field| field.tag == b"GQ")
                    .and_then(|field| match &field.values {
                        FormatFieldValues::Integer(values_per_sample) => values_per_sample
                            .first()
                            .and_then(|sample_values| sample_values.first())
                            .copied(),
                        _ => None,
                    })
                    .unwrap_or(0);

                let pl: Vec<i32> = format_fields
                    .iter()
                    .find(|field| field.tag == b"PL")
                    .and_then(|field| match &field.values {
                        FormatFieldValues::Integer(values_per_sample) => values_per_sample
                            .first()
                            .filter(|sample_values| !sample_values.is_empty())
                            .cloned(),
                        _ => None,
                    })
                    .unwrap_or_else(|| vec![MISSING_INTEGER, MISSING_INTEGER, MISSING_INTEGER]);

                let ad: Vec<i32> = format_fields
                    .iter()
                    .find(|field| field.tag == b"AD")
                    .and_then(|field| match &field.values {
                        FormatFieldValues::Integer(values_per_sample) => values_per_sample
                            .first()
                            .filter(|sample_values| !sample_values.is_empty())
                            .cloned(),
                        _ => None,
                    })
                    .map(|sample_values| {
                        if sample_values.len() == 1 && sample_values[0] == MISSING_INTEGER {
                            vec![MISSING_INTEGER, VECTOR_END_INTEGER]
                        } else {
                            sample_values
                        }
                    })
                    .unwrap_or_else(|| vec![MISSING_INTEGER, VECTOR_END_INTEGER]);

                SampleFormatData::Sv(SvSampleFormatData { gq, pl, ad })
            };

            (
                Some(VcfWriteData {
                    rid: record.rid(),
                    pos: record.pos(),
                    alleles,
                    gt: gt_alleles,
                    sample_format,
                    format_data: Some(RecordFormatData {
                        format_order,
                        sample_gts,
                        fields: format_fields,
                    }),
                }),
                support_mask,
            )
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
            id_list,
            svclaim,
            svclaims,
            bnd,
            bnd_event: None,
            index: 0,
            max_dist,
            info_hash,
            trid,
            sequence,
            support_mask,
            support_calls,
            start_mean,
            start_variance,
            svlen_mean,
            svlen_variance,
            vcf,
        })
    }

    pub fn remap_support_mask_to_output(
        &mut self,
        sample_mapping: &SampleMapping,
        total_samples: usize,
    ) -> Result<()> {
        let local_sample_count = sample_mapping
            .index_map
            .keys()
            .filter(|(vcf_idx, _)| *vcf_idx == self.vcf_id)
            .count();
        let mut local_to_merged = vec![usize::MAX; local_sample_count];
        for (local_sample_idx, merged_slot) in local_to_merged
            .iter_mut()
            .enumerate()
            .take(local_sample_count)
        {
            let merged_sample_idx = sample_mapping
                .index_map
                .get(&(self.vcf_id, local_sample_idx))
                .copied()
                .ok_or_else(|| {
                    crate::svx_error!(
                        "Missing sample mapping for VCF {} sample {} while remapping support mask for variant {}",
                        self.vcf_id,
                        local_sample_idx,
                        self.id
                    )
                })?;
            *merged_slot = merged_sample_idx;
        }

        self.remap_support_mask_to_output_with_local_map(&local_to_merged, total_samples)
    }

    pub fn remap_support_mask_to_output_with_local_map(
        &mut self,
        local_to_merged: &[usize],
        total_samples: usize,
    ) -> Result<()> {
        let mut remapped = vec![0u64; total_samples.div_ceil(u64::BITS as usize)];
        for (word_idx, word) in self.support_mask.iter().copied().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let bit_idx = remaining.trailing_zeros() as usize;
                let local_sample_idx = word_idx * (u64::BITS as usize) + bit_idx;
                let merged_sample_idx = local_to_merged
                    .get(local_sample_idx)
                    .copied()
                    .ok_or_else(|| {
                        crate::svx_error!(
                            "Missing local sample mapping for VCF {} sample {} while remapping support mask for variant {}",
                            self.vcf_id,
                            local_sample_idx,
                            self.id
                        )
                    })?;
                if merged_sample_idx >= total_samples {
                    return Err(crate::svx_error!(
                        "Mapped sample index {} is out of range (total samples: {}) while remapping support mask for variant {}",
                        merged_sample_idx,
                        total_samples,
                        self.id
                    ));
                }
                let merged_word = merged_sample_idx / (u64::BITS as usize);
                let merged_bit = merged_sample_idx % (u64::BITS as usize);
                remapped[merged_word] |= 1u64 << merged_bit;
                remaining &= remaining - 1;
            }
        }
        self.support_mask = remapped;
        Ok(())
    }

    fn remap_sample_prefixed_id(raw_id: &str, local_to_merged: &[usize]) -> Option<String> {
        let (local_sample_idx, suffix) = raw_id.split_once('_')?;
        let local_sample_idx = local_sample_idx.parse::<usize>().ok()?;
        let merged_sample_idx = local_to_merged.get(local_sample_idx).copied()?;
        Some(format!("{merged_sample_idx}_{suffix}"))
    }

    fn internal_variant_id(
        vcf_id: usize,
        raw_id_bytes: &[u8],
        local_to_merged: Option<&[usize]>,
    ) -> String {
        if let Ok(raw_id) = std::str::from_utf8(raw_id_bytes) {
            if let Some(remapped_id) = local_to_merged
                .and_then(|local_to_merged| Self::remap_sample_prefixed_id(raw_id, local_to_merged))
            {
                return remapped_id;
            }
            return format!("{vcf_id}_{raw_id}");
        }

        let stable_hash64 = stable_hash64(raw_id_bytes);
        let secondary_hash = stable_hash(raw_id_bytes) as u32;
        format!(
            "{vcf_id}__svx_nonutf8_{stable_hash64:016x}_{secondary_hash:08x}_{}",
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

    pub fn parse_info_string_values(record: &bcf::Record, key: &[u8]) -> Result<Vec<Vec<u8>>> {
        let raw = match record.info(key).string() {
            Ok(v) => v,
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("undefined in BCF/VCF header") {
                    return Ok(Vec::new());
                }
                return Err(crate::svx_error!(
                    "Error reading {} INFO from record: {e}",
                    String::from_utf8_lossy(key)
                ));
            }
        };
        let Some(values) = raw else {
            return Ok(Vec::new());
        };

        Ok(values
            .iter()
            .filter(|value| !value.is_empty() && **value != b".")
            .map(|value| value.to_vec())
            .collect())
    }

    fn parse_info_utf8_values(record: &bcf::Record, key: &[u8]) -> Result<Vec<String>> {
        Self::parse_info_string_values(record, key)?
            .iter()
            .map(|value| {
                std::str::from_utf8(value)
                    .map(str::to_owned)
                    .map_err(|error| {
                        crate::svx_error!(
                            "{} INFO entry is not valid UTF-8: {error}",
                            String::from_utf8_lossy(key)
                        )
                    })
            })
            .collect()
    }

    fn parse_header_tag_type(type_name: &str, tag: &str) -> Result<TagType> {
        match type_name {
            "Integer" => Ok(TagType::Integer),
            "Float" => Ok(TagType::Float),
            "String" => Ok(TagType::String),
            "Flag" => Ok(TagType::Flag),
            _ => Err(crate::svx_error!(
                "Unsupported FORMAT type {} for tag {} in header",
                type_name,
                tag
            )),
        }
    }

    pub fn build_header_format_cache(
        header: &bcf::header::HeaderView,
    ) -> Result<HeaderFormatCache> {
        let mut tag_by_id = HashMap::new();
        let mut type_by_tag = HashMap::new();
        for record in header.header_records() {
            let HeaderRecord::Format { values, .. } = record else {
                continue;
            };

            let tag = values.get("ID").ok_or_else(|| {
                crate::svx_error!("FORMAT header record is missing required ID field")
            })?;
            let type_name = values.get("Type").ok_or_else(|| {
                crate::svx_error!(
                    "FORMAT header record for {} is missing required Type field",
                    tag
                )
            })?;
            let tag_type = Self::parse_header_tag_type(type_name, tag)?;
            let tag_bytes = tag.as_bytes().to_vec();
            let tag_cstr = CString::new(tag.as_bytes())
                .map_err(|e| crate::svx_error!("Invalid FORMAT tag {} in header: {}", tag, e))?;
            let tag_id = unsafe {
                htslib::bcf_hdr_id2int(header.inner, htslib::BCF_DT_ID as i32, tag_cstr.as_ptr())
            };
            if tag_id < 0 {
                return Err(crate::svx_error!(
                    "FORMAT tag {} could not be resolved to a header ID",
                    tag
                ));
            }

            if tag_by_id.insert(tag_id, tag_bytes.clone()).is_some() {
                return Err(crate::svx_error!(
                    "Duplicate FORMAT tag ID {} encountered while caching header FORMAT definitions",
                    tag_id
                ));
            }
            if type_by_tag.insert(tag_bytes.clone(), tag_type).is_some() {
                return Err(crate::svx_error!(
                    "Duplicate FORMAT tag {} encountered in header",
                    tag
                ));
            }
        }

        Ok(HeaderFormatCache {
            tag_by_id,
            type_by_tag,
        })
    }

    fn parse_format_order_from_record(
        record: &bcf::Record,
        format_cache: &HeaderFormatCache,
    ) -> Result<Vec<Vec<u8>>> {
        unsafe {
            htslib::bcf_unpack(record.inner, htslib::BCF_UN_FMT as i32);
        }
        let format_count = record.inner().n_fmt() as usize;
        if format_count == 0 {
            return Err(crate::svx_error!(
                "Missing FORMAT column in VCF record; sites-only records are unsupported"
            ));
        }
        let format_ptr = record.inner().d.fmt;
        if format_ptr.is_null() {
            return Err(crate::svx_error!(
                "FORMAT column is malformed in VCF record: missing decoded FORMAT entries"
            ));
        }

        let mut format_order = Vec::with_capacity(format_count);
        for format_idx in 0..format_count {
            let format_entry = unsafe { format_ptr.add(format_idx).read() };
            let format_tag = format_cache
                .tag_by_id
                .get(&format_entry.id)
                .ok_or_else(|| {
                    crate::svx_error!(
                        "FORMAT ID {} in record is undefined in header",
                        format_entry.id
                    )
                })?
                .clone();
            format_order.push(format_tag);
        }

        if !format_order.iter().any(|tag| tag.as_slice() == b"GT") {
            return Err(crate::svx_error!(
                "FORMAT column is missing required GT field"
            ));
        }
        Ok(format_order)
    }

    pub fn info_hash_from_record(record: &bcf::Record) -> i32 {
        let vcf_line = match Self::format_record_vcf_line(record, "hashing the INFO field") {
            Ok(s) => s,
            Err(e) => {
                log::debug!("{e}");
                return 0;
            }
        };

        Self::info_hash_from_vcf_line(&vcf_line)
    }

    pub fn info_hash_from_vcf_line(vcf_line: &str) -> i32 {
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

        stable_hash(&bytes[start..end])
    }

    pub fn parse_info_string(record: &bcf::Record, key: &[u8]) -> Result<Option<String>> {
        let values = Self::parse_info_string_values(record, key)?;
        let Some(v) = values.first() else {
            return Ok(None);
        };
        let s = std::str::from_utf8(v).map_err(|e| {
            crate::svx_error!(
                "{} INFO is not valid UTF-8: {e}",
                String::from_utf8_lossy(key)
            )
        })?;
        Ok(Some(s.to_string()))
    }

    pub fn record_has_filter(record: &bcf::Record, filter_id: &str) -> Result<bool> {
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

    pub fn parse_info_i64(record: &bcf::Record, key: &[u8]) -> Result<Option<i64>> {
        Ok(Self::parse_info_i64_values(record, key)?.into_iter().next())
    }

    pub fn parse_info_i64_values(record: &bcf::Record, key: &[u8]) -> Result<Vec<i64>> {
        let raw = match record.info(key).integer() {
            Ok(v) => v,
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("undefined in BCF/VCF header") {
                    return Ok(Vec::new());
                }
                return Err(crate::svx_error!(
                    "Error reading {} INFO from record: {e}",
                    String::from_utf8_lossy(key)
                ));
            }
        };
        let Some(values) = raw else {
            return Ok(Vec::new());
        };
        let mut parsed_values = Vec::with_capacity(values.len());
        for value in values.iter().copied() {
            if value == MISSING_INTEGER || value == VECTOR_END_INTEGER {
                continue;
            }
            parsed_values.push(i64::from(value));
        }
        Ok(parsed_values)
    }

    pub fn parse_info_f64(record: &bcf::Record, key: &[u8]) -> Result<Option<f64>> {
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
        if alt.is_empty()
            || alt == b"."
            || alt == b"*"
            || (alt.first() == Some(&b'<') && alt.last() == Some(&b'>'))
            || alt.contains(&b'[')
            || alt.contains(&b']')
        {
            return None;
        }
        Some(alt.to_vec())
    }

    fn parse_required_format_scalar_f32_for_sample(
        format_fields: &[FormatField],
        key: &[u8],
        sample_idx: usize,
    ) -> Result<f32> {
        let value = Self::parse_format_scalar_f32_for_sample(format_fields, key, sample_idx)?;
        let Some(value) = value else {
            return Err(crate::svx_error!(
                "Missing {} FORMAT value in VCF record",
                String::from_utf8_lossy(key)
            ));
        };
        Ok(value)
    }

    fn parse_format_scalar_f32_for_sample(
        format_fields: &[FormatField],
        key: &[u8],
        sample_idx: usize,
    ) -> Result<Option<f32>> {
        let Some(field) = format_fields.iter().find(|field| field.tag == key) else {
            return Ok(None);
        };

        match &field.values {
            FormatFieldValues::Float(values_per_sample) => {
                let Some(sample_values) = values_per_sample.get(sample_idx) else {
                    return Ok(None);
                };
                let Some(value) = sample_values.first() else {
                    return Ok(None);
                };
                if value.is_nan() {
                    return Ok(None);
                }
                Ok(Some(*value))
            }
            FormatFieldValues::Integer(values_per_sample) => {
                let Some(sample_values) = values_per_sample.get(sample_idx) else {
                    return Ok(None);
                };
                let Some(value) = sample_values.first() else {
                    return Ok(None);
                };
                if *value == MISSING_INTEGER || *value == VECTOR_END_INTEGER {
                    return Ok(None);
                }
                Ok(Some(*value as f32))
            }
            FormatFieldValues::String(_) => Err(crate::svx_error!(
                "Error reading {} FORMAT from VCF record as float or integer",
                String::from_utf8_lossy(key)
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::containers::interval_tree::{Interval, IntervalTree};
    use crate::core::variant::test_utils::make_temp_vcf;
    use crate::io::bed_reader::TrId;
    use crate::io::vcf_reader::SampleMapping;
    use rust_htslib::bcf::{self, Read};

    #[test]
    fn test_info_hash_is_hash_of_formatted_vcf_info_field() {
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

        let expected = stable_hash(&bytes[start..end]);
        let got = VariantInternal::info_hash_from_record(&record);
        assert_eq!(got, expected);
    }

    #[test]
    fn cnv_record_parses_cn_and_cnq_format() {
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

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        let format_data = variant
            .vcf
            .as_ref()
            .and_then(VcfWriteData::cnv_format)
            .expect("expected CNV format payload");
        assert_eq!(format_data.cn, 1.0);
        assert_eq!(format_data.cnq, 8.0);
    }

    #[test]
    fn cnv_record_parses_cn_and_cnq_from_supporting_sample_when_first_sample_missing() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t10\tsawfish:CNV:0:0:10\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:CN:CNQ\t./.:.:.\t0/1:2:9
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        let format_data = variant
            .vcf
            .as_ref()
            .and_then(VcfWriteData::cnv_format)
            .expect("expected CNV format payload");
        assert_eq!(variant.support_mask, vec![0b10]);
        assert_eq!(format_data.cn, 2.0);
        assert_eq!(format_data.cnq, 9.0);
    }

    #[test]
    fn cnv_record_parses_svclaim_info() {
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

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.svclaim.as_deref(), Some("D"));
        assert_eq!(variant.svclaims, vec!["D".to_string()]);
    }

    #[test]
    fn non_cnv_del_record_parses_svclaim_info() {
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

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.svclaim.as_deref(), Some("DJ"));
        assert_eq!(variant.svclaims, vec!["DJ".to_string()]);
    }

    #[test]
    fn cnv_record_parses_svclaim_set_info() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SVCLAIM_SET,Number=.,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10\tparent\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;SVLEN=100;SVCLAIM_SET=J,D\tGT:CN:CNQ\t0/1:1:8
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.svclaim, None);
        assert_eq!(variant.svclaims, vec!["D".to_string(), "J".to_string()]);
    }

    #[test]
    fn merged_record_parses_ci_fields_from_info() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t120\tmerged\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-20;SUPP_CALLS=2;END=140;CIPOS=-20,20;CIEND=-10,30\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.support_calls, 2);
        assert!((variant.start_mean - 120.0).abs() < 1e-6);
        assert!((variant.start_variance - 104.123_281_965_847).abs() < 1e-6);
        assert!((variant.svlen_mean + 30.0).abs() < 1e-6);
        assert!((variant.svlen_variance - 208.246_563_931_694).abs() < 1e-6);
    }

    #[test]
    fn merged_record_ignores_unrelated_internal_svx_tags() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"\">
##INFO=<ID=SVX_START_SUM,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t120\tmerged\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-20;SUPP_CALLS=2;END=140;CIPOS=-20,20;CIEND=-10,30;SVX_START_SUM=200\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.support_calls, 2);
        assert!((variant.start_mean - 120.0).abs() < 1e-6);
        assert!((variant.start_variance - 104.123_281_965_847).abs() < 1e-6);
        assert!((variant.svlen_mean + 30.0).abs() < 1e-6);
        assert!((variant.svlen_variance - 208.246_563_931_694).abs() < 1e-6);
    }

    #[test]
    fn merged_record_rejects_partial_ci_fields() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t120\tmerged\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-20;SUPP_CALLS=2;CIPOS=-20,20\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let error = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .expect_err("partial CI fields should fail");
        assert!(
            error.to_string().contains("incomplete CI fields"),
            "unexpected error: {error}"
        );
    }

    #[test]
    fn cnv_record_errors_when_cnq_missing() {
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

        let err = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap_err();
        assert!(
            err.to_string()
                .contains("Missing CNQ FORMAT value in VCF record"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn multisample_record_computes_support_mask_from_all_sample_gts() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT\t0/0\t1/1\t./.
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.support_mask, vec![0b010]);
    }

    #[test]
    fn support_calls_defaults_to_one_when_info_is_missing() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.support_calls, 1);
    }

    #[test]
    fn support_calls_uses_info_value_when_present() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;SUPP_CALLS=3\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.support_calls, 3);
    }

    #[test]
    fn support_calls_rejects_non_positive_values() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;SUPP_CALLS=0\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let err = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("SUPP_CALLS=0 is invalid"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn id_list_defaults_to_normalized_id_when_info_is_missing() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.id_list, vec!["0_var1".to_string()]);
    }

    #[test]
    fn id_list_uses_input_info_values_when_present() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=IDLIST,Number=.,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t10\tmerged_parent\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;IDLIST=a,b\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.id_list, vec!["0_a".to_string(), "0_b".to_string()]);
    }

    #[test]
    fn id_list_remaps_prefixed_sample_indices_when_local_mapping_is_present() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=IDLIST,Number=.,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
chr1\t10\t0_parent\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;IDLIST=0_child_a,1_child_b\tGT\t0/1\t0/0
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();
        let format_cache = VariantInternal::build_header_format_cache(record.header()).unwrap();

        let variant = VariantInternal::from_vcf_record(
            &record,
            1,
            "chr1",
            &args,
            &None,
            &format_cache,
            Some(&[5, 6]),
        )
        .unwrap();

        assert_eq!(
            variant.id_list,
            vec!["5_child_a".to_string(), "6_child_b".to_string()]
        );
    }

    #[test]
    fn multisample_record_preserves_format_order_and_string_values() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PF,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT:GQ:PF\t0/1:60:PASS\t0/0:20:LOW
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        let format_data = variant
            .vcf
            .as_ref()
            .and_then(|vcf| vcf.format_data.as_ref())
            .expect("expected parsed record format payload");
        assert_eq!(
            format_data.format_order,
            vec![b"GT".to_vec(), b"GQ".to_vec(), b"PF".to_vec()]
        );
        assert_eq!(variant.support_mask, vec![0b01]);
        assert_eq!(format_data.sample_gts.len(), 2);

        let pf_values = format_data
            .fields
            .iter()
            .find(|field| field.tag == b"PF")
            .expect("PF tag should be present");
        match &pf_values.values {
            FormatFieldValues::String(values) => {
                assert_eq!(values, &vec![b"PASS".to_vec(), b"LOW".to_vec()]);
            }
            _ => panic!("PF should be stored as string FORMAT values"),
        }
    }

    #[test]
    fn cached_header_format_lookup_matches_uncached_record_parse() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PF,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT:GQ:PF\t0/1:60:PASS\t0/0:20:LOW
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let format_cache = VariantInternal::build_header_format_cache(record.header()).unwrap();

        let with_cache =
            VariantInternal::from_vcf_record(&record, 0, "chr1", &args, &None, &format_cache, None)
                .unwrap();
        let without_cache = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();

        assert_eq!(with_cache.support_mask, without_cache.support_mask);
        let with_cache_data = with_cache
            .vcf
            .as_ref()
            .and_then(|vcf| vcf.format_data.as_ref())
            .expect("expected parsed format data");
        let without_cache_data = without_cache
            .vcf
            .as_ref()
            .and_then(|vcf| vcf.format_data.as_ref())
            .expect("expected parsed format data");
        assert_eq!(
            with_cache_data.format_order,
            without_cache_data.format_order
        );
    }

    #[test]
    fn internal_variant_id_remaps_prefixed_sample_index_when_local_mapping_is_present() {
        let remapped =
            VariantInternal::internal_variant_id(1, b"0_sawfish:0:49906:0:0", Some(&[5, 6, 7]));
        assert_eq!(remapped, "5_sawfish:0:49906:0:0");
    }

    #[test]
    fn from_vcf_record_with_cache_uses_local_mapping_for_prefixed_ids() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
chr1\t10\t0_var1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT\t0/1\t0/0
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();
        let format_cache = VariantInternal::build_header_format_cache(record.header()).unwrap();

        let variant = VariantInternal::from_vcf_record(
            &record,
            1,
            "chr1",
            &args,
            &None,
            &format_cache,
            Some(&[5, 6]),
        )
        .unwrap();

        assert_eq!(variant.id, "5_var1");
    }

    #[test]
    fn remap_support_mask_with_local_map_matches_hash_lookup_path() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3
chr1\t10\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT\t1/1\t0/0\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();

        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 2);
        sample_mapping.index_map.insert((0, 1), 0);
        sample_mapping.index_map.insert((0, 2), 1);

        let mut hash_lookup = variant.clone();
        hash_lookup
            .remap_support_mask_to_output(&sample_mapping, 3)
            .unwrap();

        let mut local_map = variant;
        local_map
            .remap_support_mask_to_output_with_local_map(&[2, 0, 1], 3)
            .unwrap();

        assert_eq!(hash_lookup.support_mask, local_map.support_mask);
    }

    #[test]
    fn tr_annotation_is_limited_to_ins_and_del() {
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
        let format_cache = VariantInternal::build_header_format_cache(ins_record.header()).unwrap();

        let ins_variant = VariantInternal::from_vcf_record(
            &ins_record,
            0,
            "chr1",
            &args,
            &tr_it,
            &format_cache,
            None,
        )
        .unwrap();
        let del_variant = VariantInternal::from_vcf_record(
            &del_record,
            0,
            "chr1",
            &args,
            &tr_it,
            &format_cache,
            None,
        )
        .unwrap();
        let inv_variant = VariantInternal::from_vcf_record(
            &inv_record,
            0,
            "chr1",
            &args,
            &tr_it,
            &format_cache,
            None,
        )
        .unwrap();

        assert!(ins_variant.trid.is_some());
        assert!(del_variant.trid.is_some());
        assert!(inv_variant.trid.is_none());
    }

    #[test]
    fn ins_symbolic_alt_does_not_set_sequence() {
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

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert!(variant.sequence.is_none());
    }

    #[test]
    fn ins_symbolic_subtype_alt_does_not_set_sequence() {
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

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert!(variant.sequence.is_none());
    }

    #[test]
    fn ins_literal_alt_sets_sequence() {
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

        let variant = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(variant.sequence, Some(b"NAAAA".to_vec()));
    }

    #[test]
    fn svlen_type_error_is_propagated_for_non_bnd_records() {
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

        let err = VariantInternal::from_vcf_record(
            &record,
            0,
            "chr1",
            &args,
            &None,
            &VariantInternal::build_header_format_cache(record.header()).unwrap(),
            None,
        )
        .unwrap_err();
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
