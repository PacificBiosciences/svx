use crate::{
    cli::{FULL_VERSION, MergeArgs},
    core::{
        svtype::SvType,
        variant::{FormatFieldValues, VariantInternal, VcfWriteData},
        variant_block::VariantBlockResult,
    },
    io::{
        vcf_reader::{SampleMapping, VcfReaders},
        vcf_writer::VcfWriter,
    },
    utils::util::{
        MISSING_FLOAT, MISSING_INTEGER, Result, VECTOR_END_FLOAT, VECTOR_END_INTEGER, round_to_i64,
        stable_hash, to_info_i32,
    },
};
use rust_htslib::bcf;
use std::{
    collections::{BTreeSet, HashMap},
    env, iter,
};

pub fn create_output_header(
    vcf_readers: &VcfReaders,
    args: &MergeArgs,
) -> Result<(bcf::Header, SampleMapping)> {
    let mut out_header = bcf::Header::new();
    let sample_mapping = vcf_readers.merge_headers(&mut out_header)?;

    out_header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">"#);
    out_header.push_record(
        br#"##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description="Number of input calls merged into this variant">"#,
    );
    out_header.push_record(br#"##INFO=<ID=IDLIST,Number=.,Type=String,Description="Variant IDs of input calls merged to make this call">"#);
    out_header.push_record(
        br#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for merged variants">"#,
    );
    out_header.push_record(
        br#"##INFO=<ID=END,Number=1,Type=Integer,Description="Merged end coordinate">"#,
    );
    out_header.push_record(
        br#"##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for merged variants">"#,
    );
    out_header.push_record(
        br#"##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate contig for breakend pairs">"#,
    );
    out_header.push_record(br#"##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakend orientation for breakend pairs">"#);
    out_header.push_record(br#"##INFO=<ID=TR_CONTAINED,Number=1,Type=String,Description="Tandem repeat ID for merged groups containing at least one TR-annotated variant">"#);
    out_header.push_record(br#"##INFO=<ID=SVCLAIM_SET,Number=.,Type=String,Description="Unique SVCLAIM values observed across merged DEL/DUP/CNV inputs">"#);
    out_header.push_record(br#"##FILTER=<ID=InvBreakpoint,Description="Breakpoint represented as part of an inversion record (same EVENT ID)">"#);
    out_header.push_record(br#"##INFO=<ID=EVENT,Number=1,Type=String,Description="Event identifier for breakend pairs">"#);
    out_header.push_record(br#"##INFO=<ID=EVENTTYPE,Number=1,Type=String,Description="Type of associated event for breakend pairs">"#);
    out_header.push_record(
        br#"##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">"#,
    );

    if !args.no_version {
        add_version_info(&mut out_header);
    }
    Ok((out_header, sample_mapping))
}

fn add_version_info(out_header: &mut bcf::Header) {
    let version_line = format!("##{}Version={}", env!("CARGO_PKG_NAME"), FULL_VERSION);
    out_header.push_record(version_line.as_bytes());

    let command_line = env::args().collect::<Vec<String>>().join(" ");
    let command_line = format!("##{}Command={}", env!("CARGO_PKG_NAME"), command_line);
    out_header.push_record(command_line.as_bytes());
}

fn create_sample_representatives(
    group: &[VariantInternal],
    n_samples: usize,
) -> Result<Vec<Option<usize>>> {
    let mut sample_representatives = vec![None; n_samples];

    for (group_idx, variant) in group.iter().enumerate() {
        for merged_pos in support_mask_indices(&variant.support_mask, n_samples) {
            if sample_representatives[merged_pos].is_none() {
                sample_representatives[merged_pos] = Some(group_idx);
            }
        }
    }

    Ok(sample_representatives)
}

fn support_mask_indices(support_mask: &[u64], n_samples: usize) -> Vec<usize> {
    let mut indices = Vec::new();
    for (word_idx, word) in support_mask.iter().copied().enumerate() {
        let mut remaining = word;
        while remaining != 0 {
            let bit_idx = remaining.trailing_zeros() as usize;
            let sample_idx = word_idx * (u64::BITS as usize) + bit_idx;
            if sample_idx < n_samples {
                indices.push(sample_idx);
            }
            remaining &= remaining - 1;
        }
    }
    indices
}

fn first_supporting_sample_index(support_mask: &[u64], n_samples: usize) -> usize {
    support_mask_indices(support_mask, n_samples)
        .into_iter()
        .next()
        .unwrap_or(n_samples)
}

fn sample_support_count(sample_representatives: &[Option<usize>]) -> usize {
    sample_representatives
        .iter()
        .filter(|sample_variant| sample_variant.is_some())
        .count()
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum FormatFieldType {
    Integer,
    Float,
    String,
}

impl FormatFieldType {
    fn as_str(self) -> &'static str {
        match self {
            Self::Integer => "Integer",
            Self::Float => "Float",
            Self::String => "String",
        }
    }
}

fn output_sample_source(
    sample_mapping: &SampleMapping,
    output_sample_idx: usize,
) -> Result<(usize, usize)> {
    sample_mapping
        .reverse_map
        .get(&output_sample_idx)
        .copied()
        .ok_or_else(|| {
            crate::svx_error!(
                "Missing reverse sample mapping for output sample index {}",
                output_sample_idx
            )
        })
}

fn selected_group_sample_for_output(
    group: &[VariantInternal],
    sample_mapping: &SampleMapping,
    output_sample_idx: usize,
    sample_representatives: &[Option<usize>],
    singleton_group: bool,
) -> Result<Option<(usize, usize)>> {
    let (sample_vcf_id, local_sample_idx) =
        output_sample_source(sample_mapping, output_sample_idx)?;
    if singleton_group {
        return Ok((sample_vcf_id == group[0].vcf_id).then_some((0usize, local_sample_idx)));
    }

    let Some(group_idx) = sample_representatives[output_sample_idx] else {
        return Ok(None);
    };
    if sample_vcf_id != group[group_idx].vcf_id {
        return Err(crate::svx_error!(
            "Output sample {} maps to VCF {}, but supporting variant {} comes from VCF {}",
            output_sample_idx,
            sample_vcf_id,
            group[group_idx].id,
            group[group_idx].vcf_id
        ));
    }

    Ok(Some((group_idx, local_sample_idx)))
}

fn infer_format_field_type<F>(
    group: &[VariantInternal],
    payload_for_variant: &F,
    tag: &[u8],
) -> Result<FormatFieldType>
where
    F: Fn(&VariantInternal) -> Result<&VcfWriteData>,
{
    let mut inferred_type: Option<FormatFieldType> = None;
    for variant in group {
        let payload = payload_for_variant(variant)?;
        if let Some(values) = payload.format_field(tag) {
            let field_type = match values {
                FormatFieldValues::Integer(_) => FormatFieldType::Integer,
                FormatFieldValues::Float(_) => FormatFieldType::Float,
                FormatFieldValues::String(_) => FormatFieldType::String,
            };
            if let Some(expected_type) = inferred_type {
                if expected_type != field_type {
                    return Err(crate::svx_error!(
                        "FORMAT tag {} has inconsistent types across merged variants: expected {}, found {} in variant {}",
                        String::from_utf8_lossy(tag),
                        expected_type.as_str(),
                        field_type.as_str(),
                        variant.id
                    ));
                }
            } else {
                inferred_type = Some(field_type);
            }
        }
    }

    inferred_type.ok_or_else(|| {
        crate::svx_error!(
            "FORMAT tag {} is missing from all merged variants",
            String::from_utf8_lossy(tag)
        )
    })
}

fn merged_format_order<F>(
    group: &[VariantInternal],
    payload_for_variant: &F,
) -> Result<Vec<Vec<u8>>>
where
    F: Fn(&VariantInternal) -> Result<&VcfWriteData>,
{
    let mut merged_tags = BTreeSet::new();
    for variant in group {
        let payload = payload_for_variant(variant)?;
        for tag in payload.format_order() {
            if tag.as_slice() != b"GT" {
                merged_tags.insert(tag);
            }
        }
    }

    let mut format_order = vec![b"GT".to_vec()];
    format_order.extend(merged_tags);
    Ok(prioritize_sv_format_prefix(format_order))
}

fn prioritize_sv_format_prefix(format_order: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut has_gq = false;
    let mut has_pl = false;
    let mut has_ad = false;
    let mut extra_tags = Vec::new();
    for tag in format_order {
        match tag.as_slice() {
            b"GT" => {}
            b"GQ" => has_gq = true,
            b"PL" => has_pl = true,
            b"AD" => has_ad = true,
            _ => {
                if !extra_tags.iter().any(|existing| existing == &tag) {
                    extra_tags.push(tag);
                }
            }
        }
    }

    let mut ordered = Vec::with_capacity(1 + extra_tags.len() + 3);
    ordered.push(b"GT".to_vec());
    if has_gq {
        ordered.push(b"GQ".to_vec());
    }
    if has_pl {
        ordered.push(b"PL".to_vec());
    }
    if has_ad {
        ordered.push(b"AD".to_vec());
    }
    ordered.extend(extra_tags);
    ordered
}

fn flatten_integer_values(values_per_sample: &[Vec<i32>]) -> Vec<i32> {
    let max_len = values_per_sample
        .iter()
        .map(Vec::len)
        .max()
        .unwrap_or(1)
        .max(1);
    let mut flattened = Vec::with_capacity(values_per_sample.len() * max_len);
    for values in values_per_sample {
        if values.is_empty() {
            flattened.push(MISSING_INTEGER);
            flattened.extend(iter::repeat_n(VECTOR_END_INTEGER, max_len - 1));
            continue;
        }
        flattened.extend(values.iter().copied());
        if values.len() < max_len {
            flattened.extend(iter::repeat_n(VECTOR_END_INTEGER, max_len - values.len()));
        }
    }
    flattened
}

fn flatten_float_values(values_per_sample: &[Vec<f32>]) -> Vec<f32> {
    let max_len = values_per_sample
        .iter()
        .map(Vec::len)
        .max()
        .unwrap_or(1)
        .max(1);
    let mut flattened = Vec::with_capacity(values_per_sample.len() * max_len);
    for values in values_per_sample {
        if values.is_empty() {
            flattened.push(MISSING_FLOAT);
            flattened.extend(iter::repeat_n(VECTOR_END_FLOAT, max_len - 1));
            continue;
        }
        flattened.extend(values.iter().copied());
        if values.len() < max_len {
            flattened.extend(iter::repeat_n(VECTOR_END_FLOAT, max_len - values.len()));
        }
    }
    flattened
}

fn write_group_format<F>(
    out_rec: &mut bcf::Record,
    group: &[VariantInternal],
    sample_mapping: &SampleMapping,
    n_samples: usize,
    sample_representatives: &[Option<usize>],
    payload_for_variant: F,
) -> Result<()>
where
    F: Fn(&VariantInternal) -> Result<&VcfWriteData>,
{
    let singleton_group = group.len() == 1;
    let format_order = if singleton_group {
        prioritize_sv_format_prefix(payload_for_variant(&group[0])?.format_order())
    } else {
        merged_format_order(group, &payload_for_variant)?
    };

    let mut selected_group_samples = Vec::with_capacity(n_samples);
    for output_sample_idx in 0..n_samples {
        selected_group_samples.push(selected_group_sample_for_output(
            group,
            sample_mapping,
            output_sample_idx,
            sample_representatives,
            singleton_group,
        )?);
    }

    let mut max_ploidy = 2usize;
    let mut sample_gts = Vec::with_capacity(n_samples);
    for selected_group_sample in &selected_group_samples {
        if let Some((group_idx, local_sample_idx)) = selected_group_sample {
            let payload = payload_for_variant(&group[*group_idx])?;
            let gt = payload
                .sample_gt(*local_sample_idx)
                .map_or_else(Vec::new, |alleles| alleles.to_vec());
            max_ploidy = max_ploidy.max(gt.len());
            sample_gts.push(gt);
        } else {
            sample_gts.push(Vec::new());
        }
    }

    for gt in &mut sample_gts {
        if gt.is_empty() {
            gt.extend(iter::repeat_n(
                bcf::record::GenotypeAllele::UnphasedMissing,
                max_ploidy,
            ));
        } else if gt.len() < max_ploidy {
            gt.extend(iter::repeat_n(
                bcf::record::GenotypeAllele::UnphasedMissing,
                max_ploidy - gt.len(),
            ));
        }
    }
    let flattened_gts = sample_gts
        .iter()
        .flat_map(|genotype| genotype.iter().copied())
        .collect::<Vec<_>>();
    out_rec.push_genotypes(flattened_gts.as_slice())?;

    for tag in format_order {
        if tag.as_slice() == b"GT" {
            continue;
        }

        match infer_format_field_type(group, &payload_for_variant, tag.as_slice())? {
            FormatFieldType::Integer => {
                let mut values_per_sample = Vec::with_capacity(n_samples);
                for selected_group_sample in &selected_group_samples {
                    if let Some((group_idx, local_sample_idx)) = selected_group_sample {
                        let payload = payload_for_variant(&group[*group_idx])?;
                        let values = match payload.format_field(tag.as_slice()) {
                            Some(FormatFieldValues::Integer(per_sample_values)) => {
                                per_sample_values
                                    .get(*local_sample_idx)
                                    .cloned()
                                    .unwrap_or_else(|| vec![MISSING_INTEGER])
                            }
                            Some(FormatFieldValues::Float(_))
                            | Some(FormatFieldValues::String(_))
                            | None => vec![MISSING_INTEGER],
                        };
                        values_per_sample.push(values);
                    } else {
                        values_per_sample.push(vec![MISSING_INTEGER]);
                    }
                }
                let flattened = flatten_integer_values(&values_per_sample);
                out_rec.push_format_integer(tag.as_slice(), flattened.as_slice())?;
            }
            FormatFieldType::Float => {
                let mut values_per_sample = Vec::with_capacity(n_samples);
                for selected_group_sample in &selected_group_samples {
                    if let Some((group_idx, local_sample_idx)) = selected_group_sample {
                        let payload = payload_for_variant(&group[*group_idx])?;
                        let values = match payload.format_field(tag.as_slice()) {
                            Some(FormatFieldValues::Float(per_sample_values)) => per_sample_values
                                .get(*local_sample_idx)
                                .cloned()
                                .unwrap_or_else(|| vec![MISSING_FLOAT]),
                            Some(FormatFieldValues::Integer(_))
                            | Some(FormatFieldValues::String(_))
                            | None => vec![MISSING_FLOAT],
                        };
                        values_per_sample.push(values);
                    } else {
                        values_per_sample.push(vec![MISSING_FLOAT]);
                    }
                }
                let flattened = flatten_float_values(&values_per_sample);
                out_rec.push_format_float(tag.as_slice(), flattened.as_slice())?;
            }
            FormatFieldType::String => {
                let mut values_per_sample = Vec::with_capacity(n_samples);
                for selected_group_sample in &selected_group_samples {
                    if let Some((group_idx, local_sample_idx)) = selected_group_sample {
                        let payload = payload_for_variant(&group[*group_idx])?;
                        let values = match payload.format_field(tag.as_slice()) {
                            Some(FormatFieldValues::String(per_sample_values)) => per_sample_values
                                .get(*local_sample_idx)
                                .cloned()
                                .unwrap_or_else(|| b".".to_vec()),
                            Some(FormatFieldValues::Integer(_))
                            | Some(FormatFieldValues::Float(_))
                            | None => b".".to_vec(),
                        };
                        values_per_sample.push(values);
                    } else {
                        values_per_sample.push(b".".to_vec());
                    }
                }
                out_rec.push_format_string(tag.as_slice(), values_per_sample.as_slice())?;
            }
        }
    }

    Ok(())
}

fn bnd_anchor_seq_from_alt(alt: &[u8]) -> Result<Vec<u8>> {
    let first_bracket_idx = alt
        .iter()
        .position(|&b| b == b'[' || b == b']')
        .ok_or_else(|| {
            crate::svx_error!(
                "BND ALT is missing '[' or ']': {:?}",
                String::from_utf8_lossy(alt)
            )
        })?;
    let bracket = alt[first_bracket_idx];
    let second_bracket_idx = alt
        .iter()
        .skip(first_bracket_idx + 1)
        .position(|&b| b == bracket)
        .map(|rel| rel + first_bracket_idx + 1)
        .ok_or_else(|| {
            crate::svx_error!(
                "BND ALT is missing closing bracket '{}': {:?}",
                bracket as char,
                String::from_utf8_lossy(alt)
            )
        })?;
    let anchor_seq = if first_bracket_idx == 0 {
        &alt[second_bracket_idx + 1..]
    } else {
        &alt[..first_bracket_idx]
    };
    if anchor_seq.is_empty() {
        return Err(crate::svx_error!(
            "BND ALT anchor sequence is empty: {:?}",
            String::from_utf8_lossy(alt)
        ));
    }
    Ok(anchor_seq.to_vec())
}

fn bnd_alt(
    mate_contig: &str,
    mate_pos_1based: i64,
    strands: [u8; 2],
    anchor_seq: &[u8],
) -> Result<Vec<u8>> {
    let mate = format!("{mate_contig}:{mate_pos_1based}");
    let mate = mate.as_bytes();
    let mut alt = Vec::with_capacity(anchor_seq.len() + mate.len() + 2);

    match &strands {
        b"+-" => {
            alt.push(b'[');
            alt.extend_from_slice(mate);
            alt.push(b'[');
            alt.extend_from_slice(anchor_seq);
        }
        b"--" => {
            alt.push(b']');
            alt.extend_from_slice(mate);
            alt.push(b']');
            alt.extend_from_slice(anchor_seq);
        }
        b"++" => {
            alt.extend_from_slice(anchor_seq);
            alt.push(b'[');
            alt.extend_from_slice(mate);
            alt.push(b'[');
        }
        b"-+" => {
            alt.extend_from_slice(anchor_seq);
            alt.push(b']');
            alt.extend_from_slice(mate);
            alt.push(b']');
        }
        _ => {
            return Err(crate::svx_error!(
                "Invalid BND STRANDS {:?}; expected one of ++, +-, -+, --",
                std::str::from_utf8(&strands).unwrap_or("<non-utf8>")
            ));
        }
    };

    Ok(alt)
}

fn deterministic_bnd_anchor_allele<F>(
    group: &[VariantInternal],
    side_label: &str,
    candidate_for_variant: F,
) -> Result<(u8, Vec<u8>)>
where
    F: Fn(&VariantInternal) -> Result<(u8, Vec<u8>)>,
{
    let mut selected_anchor_allele: Option<(u8, Vec<u8>)> = None;

    for variant in group {
        let candidate_anchor_allele = candidate_for_variant(variant)?;
        let replace_selected = selected_anchor_allele.as_ref().is_none_or(
            |(selected_ref_base, selected_anchor_seq)| {
                candidate_anchor_allele
                    .0
                    .cmp(selected_ref_base)
                    .then_with(|| {
                        candidate_anchor_allele
                            .1
                            .as_slice()
                            .cmp(selected_anchor_seq.as_slice())
                    })
                    == std::cmp::Ordering::Less
            },
        );
        if replace_selected {
            selected_anchor_allele = Some(candidate_anchor_allele);
        }
    }

    selected_anchor_allele.ok_or_else(|| {
        crate::svx_error!(
            "Cannot select deterministic BND {side_label}-side anchor allele from empty group"
        )
    })
}

#[derive(Clone, Copy)]
struct GroupAggregates {
    support_calls: i64,
    start_sum: f64,
    end_sum: f64,
    svlen_sum: f64,
}

fn checked_add_weighted(
    total: &mut f64,
    weight: f64,
    value: f64,
    label: &str,
    variant_id: &str,
) -> Result<()> {
    if !value.is_finite() {
        return Err(crate::svx_error!(
            "{label} is not finite for variant {variant_id}"
        ));
    }

    let term = weight * value;
    if !term.is_finite() {
        return Err(crate::svx_error!(
            "{label} overflows while accumulating weighted aggregate for variant {variant_id}"
        ));
    }
    *total += term;
    if !total.is_finite() {
        return Err(crate::svx_error!(
            "{label} aggregate became non-finite while processing variant {variant_id}"
        ));
    }
    Ok(())
}

fn accumulate_group_aggregates(group: &[VariantInternal]) -> Result<GroupAggregates> {
    let mut support_calls = 0i64;
    let mut start_sum = 0.0f64;
    let mut end_sum = 0.0f64;
    let mut svlen_sum = 0.0f64;

    for variant in group {
        support_calls = support_calls
            .checked_add(variant.support_calls)
            .ok_or_else(|| {
                crate::svx_error!(
                    "SUPP_CALLS overflow while summing group support counts for output record"
                )
            })?;
        if variant.support_calls <= 0 {
            return Err(crate::svx_error!(
                "SUPP_CALLS={} is invalid in variant {}; expected a positive integer",
                variant.support_calls,
                variant.id
            ));
        }
        let weight = variant.support_calls as f64;
        let end_mean = variant.end + 1.0;

        checked_add_weighted(
            &mut start_sum,
            weight,
            variant.start_mean,
            "START_AVG",
            &variant.id,
        )?;
        checked_add_weighted(&mut end_sum, weight, end_mean, "AVG_END", &variant.id)?;
        checked_add_weighted(
            &mut svlen_sum,
            weight,
            variant.svlen_mean,
            "SVLEN_AVG",
            &variant.id,
        )?;
    }

    Ok(GroupAggregates {
        support_calls,
        start_sum,
        end_sum,
        svlen_sum,
    })
}

fn aggregate_mean(sum: f64, count: i64, label: &str) -> Result<f64> {
    if count <= 0 {
        return Err(crate::svx_error!(
            "Cannot compute {label}: aggregate count must be positive"
        ));
    }

    let count_f64 = count as f64;
    let mean = sum / count_f64;
    if !mean.is_finite() {
        return Err(crate::svx_error!(
            "Cannot compute {label}: aggregate mean is non-finite"
        ));
    }

    Ok(mean)
}

fn aggregate_variance<FMean, FVariance>(
    group: &[VariantInternal],
    count: i64,
    aggregate_mean: f64,
    mean_fn: FMean,
    variance_fn: FVariance,
    label: &str,
) -> Result<f64>
where
    FMean: Fn(&VariantInternal) -> f64,
    FVariance: Fn(&VariantInternal) -> f64,
{
    if count <= 0 {
        return Err(crate::svx_error!(
            "Cannot compute {label}: aggregate count must be positive"
        ));
    }

    let mut weighted_variance_sum = 0.0f64;
    for variant in group {
        if variant.support_calls <= 0 {
            return Err(crate::svx_error!(
                "Cannot compute {label}: SUPP_CALLS={} is invalid for variant {}",
                variant.support_calls,
                variant.id
            ));
        }

        let component_mean = mean_fn(variant);
        let component_variance = variance_fn(variant);
        if component_variance < 0.0 {
            return Err(crate::svx_error!(
                "Cannot compute {label}: component variance is negative ({component_variance}) for variant {}",
                variant.id
            ));
        }

        let delta = component_mean - aggregate_mean;
        let component_total_variance = component_variance + (delta * delta);

        checked_add_weighted(
            &mut weighted_variance_sum,
            variant.support_calls as f64,
            component_total_variance,
            label,
            &variant.id,
        )?;
    }

    Ok(weighted_variance_sum / count as f64)
}

pub const CI_Z_95: f64 = 1.96;
fn ci_offsets(mean: f64, variance: f64, anchor_1based: i64, label: &str) -> Result<[i32; 2]> {
    if !variance.is_finite() {
        return Err(crate::svx_error!(
            "Cannot compute {label}: variance is non-finite"
        ));
    }
    if variance < 0.0 {
        return Err(crate::svx_error!(
            "Cannot compute {label}: variance is negative ({variance})"
        ));
    }

    let half_width = CI_Z_95 * variance.sqrt();
    if !half_width.is_finite() {
        return Err(crate::svx_error!(
            "Cannot compute {label}: CI half-width is non-finite"
        ));
    }

    let lower = round_to_i64(mean - half_width);
    let upper = round_to_i64(mean + half_width);
    let (lower, upper) = if lower <= upper {
        (lower, upper)
    } else {
        (upper, lower)
    };

    let lower_offset = to_info_i32(lower - anchor_1based, label)?;
    let upper_offset = to_info_i32(upper - anchor_1based, label)?;
    Ok([lower_offset, upper_offset])
}

fn push_group_support_info(
    out_rec: &mut bcf::Record,
    group: &[VariantInternal],
    sample_representatives: &[Option<usize>],
) -> Result<()> {
    let supp = sample_support_count(sample_representatives) as i32;
    let supp_calls = group.iter().map(|variant| variant.support_calls).try_fold(
        0i64,
        |acc, support_calls| {
            acc.checked_add(support_calls).ok_or_else(|| {
                crate::svx_error!(
                    "SUPP_CALLS overflow while summing group support counts for output record"
                )
            })
        },
    )?;
    let supp_calls = to_info_i32(supp_calls, "SUPP_CALLS")?;
    let mut source_ids = flattened_group_source_ids(group);
    source_ids.sort_unstable();
    let id_list: Vec<&[u8]> = source_ids
        .iter()
        .map(|source_id| source_id.as_bytes())
        .collect();

    out_rec.push_info_integer(b"SUPP", &[supp])?;
    out_rec.push_info_integer(b"SUPP_CALLS", &[supp_calls])?;
    out_rec.push_info_string(b"IDLIST", &id_list)?;

    Ok(())
}

fn flattened_group_source_ids(group: &[VariantInternal]) -> Vec<&str> {
    let mut source_ids: Vec<&str> = Vec::new();
    for variant in group {
        if variant.id_list.is_empty() {
            source_ids.push(variant.id.as_str());
            continue;
        }
        source_ids.extend(variant.id_list.iter().map(String::as_str));
    }
    source_ids
}

fn merged_svlen_i64(output_svtype: SvType, mean_svlen: f64) -> i64 {
    let abs_len = round_to_i64(mean_svlen.abs());
    if output_svtype == SvType::DELETION {
        -abs_len
    } else {
        abs_len
    }
}

fn push_svclaim_info(
    out_rec: &mut bcf::Record,
    group: &[VariantInternal],
    block_svtype: SvType,
    output_svtype: SvType,
) -> Result<()> {
    if !matches!(
        block_svtype,
        SvType::CNV | SvType::DELETION | SvType::DUPLICATION
    ) {
        return Ok(());
    }

    let mut unique_claims: BTreeSet<&str> = BTreeSet::new();
    let mut n_claimed = 0usize;
    for variant in group {
        if !variant.svclaims.is_empty() {
            n_claimed += 1;
            for claim in &variant.svclaims {
                unique_claims.insert(claim);
            }
        } else if let Some(claim) = variant.svclaim.as_deref() {
            n_claimed += 1;
            unique_claims.insert(claim);
        }
    }

    if unique_claims.is_empty() {
        return Ok(());
    }

    if output_svtype == SvType::CNV {
        let claim_set: Vec<&[u8]> = unique_claims.iter().map(|claim| claim.as_bytes()).collect();
        out_rec.push_info_string(b"SVCLAIM_SET", &claim_set)?;
        return Ok(());
    }

    if unique_claims.len() == 1 && n_claimed == group.len() {
        let claim = unique_claims
            .first()
            .expect("checked non-empty unique claims set");
        out_rec.push_info_string(b"SVCLAIM", &[claim.as_bytes()])?;
    } else {
        let claim_set: Vec<&[u8]> = unique_claims.iter().map(|claim| claim.as_bytes()).collect();
        out_rec.push_info_string(b"SVCLAIM_SET", &claim_set)?;
    }

    Ok(())
}

fn push_non_bnd_group_info(
    out_rec: &mut bcf::Record,
    group: &[VariantInternal],
    sample_representatives: &[Option<usize>],
    block_svtype: SvType,
    output_svtype: SvType,
    representative_pos0: i64,
) -> Result<()> {
    let is_interval_like = block_svtype == SvType::CNV
        || matches!(
            output_svtype,
            SvType::DELETION | SvType::DUPLICATION | SvType::INVERSION | SvType::CNV
        );

    let aggregates = accumulate_group_aggregates(group)?;
    let mean_start = aggregate_mean(aggregates.start_sum, aggregates.support_calls, "AVG_START")?;
    let start_variance = aggregate_variance(
        group,
        aggregates.support_calls,
        mean_start,
        |variant| variant.start_mean,
        |variant| variant.start_variance,
        "START_VARIANCE",
    )?;
    let mean_svlen = aggregate_mean(aggregates.svlen_sum, aggregates.support_calls, "SVLEN_AVG")?;
    let svlen_variance = aggregate_variance(
        group,
        aggregates.support_calls,
        mean_svlen,
        |variant| variant.svlen_mean,
        |variant| variant.svlen_variance,
        "SVLEN_VARIANCE",
    )?;

    out_rec.set_pos(representative_pos0);
    let pos_1based = representative_pos0 + 1;

    out_rec.push_info_string(b"SVTYPE", &[output_svtype.to_string().as_bytes()])?;
    let merged_svlen = merged_svlen_i64(output_svtype, mean_svlen);
    let svlen_i32 = to_info_i32(merged_svlen, "SVLEN")?;
    out_rec.push_info_integer(b"SVLEN", &[svlen_i32])?;

    let cipos = ci_offsets(mean_start, start_variance, pos_1based, "CIPOS")?;
    out_rec.push_info_integer(b"CIPOS", &cipos)?;

    if is_interval_like {
        let mean_end = mean_start + mean_svlen.abs();
        let end_variance = (start_variance + svlen_variance).max(0.0);
        let end_1based = pos_1based + merged_svlen.abs();
        let end_i32 = to_info_i32(end_1based, "END")?;
        out_rec.push_info_integer(b"END", &[end_i32])?;
        let ciend = ci_offsets(mean_end, end_variance, end_1based, "CIEND")?;
        out_rec.push_info_integer(b"CIEND", &ciend)?;
    }

    push_group_support_info(out_rec, group, sample_representatives)?;
    push_svclaim_info(out_rec, group, block_svtype, output_svtype)?;

    if let Some(tr_id) = group
        .iter()
        .find_map(|variant| variant.trid.as_ref().map(|tr| tr.id.as_str()))
    {
        out_rec.push_info_string(b"TR_CONTAINED", &[tr_id.as_bytes()])?;
    }

    Ok(())
}

fn cnv_alt_allele_for_output_svtype(output_svtype: SvType) -> Result<&'static [u8]> {
    match output_svtype {
        SvType::DELETION => Ok(b"<DEL>"),
        SvType::DUPLICATION => Ok(b"<DUP>"),
        SvType::CNV => Ok(b"<CNV>"),
        _ => Err(crate::svx_error!(
            "Invalid CNV output type {}; expected DEL, DUP, or CNV",
            output_svtype
        )),
    }
}

struct BndRecordInfo<'a> {
    event_id: &'a str,
    mate_id: &'a str,
    mate_contig: &'a str,
    strands: [u8; 2],
    is_inversion_event: bool,
}

fn push_bnd_group_info(
    out_rec: &mut bcf::Record,
    group: &[VariantInternal],
    sample_representatives: &[Option<usize>],
    info: BndRecordInfo<'_>,
    cipos: [i32; 2],
) -> Result<()> {
    out_rec.push_info_string(b"SVTYPE", &[b"BND"])?;
    out_rec.push_info_string(b"EVENT", &[info.event_id.as_bytes()])?;
    if info.is_inversion_event {
        out_rec.push_info_string(b"EVENTTYPE", &[b"INV"])?;
        out_rec.push_filter("InvBreakpoint".as_bytes())?;
    }
    out_rec.push_info_string(b"MATEID", &[info.mate_id.as_bytes()])?;
    out_rec.push_info_string(b"CHR2", &[info.mate_contig.as_bytes()])?;
    out_rec.push_info_string(b"STRANDS", &[&info.strands[..]])?;
    out_rec.push_info_integer(b"CIPOS", &cipos)?;
    push_group_support_info(out_rec, group, sample_representatives)?;
    Ok(())
}

fn bnd_inversion_event_consensus<'a>(group: &'a [VariantInternal]) -> Option<&'a str> {
    let mut consensus_event: Option<&'a str> = None;
    for variant in group {
        let inv_event = variant
            .bnd_event
            .as_ref()
            .and_then(|event| event.inv_event_id.as_deref());
        let inv_event = inv_event?;

        if let Some(consensus_event_id) = consensus_event {
            if consensus_event_id != inv_event {
                return None;
            }
        } else {
            consensus_event = Some(inv_event);
        }
    }
    consensus_event
}

fn build_contig_rid_lookup(header: &bcf::header::HeaderView) -> Result<HashMap<String, u32>> {
    let mut lookup = HashMap::new();
    for rid in 0..header.contig_count() {
        let contig_name = std::str::from_utf8(header.rid2name(rid)?).map_err(|error| {
            crate::svx_error!("Output header contig name is not UTF-8 for RID {rid}: {error}")
        })?;
        lookup.insert(contig_name.to_string(), rid);
    }
    Ok(lookup)
}

fn resolve_contig_rid(contig: &str, contig_rid_lookup: &HashMap<String, u32>) -> Result<u32> {
    contig_rid_lookup
        .get(contig)
        .copied()
        .ok_or_else(|| crate::svx_error!("Output header missing contig {contig}"))
}

fn write_bnd_group_with_sink<F>(
    group: &[VariantInternal],
    sample_representatives: &[Option<usize>],
    sample_mapping: &SampleMapping,
    n_samples: usize,
    contig_rid_lookup: &HashMap<String, u32>,
    out_rec: &mut bcf::Record,
    write: &mut F,
) -> Result<()>
where
    F: FnMut(&bcf::Record) -> Result<()>,
{
    if group.is_empty() {
        return Ok(());
    }

    let representative_event = group[0].bnd_event.as_ref().ok_or_else(|| {
        crate::svx_error!(
            "BND group representative {} is missing BND event payload",
            group[0].id
        )
    })?;

    let (a_contig, b_contig) = (
        &representative_event.a_contig,
        &representative_event.b_contig,
    );
    let representative_a_strands = representative_event.a_strands;
    let representative_b_strands = representative_event.b_strands;
    for v in group {
        let e = v
            .bnd_event
            .as_ref()
            .ok_or_else(|| crate::svx_error!("BND event {} is missing BND event payload", v.id))?;
        if &e.a_contig != a_contig || &e.b_contig != b_contig {
            return Err(crate::svx_error!(
                "BND group contigs mismatch: expected {a_contig}/{b_contig}, got {}/{}",
                e.a_contig,
                e.b_contig
            ));
        }
        if e.a_strands != representative_a_strands || e.b_strands != representative_b_strands {
            let representative_a =
                std::str::from_utf8(&representative_a_strands).unwrap_or("<non-utf8>");
            let representative_b =
                std::str::from_utf8(&representative_b_strands).unwrap_or("<non-utf8>");
            let e_a = std::str::from_utf8(&e.a_strands).unwrap_or("<non-utf8>");
            let e_b = std::str::from_utf8(&e.b_strands).unwrap_or("<non-utf8>");
            return Err(crate::svx_error!(
                "BND group strands mismatch for contigs {a_contig}/{b_contig}: expected {representative_a}/{representative_b}, got {e_a}/{e_b} (event {})",
                v.id
            ));
        }
    }

    let aggregates = accumulate_group_aggregates(group)?;
    let mean_a_1based =
        aggregate_mean(aggregates.start_sum, aggregates.support_calls, "AVG_START")?;
    let var_a = aggregate_variance(
        group,
        aggregates.support_calls,
        mean_a_1based,
        |variant| variant.start_mean,
        |variant| variant.start_variance,
        "CIPOS",
    )?;
    let mean_b_1based = aggregate_mean(aggregates.end_sum, aggregates.support_calls, "AVG_END")?;
    let var_b = aggregate_variance(
        group,
        aggregates.support_calls,
        mean_b_1based,
        |variant| variant.end + 1.0,
        |variant| variant.svlen_variance,
        "CIEND",
    )?;
    let a_pos_1based = round_to_i64(mean_a_1based);
    let b_pos_1based = round_to_i64(mean_b_1based);
    let a_pos0 = a_pos_1based - 1;
    let b_pos0 = b_pos_1based - 1;
    let a_cipos = ci_offsets(mean_a_1based, var_a, a_pos_1based, "CIPOS")?;
    let b_cipos = ci_offsets(mean_b_1based, var_b, b_pos_1based, "CIPOS")?;

    let inversion_event_id = bnd_inversion_event_consensus(group);
    let event_id = if let Some(inv_event_id) = inversion_event_id {
        inv_event_id.to_string()
    } else {
        // Deterministic event ID across one-shot and incremental paths:
        // hash canonical flattened source-call IDs (sorted) rather than intermediate merged IDs.
        let mut source_ids = flattened_group_source_ids(group);
        source_ids.sort_unstable();
        let id_concat = source_ids.join("|");
        let h = stable_hash(id_concat.as_bytes());
        format!(
            "svx_tra_{a_contig}_{}_{}_{}_{}",
            a_pos0 + 1,
            b_contig,
            b_pos0 + 1,
            h
        )
    };
    let is_inversion_event = inversion_event_id.is_some();

    let rec_a_id = format!("{event_id}_A");
    let rec_b_id = format!("{event_id}_B");

    let rid_a = resolve_contig_rid(a_contig, contig_rid_lookup)?;
    let rid_b = resolve_contig_rid(b_contig, contig_rid_lookup)?;

    let (a_ref_base, a_anchor_seq) = deterministic_bnd_anchor_allele(group, "A", |variant| {
        let event = variant.bnd_event.as_ref().ok_or_else(|| {
            crate::svx_error!("BND event {} is missing BND event payload", variant.id)
        })?;
        let template_alt = event.a_vcf.alleles.get(1).ok_or_else(|| {
            crate::svx_error!(
                "BND event {} is missing ALT allele for breakend A",
                variant.id
            )
        })?;
        Ok((event.a_ref_base, bnd_anchor_seq_from_alt(template_alt)?))
    })?;
    let (b_ref_base, b_anchor_seq) = deterministic_bnd_anchor_allele(group, "B", |variant| {
        let event = variant.bnd_event.as_ref().ok_or_else(|| {
            crate::svx_error!("BND event {} is missing BND event payload", variant.id)
        })?;
        let template_alt = event.b_vcf.alleles.get(1).ok_or_else(|| {
            crate::svx_error!(
                "BND event {} is missing ALT allele for breakend B",
                variant.id
            )
        })?;
        Ok((event.b_ref_base, bnd_anchor_seq_from_alt(template_alt)?))
    })?;

    // Record A
    {
        out_rec.clear();
        out_rec.set_rid(Some(rid_a));
        out_rec.set_pos(a_pos0);
        out_rec.set_id(rec_a_id.as_bytes())?;

        let alt = bnd_alt(
            b_contig,
            b_pos0 + 1,
            representative_event.a_strands,
            a_anchor_seq.as_slice(),
        )?;
        out_rec.set_alleles(&[&[a_ref_base], &alt])?;

        push_bnd_group_info(
            out_rec,
            group,
            sample_representatives,
            BndRecordInfo {
                event_id: &event_id,
                mate_id: &rec_b_id,
                mate_contig: b_contig,
                strands: representative_event.a_strands,
                is_inversion_event,
            },
            a_cipos,
        )?;
        write_group_format(
            out_rec,
            group,
            sample_mapping,
            n_samples,
            sample_representatives,
            |variant| {
                variant
                    .bnd_event
                    .as_ref()
                    .map(|event| &event.a_vcf)
                    .ok_or_else(|| {
                        crate::svx_error!(
                            "BND event {} is missing breakend A VCF write data",
                            variant.id
                        )
                    })
            },
        )?;
        write(out_rec)?;
    }

    // Record B
    {
        out_rec.clear();
        out_rec.set_rid(Some(rid_b));
        out_rec.set_pos(b_pos0);
        out_rec.set_id(rec_b_id.as_bytes())?;

        let alt = bnd_alt(
            a_contig,
            a_pos0 + 1,
            representative_event.b_strands,
            b_anchor_seq.as_slice(),
        )?;
        out_rec.set_alleles(&[&[b_ref_base], &alt])?;

        push_bnd_group_info(
            out_rec,
            group,
            sample_representatives,
            BndRecordInfo {
                event_id: &event_id,
                mate_id: &rec_a_id,
                mate_contig: a_contig,
                strands: representative_event.b_strands,
                is_inversion_event,
            },
            b_cipos,
        )?;
        write_group_format(
            out_rec,
            group,
            sample_mapping,
            n_samples,
            sample_representatives,
            |variant| {
                variant
                    .bnd_event
                    .as_ref()
                    .map(|event| &event.b_vcf)
                    .ok_or_else(|| {
                        crate::svx_error!(
                            "BND event {} is missing breakend B VCF write data",
                            variant.id
                        )
                    })
            },
        )?;
        write(out_rec)?;
    }

    Ok(())
}

pub fn write_variants(
    variant_block_result: VariantBlockResult,
    sample_mapping: &SampleMapping,
    n_samples: usize,
    min_supp: usize,
    keep_monomorphic: bool,
    writer: &mut VcfWriter,
) -> Result<()> {
    let contig_rid_lookup = build_contig_rid_lookup(writer.writer.header())?;
    let inner_writer = &mut writer.writer;
    let out_rec = &mut writer.dummy_record;
    let mut write = |record: &bcf::Record| -> Result<()> {
        inner_writer.write(record)?;
        Ok(())
    };
    write_variants_with_sink(
        variant_block_result,
        sample_mapping,
        n_samples,
        min_supp,
        keep_monomorphic,
        &contig_rid_lookup,
        out_rec,
        &mut write,
    )
}

pub fn collect_variants(
    variant_block_result: VariantBlockResult,
    sample_mapping: &SampleMapping,
    n_samples: usize,
    min_supp: usize,
    keep_monomorphic: bool,
    header: &bcf::header::HeaderView,
    out_rec: &mut bcf::Record,
) -> Result<Vec<bcf::Record>> {
    let contig_rid_lookup = build_contig_rid_lookup(header)?;
    let mut records = Vec::new();
    let mut write = |record: &bcf::Record| -> Result<()> {
        records.push(record.clone());
        Ok(())
    };
    write_variants_with_sink(
        variant_block_result,
        sample_mapping,
        n_samples,
        min_supp,
        keep_monomorphic,
        &contig_rid_lookup,
        out_rec,
        &mut write,
    )?;
    Ok(records)
}

#[expect(clippy::too_many_arguments)]
fn write_variants_with_sink<F>(
    variant_block_result: VariantBlockResult,
    sample_mapping: &SampleMapping,
    n_samples: usize,
    min_supp: usize,
    keep_monomorphic: bool,
    contig_rid_lookup: &HashMap<String, u32>,
    out_rec: &mut bcf::Record,
    write: &mut F,
) -> Result<()>
where
    F: FnMut(&bcf::Record) -> Result<()>,
{
    log::debug!(
        "Write: Writing {} block for contig {}",
        variant_block_result.variant_type,
        variant_block_result.contig
    );

    let block_rid = if variant_block_result.variant_type == SvType::BND {
        None
    } else {
        Some(resolve_contig_rid(
            variant_block_result.contig.as_str(),
            contig_rid_lookup,
        )?)
    };

    for mut group in variant_block_result.groups {
        if group.is_empty() {
            continue;
        }

        group.sort_by(|a_variant, b_variant| {
            first_supporting_sample_index(&a_variant.support_mask, n_samples)
                .cmp(&first_supporting_sample_index(
                    &b_variant.support_mask,
                    n_samples,
                ))
                .then_with(|| a_variant.id.cmp(&b_variant.id))
        });
        // Representative variant index per sample in write order.
        let sample_representatives = create_sample_representatives(&group, n_samples)?;
        let supp = sample_support_count(&sample_representatives);

        if !(supp >= min_supp || (keep_monomorphic && supp == 0)) {
            continue;
        }

        let v0 = &group[0];
        if v0.svtype == SvType::BND {
            write_bnd_group_with_sink(
                &group,
                &sample_representatives,
                sample_mapping,
                n_samples,
                contig_rid_lookup,
                out_rec,
                write,
            )?;
            continue;
        }

        let v0_vcf = v0
            .vcf
            .as_ref()
            .ok_or_else(|| crate::svx_error!("Variant {} is missing VCF write data", v0.id))?;
        out_rec.clear();

        out_rec.set_rid(block_rid);
        out_rec.set_id(v0.id.as_bytes())?;
        let output_svtype = if variant_block_result.variant_type == SvType::CNV {
            SvType::CNV
        } else {
            v0.svtype
        };
        if variant_block_result.variant_type == SvType::CNV {
            let ref_allele = v0_vcf.alleles.first().ok_or_else(|| {
                crate::svx_error!(
                    "CNV variant {} is missing REF allele in VCF write payload",
                    v0.id
                )
            })?;
            let alt_allele = cnv_alt_allele_for_output_svtype(output_svtype)?;
            out_rec.set_alleles(&[ref_allele.as_slice(), alt_allele])?;
        } else {
            let alleles: Vec<&[u8]> = v0_vcf.alleles.iter().map(|v| v.as_slice()).collect();
            out_rec.set_alleles(&alleles)?;
        }
        push_non_bnd_group_info(
            out_rec,
            &group,
            &sample_representatives,
            variant_block_result.variant_type,
            output_svtype,
            v0_vcf.pos,
        )?;

        write_group_format(
            out_rec,
            &group,
            sample_mapping,
            n_samples,
            &sample_representatives,
            |variant| {
                variant.vcf.as_ref().ok_or_else(|| {
                    crate::svx_error!("Variant {} is missing VCF write data", variant.id)
                })
            },
        )?;

        write(out_rec)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{MergeArgs, MergeArgsInner, MergeSvType};
    use crate::core::svtype::SvType;
    use crate::core::variant::test_utils;
    use crate::core::variant::{
        BndEventData, CnvSampleFormatData, SampleFormatData, SvSampleFormatData, VcfWriteData,
    };
    use crate::io::bed_reader::TrId;
    use crate::io::vcf_writer::OutputType;
    use crate::{DEFAULT_SORT_MAX_MEM, DEFAULT_SORT_MAX_OPEN_FILES, DEFAULT_SORT_MERGE_FAN_IN};
    use rust_htslib::bcf::Read;
    use rust_htslib::bcf::record::GenotypeAllele;
    use std::{
        fs,
        sync::atomic::{AtomicU64, Ordering},
        time::SystemTime,
    };

    static TEMP_FILE_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn make_temp_vcf(contents: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let counter = TEMP_FILE_COUNTER.fetch_add(1, Ordering::Relaxed);
        path.push(format!(
            "svx_test_merge_writer_bnd_in_{nanos}_{counter}.vcf"
        ));
        fs::write(&path, contents).unwrap();
        path
    }

    fn make_temp_out(ext: &str) -> String {
        let mut path = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let counter = TEMP_FILE_COUNTER.fetch_add(1, Ordering::Relaxed);
        path.push(format!(
            "svx_test_merge_writer_bnd_out_{nanos}_{counter}.{ext}"
        ));
        path.to_string_lossy().to_string()
    }

    fn parse_variant_record(
        record: &rust_htslib::bcf::Record,
        vcf_id: usize,
        contig: &str,
        args: &MergeArgsInner,
    ) -> VariantInternal {
        let format_cache = VariantInternal::build_header_format_cache(record.header()).unwrap();
        VariantInternal::from_vcf_record(record, vcf_id, contig, args, &None, &format_cache, None)
            .unwrap()
    }

    fn sample_mapping_two_samples() -> SampleMapping {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.index_map.insert((1, 0), 1);
        sample_mapping.reverse_map.insert(0, (0, 0));
        sample_mapping.reverse_map.insert(1, (1, 0));
        sample_mapping
    }

    fn sample_mapping_one_sample() -> SampleMapping {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.reverse_map.insert(0, (0, 0));
        sample_mapping
    }

    fn merge_args_for_header_test(no_version: bool) -> MergeArgs {
        MergeArgs {
            vcfs: None,
            vcf_list: None,
            config: None,
            output: None,
            tr_bed_path: None,
            num_threads: 1,
            num_io_threads: 1,
            blob_queue_capacity: None,
            result_queue_capacity: None,
            output_type: None,
            print_header: false,
            sort_output: false,
            sort_max_mem: DEFAULT_SORT_MAX_MEM,
            sort_tmp_dir: None,
            sort_max_open_files: DEFAULT_SORT_MAX_OPEN_FILES,
            sort_merge_fan_in: DEFAULT_SORT_MERGE_FAN_IN,
            progress: false,
            no_progress: false,
            force_single: false,
            no_version,
            contigs: None,
            target_positions: None,
            svtypes: vec![MergeSvType::All],
            merge_args: MergeArgsInner::default(),
        }
    }

    fn header_text(header: &rust_htslib::bcf::Header) -> String {
        let out_path = make_temp_out("vcf");
        let writer =
            rust_htslib::bcf::Writer::from_path(&out_path, header, true, bcf::Format::Vcf).unwrap();
        drop(writer);
        fs::read_to_string(out_path).unwrap()
    }

    #[test]
    fn create_output_header_omits_supp_vec() {
        let vcf_readers = VcfReaders {
            readers: Vec::new(),
            n: 0,
            io_tpool: None,
        };
        let args = merge_args_for_header_test(true);

        let (out_header, _) = create_output_header(&vcf_readers, &args).unwrap();
        let header = header_text(&out_header);

        assert!(
            !header.contains("##INFO=<ID=SUPP_VEC,"),
            "expected SUPP_VEC INFO header to be omitted, got header:\n{header}"
        );
    }

    #[test]
    fn create_output_header_declares_end_chr2_strands_and_svclaim_set() {
        let vcf_readers = VcfReaders {
            readers: Vec::new(),
            n: 0,
            io_tpool: None,
        };
        let args = merge_args_for_header_test(true);

        let (out_header, _) = create_output_header(&vcf_readers, &args).unwrap();
        let header = header_text(&out_header);

        assert!(
            header.contains(
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Merged end coordinate\">"
            ),
            "expected END INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Mate contig for breakend pairs\">"),
            "expected CHR2 INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakend orientation for breakend pairs\">"),
            "expected STRANDS INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains(
                "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for merged variants\">"
            ),
            "expected CIPOS INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains(
                "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for merged variants\">"
            ),
            "expected CIEND INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains("##INFO=<ID=SVCLAIM_SET,Number=.,Type=String,Description=\"Unique SVCLAIM values observed across merged DEL/DUP/CNV inputs\">"),
            "expected SVCLAIM_SET INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains(
                "##INFO=<ID=EVENTTYPE,Number=1,Type=String,Description=\"Type of associated event for breakend pairs\">"
            ),
            "expected EVENTTYPE INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains(
                "##FILTER=<ID=InvBreakpoint,Description=\"Breakpoint represented as part of an inversion record (same EVENT ID)\">"
            ),
            "expected InvBreakpoint FILTER header record, got header:\n{header}"
        );
        for tag in [
            "SVX_START_SUM",
            "SVX_START_SQ_SUM",
            "SVX_END_SUM",
            "SVX_END_SQ_SUM",
            "SVX_SVLEN_SUM",
            "SVX_SVLEN_SQ_SUM",
        ] {
            assert!(
                !header.contains(&format!("##INFO=<ID={tag},")),
                "did not expect {tag} INFO header record, got header:\n{header}"
            );
        }
    }

    #[test]
    fn create_output_header_declares_ci_fields_and_omits_legacy_stats() {
        let vcf_readers = VcfReaders {
            readers: Vec::new(),
            n: 0,
            io_tpool: None,
        };
        let args = merge_args_for_header_test(true);

        let (out_header, _) = create_output_header(&vcf_readers, &args).unwrap();
        let header = header_text(&out_header);

        assert!(
            header.contains(
                "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for merged variants\">"
            ),
            "expected CIPOS INFO header record, got header:\n{header}"
        );
        assert!(
            header.contains(
                "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for merged variants\">"
            ),
            "expected CIEND INFO header record, got header:\n{header}"
        );
        for tag in [
            "START_AVG",
            "START_VARIANCE",
            "SVLEN_AVG",
            "SVLEN_VARIANCE",
            "AVG_START",
            "AVG_END",
        ] {
            assert!(
                !header.contains(&format!("##INFO=<ID={tag},")),
                "did not expect {tag} INFO header record, got header:\n{header}"
            );
        }
    }

    fn common_sv_info_header_records(header: &mut rust_htslib::bcf::Header) {
        header.push_record(br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=END,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=CIEND,Number=2,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=IDLIST,Number=.,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SVCLAIM,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SVCLAIM_SET,Number=.,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=TR_CONTAINED,Number=1,Type=String,Description="">"#);
    }

    fn common_sv_format_header_records(header: &mut rust_htslib::bcf::Header) {
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="">"#);
        header.push_record(br#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="">"#);
    }

    fn common_cnv_format_header_records(header: &mut rust_htslib::bcf::Header) {
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##FORMAT=<ID=CN,Number=1,Type=Float,Description="">"#);
        header.push_record(br#"##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="">"#);
    }

    fn common_bnd_info_header_records(header: &mut rust_htslib::bcf::Header) {
        header.push_record(br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=EVENT,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=EVENTTYPE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=MATEID,Number=.,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=CHR2,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=STRANDS,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=IDLIST,Number=.,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=TR_CONTAINED,Number=1,Type=String,Description="">"#);
    }

    fn base_header(contigs: &[&str], samples: &[&str]) -> rust_htslib::bcf::Header {
        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br#"##fileformat=VCFv4.2"#);
        for contig in contigs {
            let record = format!("##contig=<ID={contig}>");
            header.push_record(record.as_bytes());
        }
        for sample in samples {
            header.push_sample(sample.as_bytes());
        }
        header
    }

    fn sv_header(contigs: &[&str], samples: &[&str]) -> rust_htslib::bcf::Header {
        let mut header = base_header(contigs, samples);
        common_sv_info_header_records(&mut header);
        common_sv_format_header_records(&mut header);
        header
    }

    fn cnv_header(contigs: &[&str], samples: &[&str]) -> rust_htslib::bcf::Header {
        let mut header = base_header(contigs, samples);
        common_sv_info_header_records(&mut header);
        common_cnv_format_header_records(&mut header);
        header
    }

    fn bnd_header(contigs: &[&str], samples: &[&str]) -> rust_htslib::bcf::Header {
        let mut header = base_header(contigs, samples);
        header.push_record(br#"##FILTER=<ID=InvBreakpoint,Description="">"#);
        common_bnd_info_header_records(&mut header);
        common_sv_format_header_records(&mut header);
        header
    }

    fn single_sample_sv_header() -> rust_htslib::bcf::Header {
        sv_header(&["chr1"], &["sample1"])
    }

    fn find_sample_field<'a>(format_keys: &'a str, sample_data: &'a str, field: &str) -> &'a str {
        let format_fields: Vec<&str> = format_keys.split(':').collect();
        let sample_fields: Vec<&str> = sample_data.split(':').collect();
        let field_idx = format_fields
            .iter()
            .position(|format_field| *format_field == field)
            .unwrap();
        sample_fields[field_idx]
    }

    fn gt_from_array(gt: [i32; 2]) -> Vec<GenotypeAllele> {
        vec![
            GenotypeAllele::Unphased(gt[0]),
            GenotypeAllele::Unphased(gt[1]),
        ]
    }

    fn make_del_sv_vcf(gq: i32, gt: [i32; 2]) -> VcfWriteData {
        VcfWriteData {
            rid: Some(0),
            pos: 0,
            alleles: vec![b"A".to_vec(), b"<DEL>".to_vec()],
            gt: gt_from_array(gt),
            sample_format: SampleFormatData::Sv(SvSampleFormatData {
                gq,
                pl: vec![gq, gq + 1, gq + 2],
                ad: vec![gq, gq + 3],
            }),
            format_data: None,
        }
    }

    fn make_cnv_vcf(cn: f32, cnq: f32, gt: [i32; 2]) -> VcfWriteData {
        make_cnv_vcf_with_alt(cn, cnq, gt, b"<CNV>")
    }

    fn make_cnv_vcf_with_alt(cn: f32, cnq: f32, gt: [i32; 2], alt: &[u8]) -> VcfWriteData {
        VcfWriteData {
            rid: Some(0),
            pos: 0,
            alleles: vec![b"A".to_vec(), alt.to_vec()],
            gt: gt_from_array(gt),
            sample_format: SampleFormatData::Cnv(CnvSampleFormatData { cn, cnq }),
            format_data: None,
        }
    }

    fn output_non_header_lines(path: &str) -> Vec<String> {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_owned)
            .collect()
    }

    fn info_field(record_line: &str) -> &str {
        record_line.split('\t').nth(7).unwrap()
    }

    fn info_contains_field(info: &str, expected_field: &str) -> bool {
        info.split(';').any(|field| field == expected_field)
    }

    fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
        info.split(';').find_map(|field| {
            field
                .split_once('=')
                .and_then(|(k, v)| if k == key { Some(v) } else { None })
        })
    }

    #[test]
    fn aggregate_variance_preserves_small_input_variance_for_large_coordinates() {
        let mut group = Vec::new();
        for idx in 0..28 {
            let mut variant =
                test_utils::from_parts(idx, format!("v{idx}"), SvType::DELETION, 0.0, -100.0)
                    .unwrap();
            variant.support_calls = 1;
            variant.start_mean = 29_660_535.0;
            variant.start_variance = 0.01;
            group.push(variant);
        }

        let aggregates = accumulate_group_aggregates(&group).unwrap();
        let mean_start =
            aggregate_mean(aggregates.start_sum, aggregates.support_calls, "AVG_START").unwrap();
        let start_variance = aggregate_variance(
            &group,
            aggregates.support_calls,
            mean_start,
            |variant| variant.start_mean,
            |variant| variant.start_variance,
            "START_VARIANCE",
        )
        .unwrap();

        assert!(
            (start_variance - 0.01).abs() < 1e-12,
            "expected aggregate variance close to 0.01, got {start_variance}",
        );
    }

    #[test]
    fn bnd_output_writes_paired_bnd_records_with_mateid() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr9>
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tbnd_a\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
chr1\t200\tbnd_b\tA\t]chr9:100]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
";
        let in_path = make_temp_vcf(vcf);
        let mut reader = rust_htslib::bcf::Reader::from_path(&in_path).unwrap();
        let mut it = reader.records();
        let r0 = it.next().unwrap().unwrap();
        let r1 = it.next().unwrap().unwrap();
        let args = MergeArgsInner::default();
        let bnd_a = parse_variant_record(&r0, 0, "chr9", &args);
        let bnd_b = parse_variant_record(&r1, 0, "chr1", &args);
        let event = VariantInternal::from_bnd_pair(bnd_a, bnd_b, &args).unwrap();

        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.reverse_map.insert(0, (0, 0));

        let header = bnd_header(&["chr9", "chr1"], &["sample1"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(
            &header,
            &Some(OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }),
            Some(&out_path),
            None,
        )
        .unwrap();

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![event]],
            contig: "chr1_chr9_TRA".to_string(),
            variant_type: SvType::BND,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap();
        drop(writer);

        let out = fs::read_to_string(&out_path).unwrap();
        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 2);
        assert!(out.contains("SVTYPE=BND"));
        assert!(out.contains("MATEID="));
        assert!(out.contains("EVENT="));
        assert!(out.contains("CHR2="));
        assert!(out.contains("STRANDS="));

        let records: Vec<Vec<&str>> = non_header_lines
            .iter()
            .map(|line| line.split('\t').collect())
            .collect();
        assert_eq!(records.len(), 2);
        for idx in 0..2 {
            let mate_idx = 1 - idx;
            let id = records[idx][2];
            let mate_id = records[mate_idx][2];
            let mate_chrom = records[mate_idx][0];
            let info = records[idx][7];

            assert!(
                info.contains(&format!("MATEID={mate_id}")),
                "record {id} must point to mate {mate_id} via MATEID"
            );
            assert!(
                info.contains(&format!("CHR2={mate_chrom}")),
                "record {id} should carry CHR2 matching mate chrom"
            );
            assert!(
                !info.split(';').any(|field| field.starts_with("END=")),
                "record {id} should omit END"
            );
            assert!(
                info.contains("STRANDS=--"),
                "record {id} should carry STRANDS"
            );
        }
    }

    #[test]
    fn write_variants_non_bnd_writes_ci_fields_for_aggregated_positions() {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.index_map.insert((1, 0), 1);
        sample_mapping.index_map.insert((2, 0), 2);
        sample_mapping.reverse_map.insert(0, (0, 0));
        sample_mapping.reverse_map.insert(1, (1, 0));
        sample_mapping.reverse_map.insert(2, (2, 0));

        let header = sv_header(&["chr1"], &["sample1", "sample2", "sample3"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let vcf = make_del_sv_vcf(0, [0, 1]);

        let mut v0 =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        v0.vcf = Some(vcf.clone());
        let mut v1 =
            test_utils::from_parts(1, "v1".to_string(), SvType::DELETION, 1.0, 2.0).unwrap();
        v1.vcf = Some(vcf.clone());
        let mut v2 =
            test_utils::from_parts(2, "v2".to_string(), SvType::DELETION, 2.0, 3.0).unwrap();
        v2.vcf = Some(vcf);

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v1, v2]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 3, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let info = info_field(&non_header_lines[0]);
        assert!(
            info_contains_field(info, "CIPOS=-1,3"),
            "expected CIPOS to reflect aggregated start uncertainty, got INFO={info}"
        );
        assert!(
            info_contains_field(info, "CIEND=-1,3"),
            "expected CIEND to reflect aggregated end uncertainty, got INFO={info}"
        );
        assert!(
            info_contains_field(info, "SVLEN=-2"),
            "expected merged SVLEN to be derived from aggregated length, got INFO={info}"
        );
        assert!(
            info_contains_field(info, "END=3"),
            "expected END to be derived from aggregated breakpoints, got INFO={info}"
        );
        for tag in ["START_AVG", "START_VARIANCE", "SVLEN_AVG", "SVLEN_VARIANCE"] {
            assert!(
                info_value(info, tag).is_none(),
                "did not expect {tag} in INFO={info}"
            );
        }
    }

    #[test]
    fn write_variants_non_bnd_writes_ci_and_merged_svlen_without_legacy_stats() {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.index_map.insert((1, 0), 1);
        sample_mapping.index_map.insert((2, 0), 2);
        sample_mapping.reverse_map.insert(0, (0, 0));
        sample_mapping.reverse_map.insert(1, (1, 0));
        sample_mapping.reverse_map.insert(2, (2, 0));
        let header = sv_header(&["chr1"], &["sample0", "sample1", "sample2"]);
        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut v0 =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, -10.0).unwrap();
        v0.vcf = Some(make_del_sv_vcf(20, [0, 1]));
        let mut v1 =
            test_utils::from_parts(1, "v1".to_string(), SvType::DELETION, 0.0, -20.0).unwrap();
        v1.vcf = Some(make_del_sv_vcf(20, [0, 1]));
        let mut v2 =
            test_utils::from_parts(2, "v2".to_string(), SvType::DELETION, 0.0, -30.0).unwrap();
        v2.vcf = Some(make_del_sv_vcf(20, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v1, v2]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 3, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);

        assert!(
            info_contains_field(info, "CIPOS=0,0"),
            "expected CIPOS to be written, got INFO={info}"
        );
        assert!(
            info_contains_field(info, "CIEND=-16,16"),
            "expected CIEND to be written, got INFO={info}"
        );
        assert!(
            info_contains_field(info, "SVLEN=-20"),
            "expected merged SVLEN semantics, got INFO={info}"
        );
        for tag in [
            "START_AVG",
            "START_VARIANCE",
            "SVLEN_AVG",
            "SVLEN_VARIANCE",
            "AVG_START",
            "AVG_END",
        ] {
            assert!(
                info_value(info, tag).is_none(),
                "did not expect {tag} in INFO={info}"
            );
        }
    }

    #[test]
    fn bnd_output_writes_cipos_and_omits_avg_fields() {
        let sample_mapping = sample_mapping_two_samples();
        let header = bnd_header(&["chr1", "chr2"], &["sample1", "sample2"]);
        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();
        let make_sv_format = |gq: i32| {
            SampleFormatData::Sv(SvSampleFormatData {
                gq,
                pl: vec![gq, gq + 1, gq + 2],
                ad: vec![gq, gq + 3],
            })
        };

        let make_bnd_variant = |vcf_id: usize, id: &str, gq: i32, gt: [i32; 2], a_pos0: i64| {
            let mut v =
                test_utils::from_parts(vcf_id, id.to_string(), SvType::BND, 100.0, 100.0).unwrap();
            v.support_calls = 1;
            v.start_mean = (a_pos0 + 1) as f64;
            v.start_variance = 25.0;
            v.end = 200.0;
            v.svlen = 0.0;
            v.svlen_mean = 0.0;
            v.svlen_variance = 0.0;
            v.bnd_event = Some(BndEventData {
                a_contig: "chr1".to_string(),
                a_pos0,
                a_strands: *b"++",
                a_ref_base: b'A',
                a_vcf: VcfWriteData {
                    rid: Some(0),
                    pos: a_pos0,
                    alleles: vec![b"A".to_vec(), b"A[chr2:201[".to_vec()],
                    gt: gt_from_array(gt),
                    sample_format: make_sv_format(gq),
                    format_data: None,
                },
                b_contig: "chr2".to_string(),
                b_pos0: 200,
                b_strands: *b"-+",
                b_ref_base: b'C',
                b_vcf: VcfWriteData {
                    rid: Some(1),
                    pos: 200,
                    alleles: vec![b"C".to_vec(), b"C]chr1:101]".to_vec()],
                    gt: gt_from_array(gt),
                    sample_format: make_sv_format(gq),
                    format_data: None,
                },
                inv_event_id: None,
            });
            v
        };

        let v0 = make_bnd_variant(0, "bnd0", 10, [0, 1], 100);
        let v1 = make_bnd_variant(1, "bnd1", 20, [1, 1], 110);

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v1]],
            contig: "chr1_chr2_TRA".to_string(),
            variant_type: SvType::BND,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 2);
        for line in &non_header_lines {
            let info = info_field(line);
            assert!(
                info.contains("CIPOS="),
                "expected BND record to write CIPOS, got INFO={info}"
            );
            assert!(
                !info.contains("AVG_START=") && !info.contains("AVG_END="),
                "did not expect AVG_START/AVG_END, got INFO={info}"
            );
        }
    }

    #[test]
    fn write_variants_returns_error_when_supporting_variant_vcf_mismatches_output_mapping() {
        let sample_mapping = sample_mapping_one_sample();
        let header = single_sample_sv_header();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut v0 =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 10.0, 20.0).unwrap();
        v0.vcf = Some(make_del_sv_vcf(0, [0, 1]));
        let mut v1 =
            test_utils::from_parts(1, "v1".to_string(), SvType::DELETION, 11.0, 21.0).unwrap();
        v1.vcf = Some(make_del_sv_vcf(0, [0, 1]));
        v0.support_mask = vec![0];
        v1.support_mask = vec![1];

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 2,
        };

        let err = write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap_err();
        assert!(
            err.to_string().contains("maps to VCF"),
            "expected sample-mapping mismatch error; got: {err}"
        );
    }

    #[test]
    fn write_variants_uses_block_contig_when_input_rid_points_elsewhere() {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.reverse_map.insert(0, (0, 0));

        let header = sv_header(&["chr9", "chr1"], &["sample1"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut variant =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        variant.vcf = Some(make_del_sv_vcf(0, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let chrom = non_header_lines[0].split('\t').next().unwrap();
        assert_eq!(
            chrom, "chr1",
            "expected output CHROM to match the block contig, not the source RID"
        );
    }

    #[test]
    fn write_variants_populates_dummy_record_for_record_reuse() {
        let sample_mapping = sample_mapping_one_sample();
        let header = single_sample_sv_header();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();
        assert_eq!(writer.dummy_record.id(), b".".to_vec());

        let mut variant =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        variant.vcf = Some(make_del_sv_vcf(0, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap();
        assert_eq!(writer.dummy_record.id(), b"v0".to_vec());
    }

    #[test]
    fn tr_contained_info_is_omitted_when_group_has_no_tr_annotation() {
        let sample_mapping = sample_mapping_one_sample();
        let header = single_sample_sv_header();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut variant =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        variant.vcf = Some(make_del_sv_vcf(0, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let info = info_field(&non_header_lines[0]);
        assert!(
            !info.contains("TR_CONTAINED="),
            "expected TR_CONTAINED to be omitted when no grouped variant has TR annotation; got INFO={info}"
        );
    }

    #[test]
    fn tr_contained_info_is_written_when_group_has_tr_annotation() {
        let sample_mapping = sample_mapping_one_sample();
        let header = single_sample_sv_header();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut variant =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        variant.trid = Some(TrId {
            id: "TR1".to_string(),
            motif_len: 2,
        });
        variant.vcf = Some(make_del_sv_vcf(0, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let info = info_field(&non_header_lines[0]);
        assert!(
            info.contains("TR_CONTAINED=TR1"),
            "expected TR_CONTAINED to carry the TR ID when grouped variant has TR annotation; got INFO={info}"
        );
    }

    #[test]
    fn write_variants_intrasample_sv_uses_sample_support_and_sample_matched_format() {
        let sample_mapping = sample_mapping_two_samples();
        let header = sv_header(&["chr1"], &["sample1", "sample2"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut v0 =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 10.0, 10.0).unwrap();
        v0.vcf = Some(make_del_sv_vcf(10, [0, 1]));
        let mut v0_dup =
            test_utils::from_parts(0, "v0_dup".to_string(), SvType::DELETION, 11.0, 11.0).unwrap();
        v0_dup.vcf = Some(make_del_sv_vcf(20, [0, 1]));
        let mut v1 =
            test_utils::from_parts(1, "v1".to_string(), SvType::DELETION, 12.0, 12.0).unwrap();
        v1.vcf = Some(make_del_sv_vcf(30, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v0_dup, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let fields: Vec<&str> = non_header_lines[0].split('\t').collect();
        let info = fields[7];
        assert!(info.contains("SUPP=2"));
        assert!(info.contains("SUPP_CALLS=3"));
        assert!(!info.contains("SUPP_VEC="));

        let sample2 = fields[10];
        let sample2_gq = find_sample_field(fields[8], sample2, "GQ");
        assert_eq!(
            sample2_gq, "30",
            "expected sample2 FORMAT to come from sample2 variant, not duplicate sample1 variant"
        );
    }

    #[test]
    fn write_variants_preserves_nested_support_calls_from_input_info() {
        let sample_mapping = sample_mapping_two_samples();
        let header = sv_header(&["chr1"], &["sample1", "sample2"]);
        let args = MergeArgsInner::default();

        let parent_merged_vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=IDLIST,Number=.,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\tmerged_parent\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;SUPP_CALLS=2;IDLIST=a,b\tGT:GQ:PL:AD\t0/1:20:20,0,30:8,4
";
        let child_raw_vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample2
chr1\t100\traw_child\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10\tGT:GQ:PL:AD\t0/1:30:30,0,40:9,5
";

        let parent_path = make_temp_vcf(parent_merged_vcf);
        let child_path = make_temp_vcf(child_raw_vcf);

        let mut parent_reader = rust_htslib::bcf::Reader::from_path(&parent_path).unwrap();
        let parent_record = parent_reader.records().next().unwrap().unwrap();
        let mut parent = parse_variant_record(&parent_record, 0, "chr1", &args);
        parent
            .remap_support_mask_to_output(&sample_mapping, 2)
            .unwrap();

        let mut child_reader = rust_htslib::bcf::Reader::from_path(&child_path).unwrap();
        let child_record = child_reader.records().next().unwrap().unwrap();
        let mut child = parse_variant_record(&child_record, 1, "chr1", &args);
        child
            .remap_support_mask_to_output(&sample_mapping, 2)
            .unwrap();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![parent, child]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);
        assert!(info.contains("SUPP=2"));
        assert!(
            info.contains("SUPP_CALLS=3"),
            "expected nested support calls (2 + 1) to be preserved, got INFO={info}"
        );
        assert!(
            info.contains("IDLIST=0_a,0_b,1_raw_child"),
            "expected nested IDLIST entries from parent plus child ID, got INFO={info}"
        );
    }

    #[test]
    fn write_variants_idlist_is_stable_between_raw_and_nested_inputs() {
        let sample_mapping = sample_mapping_two_samples();
        let header = sv_header(&["chr1"], &["sample1", "sample2"]);

        let write_group_and_get_idlist = |group: Vec<VariantInternal>| {
            let out_path = make_temp_out("vcf");
            let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();
            let vbr = VariantBlockResult {
                blob_ordinal: 0,
                groups: vec![group],
                contig: "chr1".to_string(),
                variant_type: SvType::DELETION,
                n: 2,
            };
            write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
            drop(writer);

            let lines = output_non_header_lines(&out_path);
            assert_eq!(lines.len(), 1);
            let info = info_field(&lines[0]);
            info_value(info, "IDLIST").unwrap().to_string()
        };

        let mut raw_v0 =
            test_utils::from_parts(0, "raw0".to_string(), SvType::DELETION, 100.0, 100.0).unwrap();
        raw_v0.vcf = Some(make_del_sv_vcf(10, [0, 1]));
        let mut raw_v1 =
            test_utils::from_parts(1, "raw1".to_string(), SvType::DELETION, 100.0, 100.0).unwrap();
        raw_v1.vcf = Some(make_del_sv_vcf(20, [0, 1]));

        let mut nested_v0 =
            test_utils::from_parts(0, "p1_merged".to_string(), SvType::DELETION, 100.0, 100.0)
                .unwrap();
        nested_v0.id_list = vec!["raw1".to_string()];
        nested_v0.vcf = Some(make_del_sv_vcf(10, [0, 1]));

        let mut nested_v1 =
            test_utils::from_parts(1, "p2_merged".to_string(), SvType::DELETION, 100.0, 100.0)
                .unwrap();
        nested_v1.id_list = vec!["raw0".to_string()];
        nested_v1.vcf = Some(make_del_sv_vcf(20, [0, 1]));

        let raw_idlist = write_group_and_get_idlist(vec![raw_v0, raw_v1]);
        let nested_idlist = write_group_and_get_idlist(vec![nested_v0, nested_v1]);
        assert_eq!(
            raw_idlist, nested_idlist,
            "IDLIST should be stable between one-shot raw calls and incremental nested calls"
        );
    }

    #[test]
    fn write_variants_preserves_nested_ci_fields_from_input_info() {
        let sample_mapping = sample_mapping_two_samples();
        let header = sv_header(&["chr1"], &["sample1", "sample2"]);
        let args = MergeArgsInner::default();

        let parent_merged_vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=IDLIST,Number=.,Type=String,Description=\"\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t120\tmerged_parent\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-20;SUPP_CALLS=2;IDLIST=a,b;END=140;CIPOS=-20,20;CIEND=-10,30\tGT:GQ:PL:AD\t0/1:20:20,0,30:8,4
";
        let child_raw_vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample2
chr1\t130\traw_child\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-8\tGT:GQ:PL:AD\t0/1:30:30,0,40:9,5
";

        let parent_path = make_temp_vcf(parent_merged_vcf);
        let child_path = make_temp_vcf(child_raw_vcf);

        let mut parent_reader = rust_htslib::bcf::Reader::from_path(&parent_path).unwrap();
        let parent_record = parent_reader.records().next().unwrap().unwrap();
        let mut parent = parse_variant_record(&parent_record, 0, "chr1", &args);
        parent
            .remap_support_mask_to_output(&sample_mapping, 2)
            .unwrap();

        let mut child_reader = rust_htslib::bcf::Reader::from_path(&child_path).unwrap();
        let child_record = child_reader.records().next().unwrap().unwrap();
        let mut child = parse_variant_record(&child_record, 1, "chr1", &args);
        child
            .remap_support_mask_to_output(&sample_mapping, 2)
            .unwrap();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![parent, child]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);
        assert!(
            info.contains("CIPOS="),
            "expected nested output to carry CIPOS, got INFO={info}"
        );
        assert!(
            info.contains("CIEND="),
            "expected nested output to carry CIEND, got INFO={info}"
        );
        assert!(
            info.contains("END="),
            "expected nested output to carry END, got INFO={info}"
        );
        assert!(
            info.contains("SVLEN="),
            "expected nested output to carry merged SVLEN, got INFO={info}"
        );
        for tag in ["START_AVG", "START_VARIANCE", "SVLEN_AVG", "SVLEN_VARIANCE"] {
            assert!(
                info_value(info, tag).is_none(),
                "did not expect {tag} in INFO={info}"
            );
        }
        for tag in [
            "SVX_START_SUM",
            "SVX_START_SQ_SUM",
            "SVX_END_SUM",
            "SVX_END_SQ_SUM",
            "SVX_SVLEN_SUM",
            "SVX_SVLEN_SQ_SUM",
        ] {
            assert!(
                !info.contains(&format!("{tag}=")),
                "did not expect {tag} in INFO={info}"
            );
        }
    }

    #[test]
    fn write_variants_preserves_nested_svclaim_set_from_input_info() {
        let sample_mapping = sample_mapping_two_samples();
        let header = cnv_header(&["chr1"], &["sample1", "sample2"]);
        let args = MergeArgsInner::default();

        let parent_merged_vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=IDLIST,Number=.,Type=String,Description=\"\">
##INFO=<ID=SVCLAIM_SET,Number=.,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\tmerged_parent\tA\t<CNV>\t.\tPASS\tSVTYPE=CNV;SVLEN=10;SUPP_CALLS=2;IDLIST=a,b;SVCLAIM_SET=D,J\tGT:CN:CNQ\t0/1:2:10
";
        let child_raw_vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=SVCLAIM,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample2
chr1\t100\traw_child\tA\t<CNV>\t.\tPASS\tSVTYPE=CNV;SVLEN=10;SVCLAIM=D\tGT:CN:CNQ\t0/1:3:11
";

        let parent_path = make_temp_vcf(parent_merged_vcf);
        let child_path = make_temp_vcf(child_raw_vcf);

        let mut parent_reader = rust_htslib::bcf::Reader::from_path(&parent_path).unwrap();
        let parent_record = parent_reader.records().next().unwrap().unwrap();
        let mut parent = parse_variant_record(&parent_record, 0, "chr1", &args);
        parent
            .remap_support_mask_to_output(&sample_mapping, 2)
            .unwrap();

        let mut child_reader = rust_htslib::bcf::Reader::from_path(&child_path).unwrap();
        let child_record = child_reader.records().next().unwrap().unwrap();
        let mut child = parse_variant_record(&child_record, 1, "chr1", &args);
        child
            .remap_support_mask_to_output(&sample_mapping, 2)
            .unwrap();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![parent, child]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);
        assert!(
            info.contains("SVCLAIM_SET=D,J"),
            "expected nested SVCLAIM_SET values to be preserved, got INFO={info}"
        );
    }

    #[test]
    fn write_variants_cnv_stage_output_stays_canonical_and_incremental_matches_one_shot() {
        let sample_mapping = sample_mapping_one_sample();
        let header = cnv_header(&["chr1"], &["sample1"]);
        let args = MergeArgsInner::default();

        let make_del_like_cnv_variant = |id: &str, svlen: f64, cn: f32, cnq: f32| {
            let mut variant =
                test_utils::from_parts(0, id.to_string(), SvType::DELETION, 100.0, svlen).unwrap();
            variant.vcf = Some(make_cnv_vcf_with_alt(cn, cnq, [0, 1], b"<DEL>"));
            variant
        };
        let make_cnv_variant = || {
            let mut variant = test_utils::from_parts(
                0,
                "sawfish:CNV:0:0:300".to_string(),
                SvType::CNV,
                100.0,
                400.0,
            )
            .unwrap();
            variant.vcf = Some(make_cnv_vcf(3.0, 30.0, [0, 1]));
            variant
        };

        let one_shot_out_path = make_temp_out("vcf");
        let mut one_shot_writer =
            VcfWriter::new(&header, &None, Some(&one_shot_out_path), None).unwrap();
        let one_shot_group = vec![
            make_del_like_cnv_variant("sawfish:CNV:0:0:100", 100.0, 1.0, 10.0),
            make_del_like_cnv_variant("sawfish:CNV:0:0:200", 100.0, 2.0, 20.0),
            make_cnv_variant(),
        ];
        let one_shot_vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![one_shot_group],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 3,
        };
        write_variants(
            one_shot_vbr,
            &sample_mapping,
            1,
            1,
            false,
            &mut one_shot_writer,
        )
        .unwrap();
        drop(one_shot_writer);

        let one_shot_lines = output_non_header_lines(&one_shot_out_path);
        assert_eq!(one_shot_lines.len(), 1);
        let one_shot_fields: Vec<&str> = one_shot_lines[0].split('\t').collect();
        let one_shot_alt = one_shot_fields[4];
        let one_shot_info = info_field(&one_shot_lines[0]);
        let one_shot_svtype = info_value(one_shot_info, "SVTYPE").unwrap();
        let one_shot_svlen = info_value(one_shot_info, "SVLEN").unwrap();
        let one_shot_end = info_value(one_shot_info, "END").unwrap();

        let stage1_out_path = make_temp_out("vcf");
        let mut stage1_writer =
            VcfWriter::new(&header, &None, Some(&stage1_out_path), None).unwrap();
        let stage1_group = vec![
            make_del_like_cnv_variant("sawfish:CNV:0:0:100", 100.0, 1.0, 10.0),
            make_del_like_cnv_variant("sawfish:CNV:0:0:200", 100.0, 2.0, 20.0),
        ];
        let stage1_vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![stage1_group],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };
        write_variants(stage1_vbr, &sample_mapping, 1, 1, false, &mut stage1_writer).unwrap();
        drop(stage1_writer);

        let stage1_lines = output_non_header_lines(&stage1_out_path);
        assert_eq!(stage1_lines.len(), 1);
        let stage1_fields: Vec<&str> = stage1_lines[0].split('\t').collect();
        let stage1_info = info_field(&stage1_lines[0]);
        assert_eq!(
            stage1_fields[4], "<CNV>",
            "CNV block output ALT must stay canonical in incremental intermediates"
        );
        assert_eq!(
            info_value(stage1_info, "SVTYPE"),
            Some("CNV"),
            "CNV block output SVTYPE must stay canonical in incremental intermediates"
        );
        let stage1_svlen = info_value(stage1_info, "SVLEN").unwrap();
        assert!(
            !stage1_svlen.starts_with('-'),
            "CNV block output SVLEN must remain non-directional in incremental intermediates"
        );

        let mut stage1_reader = rust_htslib::bcf::Reader::from_path(&stage1_out_path).unwrap();
        let stage1_record = stage1_reader.records().next().unwrap().unwrap();
        let mut stage1_parent = parse_variant_record(&stage1_record, 0, "chr1", &args);
        stage1_parent
            .remap_support_mask_to_output(&sample_mapping, 1)
            .unwrap();

        let stage2_out_path = make_temp_out("vcf");
        let mut stage2_writer =
            VcfWriter::new(&header, &None, Some(&stage2_out_path), None).unwrap();
        let stage2_vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![stage1_parent, make_cnv_variant()]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };
        write_variants(stage2_vbr, &sample_mapping, 1, 1, false, &mut stage2_writer).unwrap();
        drop(stage2_writer);

        let stage2_lines = output_non_header_lines(&stage2_out_path);
        assert_eq!(stage2_lines.len(), 1);
        let stage2_fields: Vec<&str> = stage2_lines[0].split('\t').collect();
        let stage2_alt = stage2_fields[4];
        let stage2_info = info_field(&stage2_lines[0]);
        let stage2_svtype = info_value(stage2_info, "SVTYPE").unwrap();
        let stage2_svlen = info_value(stage2_info, "SVLEN").unwrap();
        let stage2_end = info_value(stage2_info, "END").unwrap();

        assert_eq!(
            stage2_alt, one_shot_alt,
            "incremental CNV ALT should match one-shot output"
        );
        assert_eq!(
            stage2_svtype, one_shot_svtype,
            "incremental CNV SVTYPE should match one-shot output"
        );
        assert_eq!(
            stage2_svlen, one_shot_svlen,
            "incremental CNV SVLEN should match one-shot output for shared membership"
        );
        assert_eq!(
            stage2_end, one_shot_end,
            "incremental CNV END should match one-shot output for shared membership"
        );
    }

    #[test]
    fn write_variants_intrasample_cnv_uses_sample_support_and_sample_matched_format() {
        let sample_mapping = sample_mapping_two_samples();
        let header = cnv_header(&["chr1"], &["sample1", "sample2"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut v0 =
            test_utils::from_parts(0, "cnv0".to_string(), SvType::CNV, 20.0, 20.0).unwrap();
        v0.vcf = Some(make_cnv_vcf(1.0, 11.0, [0, 1]));
        let mut v0_dup =
            test_utils::from_parts(0, "cnv0_dup".to_string(), SvType::CNV, 21.0, 21.0).unwrap();
        v0_dup.vcf = Some(make_cnv_vcf(2.0, 22.0, [0, 1]));
        let mut v1 =
            test_utils::from_parts(1, "cnv1".to_string(), SvType::CNV, 22.0, 22.0).unwrap();
        v1.vcf = Some(make_cnv_vcf(3.0, 33.0, [1, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v0_dup, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let fields: Vec<&str> = non_header_lines[0].split('\t').collect();
        let info = fields[7];
        assert!(info.contains("SUPP=2"));
        assert!(info.contains("SUPP_CALLS=3"));
        assert!(!info.contains("SUPP_VEC="));

        let sample2 = fields[10];
        let sample2_gt = find_sample_field(fields[8], sample2, "GT");
        let sample2_cn = find_sample_field(fields[8], sample2, "CN")
            .parse::<f32>()
            .unwrap();
        let sample2_cnq = find_sample_field(fields[8], sample2, "CNQ")
            .parse::<f32>()
            .unwrap();
        assert_eq!(sample2_gt, "1/1");
        assert_eq!(sample2_cn, 3.0);
        assert_eq!(sample2_cnq, 33.0);
    }
    #[test]
    fn write_variants_cnv_writes_svclaim_set_for_consensus_cnv_type() {
        let sample_mapping = sample_mapping_two_samples();
        let header = cnv_header(&["chr1"], &["sample1", "sample2"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut v0 =
            test_utils::from_parts(0, "cnv0".to_string(), SvType::CNV, 20.0, 20.0).unwrap();
        v0.svclaim = Some("D".to_string());
        v0.vcf = Some(make_cnv_vcf(1.0, 10.0, [0, 1]));
        let mut v1 =
            test_utils::from_parts(1, "cnv1".to_string(), SvType::CNV, 21.0, 21.0).unwrap();
        v1.svclaim = Some("D".to_string());
        v1.vcf = Some(make_cnv_vcf(1.0, 10.0, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);
        assert!(!info_contains_field(info, "SVCLAIM=D"));
        assert!(info_contains_field(info, "SVCLAIM_SET=D"));
    }

    #[test]
    fn write_variants_cnv_writes_svclaim_set_when_claims_disagree() {
        let sample_mapping = sample_mapping_two_samples();
        let header = cnv_header(&["chr1"], &["sample1", "sample2"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut v0 =
            test_utils::from_parts(0, "cnv0".to_string(), SvType::CNV, 20.0, 20.0).unwrap();
        v0.svclaim = Some("D".to_string());
        v0.vcf = Some(make_cnv_vcf(1.0, 10.0, [0, 1]));
        let mut v1 =
            test_utils::from_parts(1, "cnv1".to_string(), SvType::CNV, 21.0, 21.0).unwrap();
        v1.svclaim = Some("A".to_string());
        v1.vcf = Some(make_cnv_vcf(1.0, 10.0, [0, 1]));

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);
        assert!(!info.contains("SVCLAIM="));
        assert!(
            info.contains("SVCLAIM_SET=A,D"),
            "expected sorted unique mixed claim values, got INFO={info}"
        );
    }
    #[test]
    fn write_variants_intrasample_bnd_uses_sample_support_and_sample_matched_format() {
        let sample_mapping = sample_mapping_two_samples();
        let header = bnd_header(&["chr1", "chr2"], &["sample1", "sample2"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let make_sv_format = |gq: i32| {
            SampleFormatData::Sv(SvSampleFormatData {
                gq,
                pl: vec![gq, gq + 1, gq + 2],
                ad: vec![gq, gq + 3],
            })
        };

        let make_bnd_variant = |vcf_id: usize, id: &str, gq: i32, gt: [i32; 2]| {
            let mut v =
                test_utils::from_parts(vcf_id, id.to_string(), SvType::BND, 100.0, 100.0).unwrap();
            v.bnd_event = Some(BndEventData {
                a_contig: "chr1".to_string(),
                a_pos0: 100,
                a_strands: *b"+-",
                a_ref_base: b'A',
                a_vcf: VcfWriteData {
                    rid: Some(0),
                    pos: 100,
                    alleles: vec![b"A".to_vec(), b"[chr2:201[A".to_vec()],
                    gt: gt_from_array(gt),
                    sample_format: make_sv_format(gq),
                    format_data: None,
                },
                b_contig: "chr2".to_string(),
                b_pos0: 200,
                b_strands: *b"-+",
                b_ref_base: b'C',
                b_vcf: VcfWriteData {
                    rid: Some(1),
                    pos: 200,
                    alleles: vec![b"C".to_vec(), b"C]chr1:101]".to_vec()],
                    gt: gt_from_array(gt),
                    sample_format: make_sv_format(gq),
                    format_data: None,
                },
                inv_event_id: None,
            });
            v
        };

        let v0 = make_bnd_variant(0, "bnd0", 10, [0, 1]);
        let v0_dup = make_bnd_variant(0, "bnd0_dup", 20, [0, 1]);
        let v1 = make_bnd_variant(1, "bnd1", 30, [1, 1]);

        let vbr = VariantBlockResult {
            blob_ordinal: 0,
            groups: vec![vec![v0, v0_dup, v1]],
            contig: "chr1_chr2_TRA".to_string(),
            variant_type: SvType::BND,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 2);

        for line in &non_header_lines {
            let fields: Vec<&str> = line.split('\t').collect();
            let info = fields[7];
            assert!(info.contains("SUPP=2"));
            assert!(info.contains("SUPP_CALLS=3"));
            assert!(!info.contains("SUPP_VEC="));
            for tag in [
                "SVX_START_SUM",
                "SVX_START_SQ_SUM",
                "SVX_END_SUM",
                "SVX_END_SQ_SUM",
                "SVX_SVLEN_SUM",
                "SVX_SVLEN_SQ_SUM",
            ] {
                assert!(
                    !info.contains(&format!("{tag}=")),
                    "did not expect {tag} in INFO={info}"
                );
            }

            let sample2 = fields[10];
            let sample2_gt = find_sample_field(fields[8], sample2, "GT");
            let sample2_gq = find_sample_field(fields[8], sample2, "GQ");
            assert_eq!(sample2_gt, "1/1");
            assert_eq!(sample2_gq, "30");
        }
    }

    #[test]
    fn write_variants_bnd_event_id_is_stable_between_raw_and_nested_inputs() {
        let sample_mapping = sample_mapping_two_samples();
        let header = bnd_header(&["chr1", "chr2"], &["sample1", "sample2"]);

        let make_sv_format = |gq: i32| {
            SampleFormatData::Sv(SvSampleFormatData {
                gq,
                pl: vec![gq, gq + 1, gq + 2],
                ad: vec![gq, gq + 3],
            })
        };

        let make_bnd_variant = |vcf_id: usize, id: &str, gq: i32, id_list: Vec<&str>| {
            let mut v =
                test_utils::from_parts(vcf_id, id.to_string(), SvType::BND, 100.0, 100.0).unwrap();
            v.id_list = id_list.iter().map(|id| (*id).to_string()).collect();
            v.bnd_event = Some(BndEventData {
                a_contig: "chr1".to_string(),
                a_pos0: 100,
                a_strands: *b"+-",
                a_ref_base: b'A',
                a_vcf: VcfWriteData {
                    rid: Some(0),
                    pos: 100,
                    alleles: vec![b"A".to_vec(), b"[chr2:201[A".to_vec()],
                    gt: gt_from_array([0, 1]),
                    sample_format: make_sv_format(gq),
                    format_data: None,
                },
                b_contig: "chr2".to_string(),
                b_pos0: 200,
                b_strands: *b"-+",
                b_ref_base: b'C',
                b_vcf: VcfWriteData {
                    rid: Some(1),
                    pos: 200,
                    alleles: vec![b"C".to_vec(), b"C]chr1:101]".to_vec()],
                    gt: gt_from_array([0, 1]),
                    sample_format: make_sv_format(gq),
                    format_data: None,
                },
                inv_event_id: None,
            });
            v
        };

        let raw_v0 = make_bnd_variant(0, "raw0", 10, vec![]);
        let raw_v1 = make_bnd_variant(1, "raw1", 20, vec![]);
        let nested_v0 = make_bnd_variant(0, "p1_merged", 10, vec!["raw0"]);
        let nested_v1 = make_bnd_variant(1, "p2_merged", 20, vec!["raw1"]);

        let write_group_and_get_event = |group: Vec<VariantInternal>| {
            let out_path = make_temp_out("vcf");
            let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();
            let vbr = VariantBlockResult {
                blob_ordinal: 0,
                groups: vec![group],
                contig: "chr1_chr2_TRA".to_string(),
                variant_type: SvType::BND,
                n: 2,
            };
            write_variants(vbr, &sample_mapping, 2, 1, false, &mut writer).unwrap();
            drop(writer);

            let lines = output_non_header_lines(&out_path);
            assert_eq!(lines.len(), 2);
            let info = info_field(&lines[0]);
            info_value(info, "EVENT").unwrap().to_string()
        };

        let forward_event = write_group_and_get_event(vec![raw_v0, raw_v1]);
        let reverse_event = write_group_and_get_event(vec![nested_v0, nested_v1]);
        assert_eq!(
            forward_event, reverse_event,
            "BND EVENT ID should be stable between one-shot raw calls and incremental nested calls"
        );
    }

    #[test]
    fn write_variants_bnd_ref_alt_are_stable_between_one_shot_and_incremental_shared_keys() {
        let sample_mapping = sample_mapping_one_sample();
        let header = bnd_header(&["chr1", "chr2"], &["sample1"]);

        let make_sv_format = |gq: i32| {
            SampleFormatData::Sv(SvSampleFormatData {
                gq,
                pl: vec![gq, gq + 1, gq + 2],
                ad: vec![gq, gq + 3],
            })
        };

        let make_bnd_variant = |id: &str,
                                id_list: Vec<&str>,
                                a_ref_base: u8,
                                a_anchor: &str,
                                b_ref_base: u8,
                                b_anchor: &str| {
            let mut v =
                test_utils::from_parts(0, id.to_string(), SvType::BND, 100.0, 100.0).unwrap();
            v.id_list = if id_list.is_empty() {
                vec![id.to_string()]
            } else {
                id_list.into_iter().map(ToString::to_string).collect()
            };
            v.bnd_event = Some(BndEventData {
                a_contig: "chr1".to_string(),
                a_pos0: 100,
                a_strands: *b"+-",
                a_ref_base,
                a_vcf: VcfWriteData {
                    rid: Some(0),
                    pos: 100,
                    alleles: vec![
                        vec![a_ref_base],
                        format!("[chr2:201[{a_anchor}").into_bytes(),
                    ],
                    gt: gt_from_array([0, 1]),
                    sample_format: make_sv_format(10),
                    format_data: None,
                },
                b_contig: "chr2".to_string(),
                b_pos0: 200,
                b_strands: *b"-+",
                b_ref_base,
                b_vcf: VcfWriteData {
                    rid: Some(1),
                    pos: 200,
                    alleles: vec![
                        vec![b_ref_base],
                        format!("{b_anchor}]chr1:101]").into_bytes(),
                    ],
                    gt: gt_from_array([0, 1]),
                    sample_format: make_sv_format(10),
                    format_data: None,
                },
                inv_event_id: None,
            });
            v
        };

        let collect_ref_alt_by_stable_key = |group: Vec<VariantInternal>| {
            let out_path = make_temp_out("vcf");
            let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();
            let vbr = VariantBlockResult {
                blob_ordinal: 0,
                groups: vec![group],
                contig: "chr1_chr2_TRA".to_string(),
                variant_type: SvType::BND,
                n: 2,
            };
            write_variants(vbr, &sample_mapping, 1, 1, false, &mut writer).unwrap();
            drop(writer);

            let lines = output_non_header_lines(&out_path);
            assert_eq!(lines.len(), 2);

            let mut id_to_mate = HashMap::new();
            let mut by_key = HashMap::new();
            for line in lines {
                let fields: Vec<&str> = line.split('\t').collect();
                let info = info_field(&line);
                let id = fields[2].to_string();
                let mate_id = info_value(info, "MATEID").unwrap().to_string();
                id_to_mate.insert(id.clone(), mate_id);

                let stable_key = format!(
                    "{}|{}|{}|{}|{}",
                    info_value(info, "IDLIST").unwrap(),
                    info_value(info, "EVENT").unwrap(),
                    fields[0],
                    fields[1],
                    info_value(info, "STRANDS").unwrap()
                );
                by_key.insert(stable_key, (fields[3].to_string(), fields[4].to_string()));
            }

            assert_eq!(id_to_mate.len(), 2);
            for (id, mate_id) in &id_to_mate {
                assert_eq!(
                    id_to_mate.get(mate_id),
                    Some(id),
                    "BND mate reciprocity should remain consistent in output"
                );
            }
            by_key
        };

        let one_shot = collect_ref_alt_by_stable_key(vec![
            make_bnd_variant("raw_a", vec![], b'C', "C", b'T', "T"),
            make_bnd_variant("raw_b", vec![], b'G', "G", b'A', "A"),
        ]);
        let incremental = collect_ref_alt_by_stable_key(vec![
            make_bnd_variant("a_stage", vec!["raw_b"], b'G', "G", b'A', "A"),
            make_bnd_variant("z_stage", vec!["raw_a"], b'C', "C", b'T', "T"),
        ]);

        assert_eq!(one_shot.len(), incremental.len());
        for (key, one_shot_alleles) in one_shot {
            let incremental_alleles = incremental
                .get(&key)
                .unwrap_or_else(|| panic!("missing incremental BND record for stable key {key}"));
            assert_eq!(
                one_shot_alleles, *incremental_alleles,
                "shared BND key should have identical REF/ALT between one-shot and incremental paths"
            );
        }
    }
}
