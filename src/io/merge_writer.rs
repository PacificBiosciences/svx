use crate::{
    cli::MergeArgs,
    core::{svtype::SvType, variant::VariantInternal, variant_block::VariantBlockResult},
    io::{
        vcf_reader::{SampleMapping, VcfReaders},
        vcf_writer::VcfWriter,
    },
    utils::util::{MISSING_FLOAT, MISSING_INTEGER, Result},
};
use rust_htslib::bcf;
use std::{collections::BTreeSet, env};

pub fn create_output_header(
    vcf_readers: &VcfReaders,
    args: &MergeArgs,
) -> Result<(bcf::Header, SampleMapping)> {
    let mut out_header = bcf::Header::new();
    let sample_mapping = vcf_readers.merge_headers(&mut out_header)?;

    // TODO: Refine fields
    out_header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">"#);
    out_header.push_record(
        br#"##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description="Number of input calls merged into this variant">"#,
    );
    out_header.push_record(
        br#"##INFO=<ID=SUPP_VEC,Number=.,Type=Integer,Description="Vector of supporting samples (0/1 per input sample)">"#,
    );
    out_header.push_record(br#"##INFO=<ID=IDLIST,Number=.,Type=String,Description="Variant IDs of input calls merged to make this call">"#);
    out_header.push_record(br#"##INFO=<ID=START_VARIANCE,Number=1,Type=String,Description="Variance of start position for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=START_AVG,Number=1,Type=String,Description="Average start position for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=SVLEN_AVG,Number=1,Type=String,Description="Average SVLEN for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=SVLEN_VARIANCE,Number=1,Type=String,Description="Variance of SVLEN for variants merged into this one">"#);
    out_header.push_record(
        br#"##INFO=<ID=END,Number=1,Type=Integer,Description="Merged end coordinate">"#,
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
    let version_line = format!(
        "##{}Version={}",
        env!("CARGO_PKG_NAME"),
        crate::cli::FULL_VERSION
    );
    out_header.push_record(version_line.as_bytes());

    let command_line = env::args().collect::<Vec<String>>().join(" ");
    let command_line = format!("##{}Command={}", env!("CARGO_PKG_NAME"), command_line);
    out_header.push_record(command_line.as_bytes());
}

fn create_sample_representatives(
    group: &[VariantInternal],
    sample_mapping: &SampleMapping,
    n_samples: usize,
) -> Result<Vec<Option<usize>>> {
    let mut sample_representatives = vec![None; n_samples];

    for (group_idx, variant) in group.iter().enumerate() {
        let merged_pos = merged_sample_position_for_variant(
            variant,
            sample_mapping,
            "building sample representatives",
        )?;

        if sample_representatives[merged_pos].is_none() {
            sample_representatives[merged_pos] = Some(group_idx);
        }
    }

    Ok(sample_representatives)
}

fn merged_sample_position_for_variant(
    variant: &VariantInternal,
    sample_mapping: &SampleMapping,
    action: &str,
) -> Result<usize> {
    sample_mapping
        .index_map
        .get(&(variant.vcf_id, 0))
        .copied()
        .ok_or_else(|| {
            crate::svx_error!(
                "Failed to find sample mapping for variant {} while {} (vcf: {}, sample: {})",
                variant.id,
                action,
                variant.vcf_id,
                variant.sample_id
            )
        })
}

fn support_vector_values(sample_representatives: &[Option<usize>]) -> Vec<i32> {
    sample_representatives
        .iter()
        .map(|sample_variant| if sample_variant.is_some() { 1 } else { 0 })
        .collect()
}

fn sample_support_count(sample_representatives: &[Option<usize>]) -> usize {
    sample_representatives
        .iter()
        .filter(|sample_variant| sample_variant.is_some())
        .count()
}

fn stable_hash(bytes: &[u8]) -> i32 {
    let mut res: i64 = 0;
    const MOD: i64 = 1_000_000_007;
    for &b in bytes {
        res = (res * 17 + i64::from(b)) % MOD;
    }
    res as i32
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

fn round_to_i64(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

fn format_info_float(value: f64) -> String {
    let mut formatted = format!("{value:.2}");
    while formatted.ends_with('0') {
        formatted.pop();
    }
    if formatted.ends_with('.') {
        formatted.push('0');
    }
    if formatted == "-0.0" {
        return "0.0".to_string();
    }
    formatted
}

fn to_info_i32(value: i64, label: &str) -> Result<i32> {
    i32::try_from(value).map_err(|_| {
        crate::svx_error!(
            "Cannot write {label}={value} to INFO as i32: value is outside supported range"
        )
    })
}

fn push_group_support_info(
    out_rec: &mut bcf::Record,
    group: &[VariantInternal],
    sample_representatives: &[Option<usize>],
) -> Result<()> {
    let supp_vec = support_vector_values(sample_representatives);
    let supp = sample_support_count(sample_representatives) as i32;
    let supp_calls = group.len() as i32;
    let id_list: Vec<&[u8]> = group.iter().map(|v| v.id.as_bytes()).collect();

    out_rec.push_info_integer(b"SUPP", &[supp])?;
    out_rec.push_info_integer(b"SUPP_CALLS", &[supp_calls])?;
    out_rec.push_info_integer(b"SUPP_VEC", &supp_vec)?;
    out_rec.push_info_string(b"IDLIST", &id_list)?;

    Ok(())
}

fn maybe_push_end_info(
    out_rec: &mut bcf::Record,
    representative: &VariantInternal,
    representative_pos0: i64,
    variant_type: SvType,
) -> Result<()> {
    if variant_type != SvType::CNV
        && !matches!(
            representative.svtype,
            SvType::DELETION | SvType::INVERSION | SvType::DUPLICATION
        )
    {
        return Ok(());
    }

    let end_1based = representative_pos0 + round_to_i64(representative.svlen.abs()) + 1;
    let end_i32 = to_info_i32(end_1based, "END")?;
    out_rec.push_info_integer(b"END", &[end_i32])?;
    Ok(())
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
        if let Some(claim) = variant.svclaim.as_deref() {
            unique_claims.insert(claim);
            n_claimed += 1;
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
    variant_type: SvType,
    output_svtype: SvType,
    representative_pos0: i64,
) -> Result<()> {
    let representative = &group[0];
    out_rec.push_info_string(b"SVTYPE", &[output_svtype.to_string().as_bytes()])?;
    out_rec.push_info_integer(b"SVLEN", &[representative.svlen as i32])?;
    maybe_push_end_info(out_rec, representative, representative_pos0, variant_type)?;

    let n = group.len() as f64;
    let mean_start = group
        .iter()
        .map(|variant| variant.start + 1.0f64)
        .sum::<f64>()
        / n;
    let start_variance = group
        .iter()
        .map(|variant| {
            let diff = (variant.start + 1.0f64) - mean_start;
            diff * diff
        })
        .sum::<f64>()
        / n;
    let start_variance_s = format_info_float(start_variance);
    let mean_start_s = format_info_float(mean_start);
    out_rec.push_info_string(b"START_VARIANCE", &[start_variance_s.as_bytes()])?;
    out_rec.push_info_string(b"START_AVG", &[mean_start_s.as_bytes()])?;

    let mean_svlen = group.iter().map(|variant| variant.svlen).sum::<f64>() / n;
    let svlen_variance = group
        .iter()
        .map(|variant| {
            let diff = variant.svlen - mean_svlen;
            diff * diff
        })
        .sum::<f64>()
        / n;

    let svlen_variance_s = format_info_float(svlen_variance);
    let mean_svlen_s = format_info_float(mean_svlen);
    out_rec.push_info_string(b"SVLEN_VARIANCE", &[svlen_variance_s.as_bytes()])?;
    out_rec.push_info_string(b"SVLEN_AVG", &[mean_svlen_s.as_bytes()])?;

    push_group_support_info(out_rec, group, sample_representatives)?;
    push_svclaim_info(out_rec, group, variant_type, output_svtype)?;

    if let Some(tr_id) = group
        .iter()
        .find_map(|variant| variant.trid.as_ref().map(|tr| tr.id.as_str()))
    {
        out_rec.push_info_string(b"TR_CONTAINED", &[tr_id.as_bytes()])?;
    }

    Ok(())
}

fn derive_cnv_output_svtype(group: &[VariantInternal]) -> SvType {
    let mut has_del = false;
    let mut has_dup = false;
    let mut has_other = false;

    for variant in group {
        match variant.svtype {
            SvType::DELETION => has_del = true,
            SvType::DUPLICATION => has_dup = true,
            SvType::CNV | SvType::INSERTION | SvType::INVERSION | SvType::BND => has_other = true,
        }
    }

    if has_other || (has_del && has_dup) {
        SvType::CNV
    } else if has_del {
        SvType::DELETION
    } else if has_dup {
        SvType::DUPLICATION
    } else {
        SvType::CNV
    }
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

fn write_bnd_group(
    group: &[VariantInternal],
    sample_representatives: &[Option<usize>],
    writer: &mut bcf::Writer,
    out_rec: &mut bcf::Record,
) -> Result<()> {
    if group.is_empty() {
        return Ok(());
    }

    let rep = group[0].bnd_event.as_ref().ok_or_else(|| {
        crate::svx_error!(
            "BND group representative {} is missing BND event payload",
            group[0].id
        )
    })?;

    let (a_contig, b_contig) = (&rep.a_contig, &rep.b_contig);
    let rep_a_strands = rep.a_strands;
    let rep_b_strands = rep.b_strands;
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
        if e.a_strands != rep_a_strands || e.b_strands != rep_b_strands {
            let rep_a = std::str::from_utf8(&rep_a_strands).unwrap_or("<non-utf8>");
            let rep_b = std::str::from_utf8(&rep_b_strands).unwrap_or("<non-utf8>");
            let e_a = std::str::from_utf8(&e.a_strands).unwrap_or("<non-utf8>");
            let e_b = std::str::from_utf8(&e.b_strands).unwrap_or("<non-utf8>");
            return Err(crate::svx_error!(
                "BND group strands mismatch for contigs {a_contig}/{b_contig}: expected {rep_a}/{rep_b}, got {e_a}/{e_b} (event {})",
                v.id
            ));
        }
    }

    let mean_a_pos0 = group
        .iter()
        .map(|v| v.bnd_event.as_ref().unwrap().a_pos0 as f64)
        .sum::<f64>()
        / group.len() as f64;
    let mean_b_pos0 = group
        .iter()
        .map(|v| v.bnd_event.as_ref().unwrap().b_pos0 as f64)
        .sum::<f64>()
        / group.len() as f64;
    let a_pos0 = round_to_i64(mean_a_pos0);
    let b_pos0 = round_to_i64(mean_b_pos0);

    let inversion_event_id = bnd_inversion_event_consensus(group);
    let event_id = if let Some(inv_event_id) = inversion_event_id {
        inv_event_id.to_string()
    } else {
        // Deterministic event ID: hash grouped event IDs in write order.
        let id_concat = group
            .iter()
            .map(|v| v.id.as_str())
            .collect::<Vec<_>>()
            .join("|");
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

    let rid_a = writer
        .header()
        .name2rid(a_contig.as_bytes())
        .map_err(|e| crate::svx_error!("Output header missing contig {a_contig}: {e}"))?;
    let rid_b = writer
        .header()
        .name2rid(b_contig.as_bytes())
        .map_err(|e| crate::svx_error!("Output header missing contig {b_contig}: {e}"))?;

    // Record A
    {
        out_rec.clear();
        out_rec.set_rid(Some(rid_a));
        out_rec.set_pos(a_pos0);
        out_rec.set_id(rec_a_id.as_bytes())?;

        let ref_base = rep.a_ref_base;
        let a_template_alt = rep.a_vcf.alleles.get(1).ok_or_else(|| {
            crate::svx_error!(
                "BND event {} is missing ALT allele for breakend A",
                group[0].id
            )
        })?;
        let a_anchor_seq = bnd_anchor_seq_from_alt(a_template_alt.as_slice())?;
        let alt = bnd_alt(b_contig, b_pos0 + 1, rep.a_strands, &a_anchor_seq)?;
        out_rec.set_alleles(&[&[ref_base], &alt])?;

        push_bnd_group_info(
            out_rec,
            group,
            sample_representatives,
            BndRecordInfo {
                event_id: &event_id,
                mate_id: &rec_b_id,
                mate_contig: b_contig,
                strands: rep.a_strands,
                is_inversion_event,
            },
        )?;

        let mut all_gts = Vec::new();
        let mut all_gqs = Vec::new();
        let mut all_pls = Vec::new();
        let mut all_ads = Vec::new();

        for sample_variant in sample_representatives {
            if let Some(group_idx) = sample_variant {
                let event = group[*group_idx].bnd_event.as_ref().unwrap();
                let sv_format = event.a_vcf.sv_format().ok_or_else(|| {
                    crate::svx_error!(
                        "BND event {} is missing SV sample FORMAT payload",
                        group[*group_idx].id
                    )
                })?;
                all_gts.extend_from_slice(&event.a_vcf.gt);
                all_gqs.push(sv_format.gq);
                all_pls.extend_from_slice(&sv_format.pl);
                all_ads.extend_from_slice(&sv_format.ad);
            } else {
                all_gts.extend_from_slice(&[
                    bcf::record::GenotypeAllele::UnphasedMissing,
                    bcf::record::GenotypeAllele::UnphasedMissing,
                ]);
                all_gqs.push(MISSING_INTEGER);
                all_pls.extend_from_slice(&[MISSING_INTEGER, MISSING_INTEGER, MISSING_INTEGER]);
                all_ads.extend_from_slice(&[MISSING_INTEGER, MISSING_INTEGER]);
            }
        }

        out_rec.push_genotypes(&all_gts)?;
        out_rec.push_format_integer(b"GQ", &all_gqs)?;
        out_rec.push_format_integer(b"PL", &all_pls)?;
        out_rec.push_format_integer(b"AD", &all_ads)?;
        writer.write(out_rec)?;
    }

    // Record B
    {
        out_rec.clear();
        out_rec.set_rid(Some(rid_b));
        out_rec.set_pos(b_pos0);
        out_rec.set_id(rec_b_id.as_bytes())?;

        let ref_base = rep.b_ref_base;
        let b_template_alt = rep.b_vcf.alleles.get(1).ok_or_else(|| {
            crate::svx_error!(
                "BND event {} is missing ALT allele for breakend B",
                group[0].id
            )
        })?;
        let b_anchor_seq = bnd_anchor_seq_from_alt(b_template_alt.as_slice())?;
        let alt = bnd_alt(a_contig, a_pos0 + 1, rep.b_strands, &b_anchor_seq)?;
        out_rec.set_alleles(&[&[ref_base], &alt])?;

        push_bnd_group_info(
            out_rec,
            group,
            sample_representatives,
            BndRecordInfo {
                event_id: &event_id,
                mate_id: &rec_a_id,
                mate_contig: a_contig,
                strands: rep.b_strands,
                is_inversion_event,
            },
        )?;

        let mut all_gts = Vec::new();
        let mut all_gqs = Vec::new();
        let mut all_pls = Vec::new();
        let mut all_ads = Vec::new();

        for sample_variant in sample_representatives {
            if let Some(group_idx) = sample_variant {
                let event = group[*group_idx].bnd_event.as_ref().unwrap();
                let sv_format = event.b_vcf.sv_format().ok_or_else(|| {
                    crate::svx_error!(
                        "BND event {} is missing SV sample FORMAT payload",
                        group[*group_idx].id
                    )
                })?;
                all_gts.extend_from_slice(&event.b_vcf.gt);
                all_gqs.push(sv_format.gq);
                all_pls.extend_from_slice(&sv_format.pl);
                all_ads.extend_from_slice(&sv_format.ad);
            } else {
                all_gts.extend_from_slice(&[
                    bcf::record::GenotypeAllele::UnphasedMissing,
                    bcf::record::GenotypeAllele::UnphasedMissing,
                ]);
                all_gqs.push(MISSING_INTEGER);
                all_pls.extend_from_slice(&[MISSING_INTEGER, MISSING_INTEGER, MISSING_INTEGER]);
                all_ads.extend_from_slice(&[MISSING_INTEGER, MISSING_INTEGER]);
            }
        }

        out_rec.push_genotypes(&all_gts)?;
        out_rec.push_format_integer(b"GQ", &all_gqs)?;
        out_rec.push_format_integer(b"PL", &all_pls)?;
        out_rec.push_format_integer(b"AD", &all_ads)?;
        writer.write(out_rec)?;
    }

    Ok(())
}

pub fn write_variants(
    variant_block_result: VariantBlockResult,
    sample_mapping: &SampleMapping,
    n_samples: usize,
    writer: &mut VcfWriter,
) -> Result<()> {
    log::debug!(
        "Write: Writing {} block for contig {}",
        variant_block_result.variant_type,
        variant_block_result.contig
    );

    let block_rid = if variant_block_result.variant_type == crate::core::svtype::SvType::BND {
        None
    } else {
        Some(
            writer
                .writer
                .header()
                .name2rid(variant_block_result.contig.as_bytes())
                .map_err(|e| {
                    crate::svx_error!(
                        "Output header missing contig {} while writing {} block: {e}",
                        variant_block_result.contig,
                        variant_block_result.variant_type
                    )
                })?,
        )
    };

    for mut group in variant_block_result.groups {
        if group.is_empty() {
            continue;
        }

        let mut keyed_group = group
            .into_iter()
            .map(|variant| {
                let merged_sample_pos = merged_sample_position_for_variant(
                    &variant,
                    sample_mapping,
                    "sorting write group",
                )?;
                Ok((merged_sample_pos, variant))
            })
            .collect::<Result<Vec<_>>>()?;
        keyed_group.sort_by(|(a_sample, a_variant), (b_sample, b_variant)| {
            a_sample
                .cmp(b_sample)
                .then_with(|| a_variant.id.cmp(&b_variant.id))
        });
        group = keyed_group
            .into_iter()
            .map(|(_, variant)| variant)
            .collect();
        // Representative variant index per sample in write order.
        let sample_representatives =
            create_sample_representatives(&group, sample_mapping, n_samples)?;

        let v0 = &group[0];
        if v0.svtype == SvType::BND {
            write_bnd_group(
                &group,
                &sample_representatives,
                &mut writer.writer,
                &mut writer.dummy_record,
            )?;
            continue;
        }

        let v0_vcf = v0
            .vcf
            .as_ref()
            .ok_or_else(|| crate::svx_error!("Variant {} is missing VCF write data", v0.id))?;
        let out_rec = &mut writer.dummy_record;
        out_rec.clear();

        out_rec.set_rid(block_rid);
        out_rec.set_pos(v0_vcf.pos);
        out_rec.set_id(v0.id.as_bytes())?;
        let output_svtype = if variant_block_result.variant_type == SvType::CNV {
            derive_cnv_output_svtype(&group)
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

        // FORMAT field specific work: this is (should be) sample dependent
        if variant_block_result.variant_type == SvType::CNV {
            let mut all_gts = Vec::new();
            let mut all_cns = Vec::new();
            let mut all_cnqs = Vec::new();

            for sample_variant in &sample_representatives {
                if let Some(group_idx) = sample_variant {
                    let variant = &group[*group_idx];
                    let variant_vcf = variant.vcf.as_ref().ok_or_else(|| {
                        crate::svx_error!("Variant {} is missing VCF write data", variant.id)
                    })?;
                    let cnv_format = variant_vcf.cnv_format().ok_or_else(|| {
                        crate::svx_error!(
                            "Variant {} is missing CNV sample FORMAT payload",
                            variant.id
                        )
                    })?;
                    all_gts.extend_from_slice(&variant_vcf.gt);
                    all_cns.push(cnv_format.cn);
                    all_cnqs.push(cnv_format.cnq);
                } else {
                    all_gts.extend_from_slice(&[
                        bcf::record::GenotypeAllele::UnphasedMissing,
                        bcf::record::GenotypeAllele::UnphasedMissing,
                    ]);
                    all_cns.push(MISSING_FLOAT);
                    all_cnqs.push(MISSING_FLOAT);
                }
            }

            out_rec.push_genotypes(&all_gts)?;
            out_rec.push_format_float(b"CN", &all_cns)?;
            out_rec.push_format_float(b"CNQ", &all_cnqs)?;
        } else {
            let mut all_gts = Vec::new();
            let mut all_gqs = Vec::new();
            let mut all_pls = Vec::new();
            let mut all_ads = Vec::new();

            for sample_variant in &sample_representatives {
                if let Some(group_idx) = sample_variant {
                    let variant = &group[*group_idx];
                    let variant_vcf = variant.vcf.as_ref().ok_or_else(|| {
                        crate::svx_error!("Variant {} is missing VCF write data", variant.id)
                    })?;
                    let sv_format = variant_vcf.sv_format().ok_or_else(|| {
                        crate::svx_error!(
                            "Variant {} is missing SV sample FORMAT payload",
                            variant.id
                        )
                    })?;

                    all_gts.extend_from_slice(&variant_vcf.gt);
                    all_gqs.push(sv_format.gq);
                    all_pls.extend_from_slice(&sv_format.pl);
                    all_ads.extend_from_slice(&sv_format.ad);
                } else {
                    // TODO: Deal with ploidy
                    all_gts.extend_from_slice(&[
                        bcf::record::GenotypeAllele::UnphasedMissing,
                        bcf::record::GenotypeAllele::UnphasedMissing,
                    ]);
                    all_gqs.push(MISSING_INTEGER);
                    all_pls.extend_from_slice(&[MISSING_INTEGER, MISSING_INTEGER, MISSING_INTEGER]);
                    all_ads.extend_from_slice(&[MISSING_INTEGER, MISSING_INTEGER]);
                }
            }

            out_rec.push_genotypes(&all_gts)?;
            out_rec.push_format_integer(b"GQ", &all_gqs)?;
            out_rec.push_format_integer(b"PL", &all_pls)?;
            out_rec.push_format_integer(b"AD", &all_ads)?;
        }

        writer.writer.write(out_rec)?;
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
    use crate::utils::util::init_logger;
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
    fn create_output_header_declares_supp_vec_as_variable_length_integer_list() {
        init_logger();

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
                "##INFO=<ID=SUPP_VEC,Number=.,Type=Integer,Description=\"Vector of supporting samples (0/1 per input sample)\">"
            ),
            "expected SUPP_VEC INFO header to be variable-length integer list, got header:\n{header}"
        );
    }

    #[test]
    fn create_output_header_declares_end_chr2_strands_and_svclaim_set() {
        init_logger();

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
    }

    fn common_sv_info_header_records(header: &mut rust_htslib::bcf::Header) {
        header.push_record(br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=END,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=START_VARIANCE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=START_AVG,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SVLEN_VARIANCE,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SVLEN_AVG,Number=1,Type=String,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP_VEC,Number=.,Type=Integer,Description="">"#);
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
        header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description="">"#);
        header.push_record(br#"##INFO=<ID=SUPP_VEC,Number=.,Type=Integer,Description="">"#);
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

    #[test]
    fn bnd_output_emits_paired_bnd_records_with_mateid() {
        init_logger();

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
        let bnd_a = VariantInternal::from_vcf_record(&r0, 0, "chr9", &args, &None).unwrap();
        let bnd_b = VariantInternal::from_vcf_record(&r1, 0, "chr1", &args, &None).unwrap();
        let event = VariantInternal::from_bnd_pair(bnd_a, bnd_b, &args).unwrap();

        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((0, 0), 0);
        sample_mapping.reverse_map.insert(0, (0, 0));

        let header = bnd_header(&["chr9", "chr1"], &["sample1"]);

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(
            &header,
            &Some(crate::io::vcf_writer::OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }),
            Some(&out_path),
            None,
        )
        .unwrap();

        let vbr = VariantBlockResult {
            groups: vec![vec![event]],
            contig: "chr1_chr9_TRA".to_string(),
            variant_type: SvType::BND,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, &mut writer).unwrap();
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
            assert!(!info.contains("END="), "record {id} should omit END");
            assert!(
                info.contains("STRANDS=--"),
                "record {id} should carry STRANDS"
            );
        }
    }

    #[test]
    fn vcf_info_float_fields_use_2_decimals_and_compact_zero_fraction() {
        init_logger();

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
            groups: vec![vec![v0, v1, v2]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 3, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let info = info_field(&non_header_lines[0]);
        assert!(
            info.contains("START_AVG=2.0"),
            "expected START_AVG to be formatted with compact trailing-zero fraction; got INFO={info}"
        );
        assert!(
            info.contains("START_VARIANCE=0.67"),
            "expected START_VARIANCE to be formatted to 2 decimals; got INFO={info}"
        );
        assert!(
            info.contains("SVLEN_AVG=2.0"),
            "expected SVLEN_AVG to be formatted with compact trailing-zero fraction; got INFO={info}"
        );
        assert!(
            info.contains("SVLEN_VARIANCE=0.67"),
            "expected SVLEN_VARIANCE to be formatted to 2 decimals; got INFO={info}"
        );
        assert!(
            !info.contains("START_AVG=2.00"),
            "expected START_AVG to avoid redundant trailing zeros; got INFO={info}"
        );
        assert!(
            !info.contains("SVLEN_AVG=2.00"),
            "expected SVLEN_AVG to avoid redundant trailing zeros; got INFO={info}"
        );
    }

    #[test]
    fn format_info_float_limits_to_2_decimals_and_compacts_zero_fraction() {
        assert_eq!(format_info_float(10000.0), "10000.0");
        assert_eq!(format_info_float(0.0), "0.0");
        assert_eq!(format_info_float(0.666_666), "0.67");
        assert_eq!(format_info_float(12.3), "12.3");
    }

    #[test]
    fn write_variants_returns_error_when_sort_sample_mapping_is_missing() {
        init_logger();

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

        let vbr = VariantBlockResult {
            groups: vec![vec![v0, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 2,
        };

        let err = write_variants(vbr, &sample_mapping, 1, &mut writer).unwrap_err();
        assert!(
            err.to_string().contains("while sorting write group"),
            "expected missing sample mapping error from sort path; got: {err}"
        );
    }

    #[test]
    fn write_variants_uses_block_contig_when_input_rid_points_elsewhere() {
        init_logger();

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
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, &mut writer).unwrap();
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
        init_logger();

        let sample_mapping = sample_mapping_one_sample();
        let header = single_sample_sv_header();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();
        assert_eq!(writer.dummy_record.id(), b".".to_vec());

        let mut variant =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        variant.vcf = Some(make_del_sv_vcf(0, [0, 1]));

        let vbr = VariantBlockResult {
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, &mut writer).unwrap();
        assert_eq!(writer.dummy_record.id(), b"v0".to_vec());
    }

    #[test]
    fn tr_contained_info_is_omitted_when_group_has_no_tr_annotation() {
        init_logger();

        let sample_mapping = sample_mapping_one_sample();
        let header = single_sample_sv_header();

        let out_path = make_temp_out("vcf");
        let mut writer = VcfWriter::new(&header, &None, Some(&out_path), None).unwrap();

        let mut variant =
            test_utils::from_parts(0, "v0".to_string(), SvType::DELETION, 0.0, 1.0).unwrap();
        variant.vcf = Some(make_del_sv_vcf(0, [0, 1]));

        let vbr = VariantBlockResult {
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, &mut writer).unwrap();
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
        init_logger();

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
            groups: vec![vec![variant]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 1,
        };

        write_variants(vbr, &sample_mapping, 1, &mut writer).unwrap();
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
        init_logger();

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
            groups: vec![vec![v0, v0_dup, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::DELETION,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 2, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let fields: Vec<&str> = non_header_lines[0].split('\t').collect();
        let info = fields[7];
        assert!(info.contains("SUPP=2"));
        assert!(info.contains("SUPP_CALLS=3"));
        assert!(info.contains("SUPP_VEC=1,1"));

        let sample2 = fields[10];
        let sample2_gq = find_sample_field(fields[8], sample2, "GQ");
        assert_eq!(
            sample2_gq, "30",
            "expected sample2 FORMAT to come from sample2 variant, not duplicate sample1 variant"
        );
    }

    #[test]
    fn write_variants_intrasample_cnv_uses_sample_support_and_sample_matched_format() {
        init_logger();

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
            groups: vec![vec![v0, v0_dup, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 2, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);

        let fields: Vec<&str> = non_header_lines[0].split('\t').collect();
        let info = fields[7];
        assert!(info.contains("SUPP=2"));
        assert!(info.contains("SUPP_CALLS=3"));
        assert!(info.contains("SUPP_VEC=1,1"));

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
    fn write_variants_cnv_emits_svclaim_set_for_consensus_cnv_type() {
        init_logger();

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
            groups: vec![vec![v0, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 1);
        let info = info_field(&non_header_lines[0]);
        assert!(!info_contains_field(info, "SVCLAIM=D"));
        assert!(info_contains_field(info, "SVCLAIM_SET=D"));
    }

    #[test]
    fn write_variants_cnv_emits_svclaim_set_when_claims_disagree() {
        init_logger();

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
            groups: vec![vec![v0, v1]],
            contig: "chr1".to_string(),
            variant_type: SvType::CNV,
            n: 2,
        };

        write_variants(vbr, &sample_mapping, 2, &mut writer).unwrap();
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
        init_logger();

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
                },
                inv_event_id: None,
            });
            v
        };

        let v0 = make_bnd_variant(0, "bnd0", 10, [0, 1]);
        let v0_dup = make_bnd_variant(0, "bnd0_dup", 20, [0, 1]);
        let v1 = make_bnd_variant(1, "bnd1", 30, [1, 1]);

        let vbr = VariantBlockResult {
            groups: vec![vec![v0, v0_dup, v1]],
            contig: "chr1_chr2_TRA".to_string(),
            variant_type: SvType::BND,
            n: 3,
        };

        write_variants(vbr, &sample_mapping, 2, &mut writer).unwrap();
        drop(writer);

        let non_header_lines = output_non_header_lines(&out_path);
        assert_eq!(non_header_lines.len(), 2);

        for line in &non_header_lines {
            let fields: Vec<&str> = line.split('\t').collect();
            let info = fields[7];
            assert!(info.contains("SUPP=2"));
            assert!(info.contains("SUPP_CALLS=3"));
            assert!(info.contains("SUPP_VEC=1,1"));

            let sample2 = fields[10];
            let sample2_gt = find_sample_field(fields[8], sample2, "GT");
            let sample2_gq = find_sample_field(fields[8], sample2, "GQ");
            assert_eq!(sample2_gt, "1/1");
            assert_eq!(sample2_gq, "30");
        }
    }
}
