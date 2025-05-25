use crate::{
    cli::MergeArgs,
    core::{variant_block::VariantBlockResult, variant_internal::VariantInternal},
    io::{
        // readers::open_genome_reader,
        vcf_reader::{SampleMapping, VcfReaders},
        vcf_writer::VcfWriter,
    },
    utils::util::{Result, MISSING_INTEGER},
};
use anyhow::anyhow;
use rust_htslib::bcf;
use std::env;

pub fn create_output_header(
    vcf_readers: &VcfReaders,
    args: &MergeArgs,
) -> Result<(bcf::Header, SampleMapping)> {
    let mut out_header = bcf::Header::new();
    let sample_mapping = vcf_readers.merge_headers(&mut out_header, args.force_samples)?;

    // TODO: Add intrasample specific fields
    // TODO: Refine fields
    out_header.push_record(br#"##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">"#);
    out_header.push_record(
        br#"##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of supporting samples">"#,
    );
    out_header.push_record(br#"##INFO=<ID=IDLIST,Number=.,Type=String,Description="Variant IDs of variants merged to make this call (at most 1 per sample)">"#);
    out_header.push_record(br#"##INFO=<ID=START_VARIANCE,Number=1,Type=String,Description="Variance of start position for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=END_VARIANCE,Number=1,Type=String,Description="Variance of end position for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=START_AVG,Number=1,Type=String,Description="Average start position for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=END_AVG,Number=1,Type=String,Description="Average end position for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=LEN_AVG,Number=1,Type=String,Description="Average length for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=SVLEN_AVG,Number=1,Type=String,Description="Average SVLEN for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=SVLEN_VARIANCE,Number=1,Type=String,Description="Variance of SVLEN for variants merged into this one">"#);
    out_header.push_record(br#"##INFO=<ID=TR_CONTAINED,Number=1,Type=Integer,Description="At least one variant in this merged group was contained in a tandem repeat region">"#);

    if !args.no_version {
        add_version_info(&mut out_header);
    }
    Ok((out_header, sample_mapping))
}

fn add_version_info(out_header: &mut bcf::Header) {
    let version_line = format!(
        "##{}Version={}",
        env!("CARGO_PKG_NAME"),
        *crate::cli::FULL_VERSION
    );
    out_header.push_record(version_line.as_bytes());

    let command_line = env::args().collect::<Vec<String>>().join(" ");
    let command_line = format!("##{}Command={}", env!("CARGO_PKG_NAME"), command_line);
    out_header.push_record(command_line.as_bytes());
}

fn create_sample_mask(
    group: &[VariantInternal],
    sample_mapping: &SampleMapping,
    n_samples: usize,
) -> Result<Vec<bool>> {
    let mut sample_mask = vec![false; n_samples];

    for variant in group {
        // TODO: Currently assumes single sample per VCF (sample_id = 0)
        if let Some(&merged_pos) = sample_mapping.index_map.get(&(variant.vcf_id, 0)) {
            sample_mask[merged_pos] = true;
        } else {
            return Err(anyhow!(
                "Failed to find sample mapping for variant {} (vcf: {}, sample: {})",
                variant.id,
                variant.vcf_id,
                variant.sample_id
            ));
        }
    }
    Ok(sample_mask)
}

pub fn write_variants(
    variant_block_result: VariantBlockResult,
    sample_mapping: &SampleMapping,
    n_samples: usize,
    writer: &mut VcfWriter,
) -> Result<()> {
    log::info!(
        "Write: Writing {} block for contig {}",
        variant_block_result.variant_type,
        variant_block_result.contig
    );

    for mut group in variant_block_result.groups {
        if group.is_empty() {
            continue;
        }

        // A group is created in arbitrary order, sort it to the write order
        group.sort_by_key(|variant| sample_mapping.index_map.get(&(variant.vcf_id, 0)).unwrap());
        // Mask to indicate which samples should have null values written and from which we want to pull values
        let sample_mask = create_sample_mask(&group, sample_mapping, n_samples)?;

        let v0 = &group[0];
        let mut out_rec = writer.writer.empty_record();

        out_rec.set_rid(v0.rid);
        out_rec.set_pos(v0.pos);
        out_rec.set_id(v0.id.as_bytes())?;
        let alleles: Vec<&[u8]> = v0.alleles.iter().map(|v| v.as_slice()).collect();
        out_rec.set_alleles(&alleles)?;

        // INFO field specific work: this is aggregate level (dependent)
        let supp_vec: String = sample_mask
            .iter()
            .map(|&present| if present { "1" } else { "0" })
            .collect::<Vec<_>>()
            .join(",");

        out_rec.push_info_string(b"SVTYPE", &[v0.svtype.to_string().as_bytes()])?;
        out_rec.push_info_integer(b"SVLEN", &[v0.svlen as i32])?;

        let n = group.len() as f32;
        let mean_start = group.iter().map(|v| v.start + 1.0).sum::<f32>() / n;
        let start_variance = group
            .iter()
            .map(|v| {
                let diff = (v.start + 1.0) - mean_start;
                diff * diff
            })
            .sum::<f32>()
            / n;
        out_rec.push_info_string(b"START_VARIANCE", &[start_variance.to_string().as_bytes()])?;
        out_rec.push_info_string(b"START_AVG", &[mean_start.to_string().as_bytes()])?;

        // Calculate SVLEN average and variance
        let mean_svlen = group.iter().map(|v| v.svlen).sum::<f32>() / n;
        let svlen_variance = group
            .iter()
            .map(|v| {
                let diff = v.svlen - mean_svlen;
                diff * diff
            })
            .sum::<f32>()
            / n;

        out_rec.push_info_string(b"SVLEN_VARIANCE", &[svlen_variance.to_string().as_bytes()])?;
        out_rec.push_info_string(b"SVLEN_AVG", &[mean_svlen.to_string().as_bytes()])?;

        out_rec.push_info_integer(b"SUPP", &[group.len() as i32])?;
        out_rec.push_info_string(b"SUPP_VEC", &[supp_vec.as_bytes()])?;

        let id_list: Vec<&[u8]> = group.iter().map(|v| v.id.as_bytes()).collect();
        out_rec.push_info_string(b"IDLIST", &id_list)?;

        let tr_contained_value = match group.iter().find(|v| v.trid.is_some()) {
            Some(v) => i32::from(v.trid.is_some()),
            None => -1,
        };
        out_rec.push_info_integer(b"TR_CONTAINED", &[tr_contained_value])?;

        // TODO: Clean this up
        let mut current_group_idx = 0;
        let mut all_gts = Vec::new();
        let mut all_gqs = Vec::new();
        let mut all_pls = Vec::new();
        let mut all_ads = Vec::new();

        // FORMAT field specific work: this is (should be) sample dependent
        for is_present in &sample_mask {
            if *is_present {
                let variant = &group[current_group_idx];

                all_gts.extend_from_slice(&variant.gt);
                all_gqs.push(variant.gq);
                all_pls.extend_from_slice(&variant.pl);
                all_ads.extend_from_slice(&variant.ad);
                current_group_idx += 1;
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

        writer.writer.write(&out_rec)?;
    }

    Ok(())
}
