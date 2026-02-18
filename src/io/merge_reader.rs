use std::collections::HashMap;

use crate::{
    cli::MergeArgsInner,
    core::{containers::interval_tree::IntervalTree, svtype::SvType, variant::VariantInternal},
    io::{
        bed_reader::BedMap,
        vcf_reader::{SampleMapping, VcfReader},
    },
    utils::util::{Result, format_number_with_commas},
};
use rust_htslib::{bcf::Read, errors::Error as HtsError};

macro_rules! format_with_commas {
    ($fmt:expr, $($args:expr),*) => {
        format!(
            $fmt,
            $(format_number_with_commas($args)),*
        )
    };
}

#[derive(Debug, Clone, Copy)]
struct SvTypeMask {
    selected: [bool; SvType::BUCKET_COUNT],
}

impl SvTypeMask {
    fn from_selected(svtypes: &[SvType]) -> Self {
        let mut selected = [false; SvType::BUCKET_COUNT];
        for &svtype in svtypes {
            selected[svtype.bucket_index()] = true;
        }
        Self { selected }
    }

    fn contains(self, svtype: SvType) -> bool {
        self.selected[svtype.bucket_index()]
    }
}

#[inline]
fn position_is_in_target_positions(interval_tree: &IntervalTree<i64, ()>, position: i64) -> bool {
    let end = position.saturating_add(1);
    interval_tree.any_containing(position, end)
}

fn build_local_to_merged_sample_map(
    sample_mapping: &SampleMapping,
    reader_index: usize,
    sample_count: usize,
) -> Result<Vec<usize>> {
    let mut local_to_merged = Vec::with_capacity(sample_count);
    for local_sample_idx in 0..sample_count {
        let merged_sample_idx = sample_mapping
            .index_map
            .get(&(reader_index, local_sample_idx))
            .copied()
            .ok_or_else(|| {
                crate::svx_error!(
                    "Missing sample mapping for VCF {} sample {}",
                    reader_index,
                    local_sample_idx
                )
            })?;
        local_to_merged.push(merged_sample_idx);
    }

    Ok(local_to_merged)
}

pub fn load_contig(
    contig: &str,
    args: &MergeArgsInner,
    selected_svtypes: &[SvType],
    tr_map: &Option<BedMap>,
    pos_map: &Option<HashMap<String, IntervalTree<i64, ()>>>,
    sample_mapping: &SampleMapping,
    readers: &mut [VcfReader],
) -> Result<Vec<Vec<VariantInternal>>> {
    log::debug!("Reader: Processing contig: {}", contig);
    let svtype_mask = SvTypeMask::from_selected(selected_svtypes);
    let mut variants_by_type: Vec<Vec<VariantInternal>> = vec![Vec::new(); SvType::BUCKET_COUNT];

    let mut tr_containment = match tr_map {
        Some(_) => 0,
        None => -1,
    };

    let bed_it = if let Some(tr_map) = tr_map {
        tr_map.interval_map.get(contig)
    } else {
        None
    };

    let pos_it = if let Some(pos_map) = pos_map {
        pos_map.get(contig)
    } else {
        None
    };

    let contig_b = contig.as_bytes();
    for (reader_index, reader) in readers.iter_mut().enumerate() {
        let rid = match reader.reader.header().name2rid(contig_b) {
            Ok(id) => id,
            Err(_) => {
                log::trace!(
                    "VCF[{}]: Contig {} not found, skipping",
                    reader_index,
                    contig,
                );
                continue;
            }
        };

        if let Err(e) = reader.reader.fetch(rid as u32, 0, None) {
            if matches!(e, HtsError::GenomicSeek { .. }) {
                log::debug!(
                    "VCF[{}]: Contig {} declared in header but absent from index; skipping ({})",
                    reader_index,
                    contig,
                    e
                );
                continue;
            }
            return Err(crate::svx_error!(
                "VCF[{}]: Error fetching contig {}: {}",
                reader_index,
                contig,
                e
            ));
        }

        let mut variant_count = 0;
        let format_cache = VariantInternal::build_header_format_cache(reader.reader.header())?;
        let local_to_merged = build_local_to_merged_sample_map(
            sample_mapping,
            reader_index,
            reader.sample_n as usize,
        )?;
        while reader.advance()? {
            let record = &reader.current_record;

            if let Some(it) = &pos_it {
                if !position_is_in_target_positions(it, record.pos()) {
                    continue;
                }
            }
            log::trace!("VCF[{}]: {}:{}", reader_index, contig, record.pos());

            let bucket_svtype = SvType::classify_vcf_record(record)?;
            if !svtype_mask.contains(bucket_svtype) {
                continue;
            }

            let mut v = VariantInternal::from_vcf_record(
                record,
                reader_index,
                contig,
                args,
                &bed_it,
                &format_cache,
                Some(local_to_merged.as_slice()),
            )?;
            v.remap_support_mask_to_output_with_local_map(
                &local_to_merged,
                sample_mapping.index_map.len(),
            )?;

            if !args.keep_monomorphic && !v.support_mask.iter().any(|mask_word| *mask_word != 0) {
                continue;
            }

            if args.filter_tr_contained && v.trid.is_some() {
                continue;
            }

            let type_index = bucket_svtype.bucket_index();

            if v.trid.is_some() {
                tr_containment += 1;
            }

            variants_by_type[type_index].push(v);
            variant_count += 1;
        }

        log::debug!(
            "VCF[{}] contig={}: {}",
            reader_index,
            contig,
            format_with_commas!(
                "Extracted {} variants (INS: {}, DEL: {}, INV: {}, DUP: {}, CNV: {}, BND: {})",
                variant_count,
                variants_by_type[0].len(),
                variants_by_type[1].len(),
                variants_by_type[2].len(),
                variants_by_type[3].len(),
                variants_by_type[4].len(),
                variants_by_type[5].len()
            )
        );
    }

    let total_variants: usize = variants_by_type.iter().map(|v| v.len()).sum();
    if total_variants > 0 {
        log::debug!(
            "Contig {}: {}",
            contig,
            format_with_commas!(
                "Extracted total of {} variants (INS: {}, DEL: {}, INV: {}, DUP: {}, CNV: {}, BND: {}), TR containment: {}",
                total_variants,
                variants_by_type[0].len(),
                variants_by_type[1].len(),
                variants_by_type[2].len(),
                variants_by_type[3].len(),
                variants_by_type[4].len(),
                variants_by_type[5].len(),
                tr_containment
            )
        );
    }

    Ok(variants_by_type)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::containers::interval_tree::Interval;

    #[test]
    fn svtype_mask_contains_only_selected_types() {
        let mask = SvTypeMask::from_selected(&[SvType::DELETION, SvType::BND]);
        assert!(mask.contains(SvType::DELETION));
        assert!(mask.contains(SvType::BND));
        assert!(!mask.contains(SvType::INSERTION));
        assert!(!mask.contains(SvType::INVERSION));
        assert!(!mask.contains(SvType::DUPLICATION));
        assert!(!mask.contains(SvType::CNV));
    }

    #[test]
    fn position_selection_uses_half_open_upper_bound() {
        let interval_tree = IntervalTree::new(vec![Interval::new(100, 105, ())]);
        assert!(position_is_in_target_positions(&interval_tree, 100));
        assert!(position_is_in_target_positions(&interval_tree, 104));
        assert!(!position_is_in_target_positions(&interval_tree, 105));
    }

    #[test]
    fn build_local_to_merged_sample_map_returns_dense_mapping() {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((2, 0), 5);
        sample_mapping.index_map.insert((2, 1), 3);
        sample_mapping.index_map.insert((2, 2), 8);
        let local_to_merged = build_local_to_merged_sample_map(&sample_mapping, 2, 3).unwrap();
        assert_eq!(local_to_merged, vec![5, 3, 8]);
    }

    #[test]
    fn build_local_to_merged_sample_map_errors_on_missing_slot() {
        let mut sample_mapping = SampleMapping::new();
        sample_mapping.index_map.insert((1, 0), 4);

        let err = build_local_to_merged_sample_map(&sample_mapping, 1, 2).unwrap_err();
        assert!(
            err.to_string()
                .contains("Missing sample mapping for VCF 1 sample 1"),
            "unexpected error: {err}"
        );
    }
}
