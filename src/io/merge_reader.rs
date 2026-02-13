use std::collections::HashMap;

use crate::{
    cli::MergeArgsInner,
    core::{containers::interval_tree::IntervalTree, svtype::SvType, variant::VariantInternal},
    io::{
        bed_reader::{BedMap, TrId},
        vcf_reader::VcfReader,
    },
    utils::util::{Result, format_number_with_commas},
};
use rust_htslib::bcf::Read;

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

pub fn load_contig(
    contig: &str,
    args: &MergeArgsInner,
    selected_svtypes: &[SvType],
    tr_map: &Option<BedMap>,
    pos_map: &Option<HashMap<String, IntervalTree<i64, ()>>>,
    readers: &mut [VcfReader],
) -> Result<Vec<Vec<VariantInternal>>> {
    log::debug!("Reader: Processing contig: {}", contig);
    let contig_b = contig.as_bytes(); // TODO: Just keep it as bytes
    let svtype_mask = SvTypeMask::from_selected(selected_svtypes);

    // Each variant type has a corresponding vector: 0:INS, 1:DEL, 2:INV, 3:DUP, 4:CNV, 5:BND
    let mut variants_by_type: Vec<Vec<VariantInternal>> = vec![Vec::new(); SvType::BUCKET_COUNT];

    let mut tr_containment = match tr_map {
        Some(_) => 0,
        None => -1,
    };

    let bed_it: Option<&IntervalTree<u32, TrId>> = if let Some(tr_map) = tr_map {
        tr_map.interval_map.get(contig)
    } else {
        None
    };

    let pos_it = if let Some(pos_map) = pos_map {
        pos_map.get(contig)
    } else {
        None
    };

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
            // Tabix/CSI indices typically include only contigs that have at least one record.
            // It is common for headers to declare contigs (e.g. chrM) that don't appear in the
            // data, in which case fetch may fail even though name2rid succeeded.
            let msg = e.to_string();
            if msg.contains("error seeking to") {
                log::trace!(
                    "VCF[{}]: Contig {} not present in index (fetch failed: {}), skipping",
                    reader_index,
                    contig,
                    msg
                );
                continue;
            }
            return Err(crate::svx_error!(
                "VCF[{}]: Error fetching contig {}: {}",
                reader_index,
                contig,
                msg
            ));
        }

        let mut variant_count = 0;
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

            let v = VariantInternal::from_vcf_record(record, reader_index, contig, args, &bed_it)?;

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
            "{}",
            format_with_commas!(
                "VCF[{}]: Extracted {} variants (INS: {}, DEL: {}, INV: {}, DUP: {}, CNV: {}, BND: {})",
                reader_index,
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
            "{}",
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
}
