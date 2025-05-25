use std::collections::HashMap;

use crate::{
    cli::MergeArgsInner,
    core::{
        containers::interval_tree::IntervalTree, svtype::SvType, variant_internal::VariantInternal,
    },
    io::{
        bed_reader::{BedMap, TrId},
        vcf_reader::VcfReader,
    },
    utils::util::{format_number_with_commas, Result},
};

macro_rules! format_with_commas {
    ($fmt:expr, $($args:expr),*) => {
        format!(
            $fmt,
            $(format_number_with_commas($args)),*
        )
    };
}

pub fn load_contig(
    contig: &str,
    args: &MergeArgsInner,
    tr_map: &Option<BedMap>,
    pos_map: &Option<HashMap<String, IntervalTree<i64, ()>>>,
    readers: &mut [VcfReader],
) -> Result<Vec<Vec<VariantInternal>>> {
    log::info!("Reader: Processing contig: {}", contig);
    let contig_b = contig.as_bytes(); // TODO: Just keep it as bytes

    // Each variant type has a corresponding vector: 0:INS, 1:DEL, 2:INV, 3:DUP, 4:BND
    let mut variants_by_type: Vec<Vec<VariantInternal>> = vec![Vec::new(); 5];

    let mut tr_containment = match tr_map {
        Some(_) => 0,
        None => -1,
    };

    let bed_it: Option<&IntervalTree<u32, TrId>> = if let Some(ref tr_map) = tr_map {
        tr_map.interval_map.get(contig)
    } else {
        None
    };

    let pos_it = if let Some(ref pos_map) = pos_map {
        pos_map.get(contig)
    } else {
        None
    };

    for (reader_index, reader) in readers.iter_mut().enumerate() {
        let rid = match reader.header.name2rid(contig_b) {
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
            log::trace!(
                "VCF[{}]: Error fetching contig {}: {}. Skipping.",
                reader_index,
                contig,
                e
            );
            continue;
        }

        let mut variant_count = 0;
        while reader.advance() {
            let record = &reader.current_record;

            if let Some(it) = &pos_it {
                let overlaps = it.find_containing(record.pos(), record.pos());
                if overlaps.is_empty() {
                    continue;
                }
            }
            log::trace!("VCF[{}]: {}:{}", reader_index, contig, record.pos());

            let v = VariantInternal::from_vcf_record(record, reader_index, args, &bed_it)?;

            // Skip cns for now, they should be put into their own bucket
            if v.is_cn {
                continue;
            }

            if args.filter_tr_contained && v.trid.is_some() {
                continue;
            }

            let type_index = match v.svtype {
                SvType::INSERTION => 0,
                SvType::DELETION => 1,
                SvType::INVERSION => 2,
                SvType::DUPLICATION => 3,
                SvType::BND => 4,
            };

            if v.trid.is_some() {
                tr_containment += 1;
            }

            variants_by_type[type_index].push(v);
            variant_count += 1;
        }

        log::debug!(
            "{}",
            format_with_commas!(
                "VCF[{}]: Extracted {} variants (INS: {}, DEL: {}, INV: {}, DUP: {}, BND: {})",
                reader_index,
                variant_count,
                variants_by_type[0].len(),
                variants_by_type[1].len(),
                variants_by_type[2].len(),
                variants_by_type[3].len(),
                variants_by_type[4].len()
            )
        );
    }
    let total_variants: usize = variants_by_type.iter().map(|v| v.len()).sum();
    if total_variants > 0 {
        log::debug!(
                "{}", format_with_commas!(
                    "Extracted total of {} variants (INS: {}, DEL: {}, INV: {}, DUP: {}, BND: {}), TR containment: {}",
                    total_variants,
                    variants_by_type[0].len(),
                    variants_by_type[1].len(),
                    variants_by_type[2].len(),
                    variants_by_type[3].len(),
                    variants_by_type[4].len(),
                    tr_containment
                )
            );
    }

    Ok(variants_by_type)
}
