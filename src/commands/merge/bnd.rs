use crate::{
    cli::MergeArgsInner,
    core::variant::{BndData, VariantInternal},
    io::positions_reader::PositionTuple,
    utils::util::Result,
};
use std::collections::{BTreeMap, HashMap, HashSet};

pub(crate) fn bucket_bnd_variants(
    all_bnd_variants: Vec<VariantInternal>,
) -> Result<BTreeMap<String, Vec<VariantInternal>>> {
    let mut by_graph: BTreeMap<String, Vec<VariantInternal>> = BTreeMap::new();
    for v in all_bnd_variants {
        let graph_id = v.bnd_event_graph_id()?;
        by_graph.entry(graph_id).or_default().push(v);
    }
    Ok(by_graph)
}

#[derive(Debug, Clone)]
pub(crate) struct BndMateSelectionScope {
    selected_contigs: HashSet<String>,
    selected_ranges_by_contig: Option<HashMap<String, Vec<(i64, i64)>>>,
}

impl BndMateSelectionScope {
    pub(crate) fn new(
        selected_contigs: &[String],
        target_positions: Option<&[PositionTuple]>,
    ) -> Self {
        let selected_contigs = selected_contigs.iter().cloned().collect();
        let selected_ranges_by_contig = target_positions.map(|positions| {
            let mut ranges_by_contig: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
            for pos_tuple in positions {
                let end = pos_tuple.end.map_or(pos_tuple.start + 1, |end| end + 1);
                ranges_by_contig
                    .entry(pos_tuple.contig.clone())
                    .or_default()
                    .push((pos_tuple.start, end));
            }
            ranges_by_contig
        });
        Self {
            selected_contigs,
            selected_ranges_by_contig,
        }
    }

    pub(crate) fn mate_is_in_scope(&self, mate_contig: &str, mate_pos0: i64) -> bool {
        if !self.selected_contigs.contains(mate_contig) {
            return false;
        }

        if let Some(ranges_by_contig) = &self.selected_ranges_by_contig {
            let Some(ranges) = ranges_by_contig.get(mate_contig) else {
                return false;
            };
            return ranges
                .iter()
                .any(|(start, stop)| mate_pos0 >= *start && mate_pos0 < *stop);
        }

        true
    }
}

fn should_warn_for_missing_bnd_mate(
    bnd: &BndData,
    selection_scope: Option<&BndMateSelectionScope>,
) -> bool {
    selection_scope.is_none_or(|scope| scope.mate_is_in_scope(&bnd.mate_contig, bnd.mate_pos))
}

pub(crate) fn pair_bnd_breakends(
    all_bnd_breakends: Vec<VariantInternal>,
    args: &MergeArgsInner,
) -> Result<Vec<VariantInternal>> {
    pair_bnd_breakends_with_scope(all_bnd_breakends, args, None)
}

pub(crate) fn pair_bnd_breakends_with_scope(
    all_bnd_breakends: Vec<VariantInternal>,
    args: &MergeArgsInner,
    selection_scope: Option<&BndMateSelectionScope>,
) -> Result<Vec<VariantInternal>> {
    let mut by_sample: BTreeMap<usize, Vec<VariantInternal>> = BTreeMap::new();
    for v in all_bnd_breakends {
        by_sample.entry(v.vcf_id).or_default().push(v);
    }

    let mut out = Vec::new();
    for (_sample_id, mut breakends) in by_sample {
        let mut by_id: BTreeMap<String, VariantInternal> = BTreeMap::new();
        for v in breakends.drain(..) {
            let id = v.id.clone();
            if by_id.insert(id.clone(), v).is_some() {
                return Err(crate::svx_error!("Duplicate BND ID within one VCF: {id}"));
            }
        }

        while let Some((id, v)) = by_id.pop_first() {
            let mate_id = v
                .bnd
                .as_ref()
                .ok_or_else(|| crate::svx_error!("BND variant {id} missing breakend data"))?
                .mate_id
                .clone();

            let mate = match by_id.remove(&mate_id) {
                Some(mate) => mate,
                None => {
                    let bnd = v.bnd.as_ref().expect("checked BND breakend data earlier");
                    let is_self_pointing =
                        bnd.contig == bnd.mate_contig && bnd.pos0 == bnd.mate_pos;
                    let has_alt_gt = v.vcf.as_ref().is_some_and(|vcf| {
                        vcf.gt.iter().any(|a| a.index().is_some_and(|idx| idx > 0))
                    });
                    let msg = if is_self_pointing {
                        "Dropping orphan self-pointing BND"
                    } else if has_alt_gt {
                        "Dropping orphan called BND"
                    } else {
                        "Dropping orphan uncalled BND"
                    };
                    if should_warn_for_missing_bnd_mate(bnd, selection_scope) {
                        log::warn!(
                            "{msg} {} (MATEID {} missing); possibly caused when one breakend of a BND pair gets removed by deduplication in sawfish ({}:{} -> {}:{})",
                            id,
                            mate_id,
                            bnd.contig,
                            bnd.pos0 + 1,
                            bnd.mate_contig,
                            bnd.mate_pos + 1
                        );
                    }
                    continue;
                }
            };

            let event = VariantInternal::from_bnd_pair(v, mate, args)?;
            out.push(event);
        }
    }

    Ok(out)
}
