use crate::{
    cli::MergeArgsInner, core::variant::VariantInternal, io::positions_reader::PositionTuple,
    utils::util::Result,
};
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct BndMateSelectionScope {
    selected_contigs: HashSet<String>,
    selected_ranges_by_contig: Option<HashMap<String, Vec<(i64, i64)>>>,
}

impl BndMateSelectionScope {
    pub fn new(selected_contigs: &[String], target_positions: Option<&[PositionTuple]>) -> Self {
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

    pub fn mate_is_in_scope(&self, mate_contig: &str, mate_pos0: i64) -> bool {
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

pub fn bucket_bnd_variants(
    all_bnd_variants: Vec<VariantInternal>,
) -> Result<BTreeMap<String, Vec<VariantInternal>>> {
    let mut by_graph: BTreeMap<String, Vec<VariantInternal>> = BTreeMap::new();
    for v in all_bnd_variants {
        let graph_id = v.bnd_event_graph_id()?;
        by_graph.entry(graph_id).or_default().push(v);
    }
    Ok(by_graph)
}

pub fn pair_bnd_breakends(
    all_bnd_breakends: Vec<VariantInternal>,
    args: &MergeArgsInner,
    selection_scope: Option<&BndMateSelectionScope>,
) -> Result<Vec<VariantInternal>> {
    let mut by_sample: BTreeMap<usize, Vec<VariantInternal>> = BTreeMap::new();
    for v in all_bnd_breakends {
        by_sample.entry(v.vcf_id).or_default().push(v);
    }

    let mut out = Vec::new();
    for (_sample_id, breakends) in by_sample {
        let mut by_id: BTreeMap<String, Vec<VariantInternal>> = BTreeMap::new();
        for v in breakends {
            let id = v.id.clone();
            by_id.entry(id).or_default().push(v);
        }

        while let Some(id) = by_id.first_key_value().map(|(id, _)| id.clone()) {
            let mut variants_with_id = by_id
                .remove(&id)
                .ok_or_else(|| crate::svx_error!("Missing BND bucket for ID {}", id))?;
            let v = variants_with_id
                .pop()
                .ok_or_else(|| crate::svx_error!("Missing BND variant for ID {}", id))?;
            if !variants_with_id.is_empty() {
                by_id.insert(id.clone(), variants_with_id);
            }

            let mate_id = v
                .bnd
                .as_ref()
                .ok_or_else(|| crate::svx_error!("BND variant {id} missing breakend data"))?
                .mate_id
                .clone();

            let mut paired_event = None;
            let mut pairing_error = None;
            let candidate_attempts = by_id
                .get(&mate_id)
                .map(|candidates| {
                    candidates
                        .iter()
                        .enumerate()
                        .filter(|(_, candidate)| {
                            candidate
                                .bnd
                                .as_ref()
                                .is_some_and(|breakend| breakend.mate_id == id)
                        })
                        .map(|(candidate_idx, candidate)| (candidate_idx, candidate.clone()))
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();

            let mut selected_mate_idx = None;
            for (candidate_idx, candidate) in &candidate_attempts {
                match VariantInternal::from_bnd_pair(v.clone(), candidate.clone(), args) {
                    Ok(event) => {
                        if paired_event.is_some() {
                            return Err(crate::svx_error!(
                                "Ambiguous BND pairing for {} (MATEID {}): multiple reciprocal candidates can be paired",
                                id,
                                mate_id
                            ));
                        }
                        paired_event = Some(event);
                        selected_mate_idx = Some(*candidate_idx);
                    }
                    Err(error) => {
                        pairing_error = Some(error);
                    }
                }
            }
            if let Some(mate_idx) = selected_mate_idx {
                if let Some(candidates) = by_id.get_mut(&mate_id) {
                    candidates.swap_remove(mate_idx);
                    if candidates.is_empty() {
                        by_id.remove(&mate_id);
                    }
                }
            } else if !candidate_attempts.is_empty() {
                if let Some(error) = pairing_error {
                    return Err(error);
                }
            }

            if let Some(event) = paired_event {
                out.push(event);
                continue;
            }

            let bnd = v.bnd.as_ref().expect("checked BND breakend data earlier");
            let is_self_pointing = bnd.contig == bnd.mate_contig && bnd.pos0 == bnd.mate_pos;
            let has_alt_gt = v.support_mask.iter().any(|mask_word| *mask_word != 0);
            let mate_in_scope = selection_scope
                .is_none_or(|scope| scope.mate_is_in_scope(&bnd.mate_contig, bnd.mate_pos));
            let msg = if is_self_pointing {
                "Dropping orphan self-pointing BND"
            } else if has_alt_gt {
                "Dropping orphan called BND"
            } else {
                "Dropping orphan uncalled BND"
            };
            if mate_in_scope {
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
        }
    }

    Ok(out)
}
