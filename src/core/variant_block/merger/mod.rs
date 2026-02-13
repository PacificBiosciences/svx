mod constraints;
mod neighbor_search;

use self::{constraints::CentroidState, neighbor_search::NeighborSearchState};
use super::{edge::Edge, label::VariantLabel};
use crate::{
    cli::{MergeArgsInner, MergeConstraint},
    core::aligner::{AlignmentScope, MemoryModel},
    core::{
        aligner::WFAligner,
        containers::{
            forest::Forest,
            kd_tree::{KnnScratch, VariantKdTree},
        },
        variant::VariantInternal,
    },
};
use std::collections::{BinaryHeap, HashMap};

pub(super) struct VariantMerger<'a> {
    variants: &'a [VariantInternal],
    tree: &'a VariantKdTree,
    forest: &'a mut Forest,
    args: &'a MergeArgsInner,
    merge_constraint: MergeConstraint,
    ed_aligner: WFAligner,
    id_ranks: Vec<u32>,
    knn_scratch: KnnScratch,
    pub(super) search_state: Vec<NeighborSearchState>,
    pub(super) heap: BinaryHeap<Edge>,
    members: Option<Vec<Vec<usize>>>,
    centroid_state: Option<CentroidState>,
    comparison_cache: HashMap<(usize, usize), bool>, // Cache for comparison results (true=pass, false=fail)
}

impl<'a> VariantMerger<'a> {
    pub(super) fn new(
        variants: &'a [VariantInternal],
        tree: &'a VariantKdTree,
        forest: &'a mut Forest,
        args: &'a MergeArgsInner,
    ) -> Self {
        let n = variants.len();
        let ed_aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
            .edit()
            .build();

        let merge_constraint = args.merge_constraint;
        let members = matches!(
            merge_constraint,
            MergeConstraint::Clique | MergeConstraint::Centroid
        )
        .then(|| (0..n).map(|i| vec![i]).collect());
        let centroid_state =
            (merge_constraint == MergeConstraint::Centroid).then(|| CentroidState::new(variants));

        let id_ranks = Self::compute_id_ranks(variants);

        let mut knn_scratch = KnnScratch::new(args.knn_search_k);
        let mut search_state: Vec<NeighborSearchState> = Vec::with_capacity(n);
        for variant in variants {
            let mut initial_neighbors: Vec<usize> = Vec::new();
            Self::fill_neighbor_indices(
                tree,
                variant,
                args.knn_search_k,
                &mut knn_scratch,
                &mut initial_neighbors,
            );
            search_state.push(NeighborSearchState {
                neighbors: initial_neighbors,
                next_neighbor_idx: 0,
            });
        }
        let mut merger = Self {
            variants,
            tree,
            forest,
            args,
            merge_constraint,
            ed_aligner,
            id_ranks,
            knn_scratch,
            search_state,
            heap: BinaryHeap::new(),
            members,
            centroid_state,
            comparison_cache: HashMap::new(),
        };
        merger.initialize_heap();
        merger
    }

    fn compute_id_ranks(variants: &[VariantInternal]) -> Vec<u32> {
        let mut by_id: Vec<(usize, &str)> = variants
            .iter()
            .enumerate()
            .map(|(idx, v)| (idx, v.id.as_str()))
            .collect();
        by_id.sort_by(|a, b| a.1.cmp(b.1).then_with(|| a.0.cmp(&b.0)));

        let mut ranks = vec![0u32; variants.len()];
        for (rank, (idx, _)) in by_id.into_iter().enumerate() {
            ranks[idx] = rank as u32;
        }
        ranks
    }

    fn initialize_heap(&mut self) {
        for i in 0..self.variants.len() {
            let v = &self.variants[i];
            self.heap.push(Edge::new(
                i,
                i,
                0.0,
                v.info_hash,
                v.info_hash,
                self.id_ranks[i],
                self.id_ranks[i],
            ));
        }
    }

    pub(super) fn execute(&mut self) {
        let n = self.variants.len();
        if n <= 1 {
            return;
        }

        while let Some(edge) = self.heap.pop() {
            let unioned = self.try_union_with_constraint(edge.from, edge.to);
            {
                let from_id = self
                    .variants
                    .get(edge.from)
                    .map(|v| v.id.as_str())
                    .unwrap_or("<invalid-variant-index>");
                let to_id = self
                    .variants
                    .get(edge.to)
                    .map(|v| v.id.as_str())
                    .unwrap_or("<invalid-variant-index>");
                let from_label = VariantLabel::new(edge.from, from_id);
                let to_label = VariantLabel::new(edge.to, to_id);

                if unioned {
                    log::trace!("United {} and {}", from_label, to_label);
                } else {
                    log::trace!(
                        "Cannot union {} and {} (already connected or constraint violation)",
                        from_label,
                        to_label
                    );
                }
            }

            if let Some(next_edge) = self.find_next_edge_for_variant(edge.from) {
                {
                    let from_id = self
                        .variants
                        .get(next_edge.from)
                        .map(|v| v.id.as_str())
                        .unwrap_or("<invalid-variant-index>");
                    let to_id = self
                        .variants
                        .get(next_edge.to)
                        .map(|v| v.id.as_str())
                        .unwrap_or("<invalid-variant-index>");
                    let from_label = VariantLabel::new(next_edge.from, from_id);
                    let to_label = VariantLabel::new(next_edge.to, to_id);
                    log::trace!(
                        "Pushing next edge: {} -> {} (dist: {})",
                        from_label,
                        to_label,
                        next_edge.dist
                    );
                }
                self.heap.push(next_edge);
            } else {
                let from_id = self
                    .variants
                    .get(edge.from)
                    .map(|v| v.id.as_str())
                    .unwrap_or("<invalid-variant-index>");
                log::trace!(
                    "No more valid edges found for {}",
                    VariantLabel::new(edge.from, from_id)
                );
            }
        }
        log::trace!("Merging heap empty, process finished.");
    }
}
