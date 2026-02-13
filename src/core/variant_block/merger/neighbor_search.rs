use super::VariantMerger;
use crate::{
    DISTANCE_OFFSET,
    core::{
        containers::kd_tree::{KnnScratch, VariantKdTree},
        variant::VariantInternal,
        variant_block::merger::{Edge, VariantLabel},
    },
};

pub(in crate::core::variant_block) struct NeighborSearchState {
    pub(in crate::core::variant_block) neighbors: Vec<usize>,
    pub(in crate::core::variant_block) next_neighbor_idx: usize,
}

impl VariantMerger<'_> {
    pub(super) fn fill_neighbor_indices(
        tree: &VariantKdTree,
        variant: &VariantInternal,
        k: usize,
        scratch: &mut KnnScratch,
        out: &mut Vec<usize>,
    ) {
        out.clear();
        out.reserve(k);
        out.extend(tree.knn_into(variant, k, scratch).iter().map(|n| n.item));
    }

    // Finds the next valid edge originating from `from_idx` based on its search state
    pub(in crate::core::variant_block) fn find_next_edge_for_variant(
        &mut self,
        from_idx: usize,
    ) -> Option<Edge> {
        let from_variant = &self.variants[from_idx];
        let from_label = VariantLabel::new(from_idx, from_variant.id.as_str());
        let require_mutual_distance = self.args.require_mutual_distance;
        let tr_max_dist = self.args.tr_max_dist;
        let state = &mut self.search_state[from_idx];
        loop {
            if state.next_neighbor_idx >= state.neighbors.len() {
                let current_k = state.neighbors.len();
                if current_k == 0 || current_k >= self.variants.len() {
                    log::trace!(
                        "Stopping neighbor search for {}: no more potential neighbors (current_k={}, total_variants={}).",
                        from_label,
                        current_k,
                        self.variants.len()
                    );
                    return None;
                }
                let next_k = (current_k * 2).min(self.variants.len());
                log::trace!(
                    "Finding more neighbors for {}: current_k={}, next_k={}",
                    from_label,
                    current_k,
                    next_k
                );
                Self::fill_neighbor_indices(
                    self.tree,
                    from_variant,
                    next_k,
                    &mut self.knn_scratch,
                    &mut state.neighbors,
                );

                // Check again after potentially fetching more neighbors
                if state.next_neighbor_idx >= state.neighbors.len() {
                    log::trace!(
                        "Stopping neighbor search for {}: doubling k ({}) did not yield new neighbors beyond index {}.",
                        from_label,
                        next_k,
                        state.next_neighbor_idx
                    );
                    return None;
                }
            }

            let to_idx = state.neighbors[state.next_neighbor_idx];
            state.next_neighbor_idx += 1; // Consume this neighbor *before* any continue or return statements for this iteration

            let candidate_to = &self.variants[to_idx];
            let to_label = VariantLabel::new(to_idx, candidate_to.id.as_str());
            if to_idx == from_idx {
                continue;
            }

            // Start validation checks
            let distance = from_variant.distance(candidate_to);

            // Condition: Stop searching further for `from_idx` if the candidate is too far (exceeds own max_dist threshold)
            // This check is NOT symmetric and determines if we stop searching *for from_idx*, so it comes before caching.
            let stop_threshold = if from_variant.trid.is_some() {
                from_variant.max_dist.max(f64::from(tr_max_dist))
            } else {
                from_variant.max_dist
            };
            if distance > stop_threshold + DISTANCE_OFFSET {
                log::trace!(
                    "Stopping search for {}: candidate {} is too far ({} > {} + offset).",
                    from_label,
                    to_label,
                    distance,
                    stop_threshold
                );
                return None;
            }

            // Cache check
            let cache_key = (from_idx.min(to_idx), from_idx.max(to_idx));
            if let Some(&cached_result) = self.comparison_cache.get(&cache_key) {
                if cached_result {
                    log::trace!(
                        "Found valid edge {} -> {} (dist: {}, cached pass).",
                        from_label,
                        to_label,
                        distance
                    );
                    return Some(Edge::new(
                        from_idx,
                        to_idx,
                        distance,
                        from_variant.info_hash,
                        candidate_to.info_hash,
                        self.id_ranks[from_idx],
                        self.id_ranks[to_idx],
                    ));
                }

                log::trace!(
                    "Skipping edge {} -> {}: failed cached checks.",
                    from_label,
                    to_label
                );
                continue;
            }

            let v1_is_tr = from_variant.trid.is_some();
            let v2_is_tr = candidate_to.trid.is_some();
            let mut apply_tr_rules = false;
            if v1_is_tr && v2_is_tr {
                if from_variant.trid.as_ref().unwrap().id != candidate_to.trid.as_ref().unwrap().id
                {
                    // Both are TRs but different TR IDs, do not merge.
                    log::trace!(
                        "Skipping edge {} -> {}: both TR-contained but different TR IDs ('{:?}' vs '{:?}'). Caching fail.",
                        from_label,
                        to_label,
                        from_variant.trid.as_ref().unwrap().id,
                        candidate_to.trid.as_ref().unwrap().id
                    );
                    self.comparison_cache.insert(cache_key, false);
                    continue;
                }
                apply_tr_rules = true;
            }

            // If not in the cache, perform symmetric checks and cache the result
            let effective_max_dist_for_pair = if apply_tr_rules {
                f64::from(tr_max_dist)
            } else if require_mutual_distance {
                from_variant.max_dist.min(candidate_to.max_dist)
            } else {
                from_variant.max_dist.max(candidate_to.max_dist)
            };

            if distance > effective_max_dist_for_pair + DISTANCE_OFFSET {
                log::trace!(
                    "Skipping edge {} -> {}: exceeds effective max distance for pair ({} > {} + offset). TR lenient: {}. Caching fail.",
                    from_label,
                    to_label,
                    distance,
                    effective_max_dist_for_pair,
                    apply_tr_rules
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            // When intrasample merging is disabled, we still schedule the edge and defer
            // rejection to the union-find sample constraint (`Forest::valid_edge`).
            if !self.args.allow_intrasample && from_variant.sample_id == candidate_to.sample_id {
                log::trace!(
                    "Deferring intrasample rejection for edge {} -> {} (sample_id: {}); sample constraint is enforced during union.",
                    from_label,
                    to_label,
                    from_variant.sample_id
                );
                return Some(Edge::new(
                    from_idx,
                    to_idx,
                    distance,
                    from_variant.info_hash,
                    candidate_to.info_hash,
                    self.id_ranks[from_idx],
                    self.id_ranks[to_idx],
                ));
            }

            if !from_variant.passes_size_similarity(candidate_to, self.args.min_size_similarity) {
                log::trace!(
                    "Skipping edge {} -> {}: fails size similarity (threshold: {}). Caching fail.",
                    from_label,
                    to_label,
                    self.args.min_size_similarity
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            // Condition: Skip if sequence similarity check fails
            let similarity_passed;
            let threshold_for_check;
            if apply_tr_rules {
                threshold_for_check = self.args.tr_min_sequence_similarity;
                similarity_passed = from_variant.passes_seq_similarity(
                    candidate_to,
                    &mut self.ed_aligner,
                    threshold_for_check,
                );
            } else {
                threshold_for_check = self.args.min_sequence_similarity;
                similarity_passed = from_variant.passes_seq_similarity(
                    candidate_to,
                    &mut self.ed_aligner,
                    threshold_for_check,
                );
            }

            if !similarity_passed {
                log::trace!(
                    "Skipping edge {} -> {}: fails sequence similarity (threshold: {}). TR lenient: {}. Caching fail.",
                    from_label,
                    to_label,
                    threshold_for_check,
                    apply_tr_rules
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            // Condition: Skip if overlap check fails
            let overlap_threshold = if apply_tr_rules {
                self.args.tr_min_recip_overlap
            } else {
                self.args.min_recip_overlap
            };
            if !from_variant.passes_overlap(candidate_to, overlap_threshold) {
                log::trace!(
                    "Skipping edge {} -> {}: fails overlap check (threshold: {}). TR lenient: {}. Caching fail.",
                    from_label,
                    to_label,
                    overlap_threshold,
                    apply_tr_rules
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            log::trace!(
                "Found valid edge {} -> {} (dist: {}). Caching pass.",
                from_label,
                to_label,
                distance
            );
            self.comparison_cache.insert(cache_key, true);
            return Some(Edge::new(
                from_idx,
                to_idx,
                distance,
                from_variant.info_hash,
                candidate_to.info_hash,
                self.id_ranks[from_idx],
                self.id_ranks[to_idx],
            ));
        }
    }
}
