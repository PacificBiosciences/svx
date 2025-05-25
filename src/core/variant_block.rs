use super::{
    aligner::WFAligner,
    containers::{forest::Forest, kd_tree::VariantKdTree},
    svtype::SvType,
    variant_internal::VariantInternal,
};
use crate::{
    cli::MergeArgsInner,
    commands::merge::VariantBlob,
    core::aligner::{AlignmentScope, MemoryModel},
    utils::util::format_number_with_commas,
    DISTANCE_OFFSET,
};
use std::{
    cmp::Ordering,
    collections::{BinaryHeap, HashMap},
    fmt,
};

pub struct VariantBlockResult {
    pub groups: Vec<Vec<VariantInternal>>,
    pub contig: String,
    pub variant_type: SvType,
    pub n: usize,
}

pub struct VariantBlock {
    pub variants: Vec<VariantInternal>,
    pub tree: VariantKdTree,
    pub forest: Forest,
    pub contig: String,
    pub variant_type: SvType,
    pub n: usize,
    pub args: MergeArgsInner,
}

impl VariantBlock {
    pub fn new(variant_blob: VariantBlob) -> Self {
        let n = variant_blob.variants.len();
        log::debug!(
            "Creating VariantBlock for {} {} variants on contig {}",
            format_number_with_commas(n),
            variant_blob.variant_type,
            variant_blob.contig
        );

        let mut variants_with_indices = variant_blob.variants;
        for (i, variant) in variants_with_indices.iter_mut().enumerate() {
            variant.index = i;
        }

        let tree = VariantKdTree::new(&variants_with_indices);
        if variant_blob.args.dump {
            tree.dump_coordinates(&variants_with_indices);
        }

        let forest = Forest::new(&variants_with_indices, variant_blob.args.allow_intrasample);

        Self {
            variants: variants_with_indices,
            tree,
            forest,
            variant_type: variant_blob.variant_type,
            contig: variant_blob.contig,
            n,
            args: variant_blob.args,
        }
    }

    fn top_n_counts<K: Eq + std::hash::Hash + Clone>(
        map: &HashMap<K, usize>,
        n: usize,
    ) -> Vec<(K, usize)> {
        let mut vec: Vec<(K, usize)> = map.iter().map(|(k, &v)| (k.clone(), v)).collect();
        vec.sort_by(|a, b| b.1.cmp(&a.1));
        vec.into_iter().take(n).collect()
    }

    fn log_merge_statistics(
        candidate_counts: &HashMap<usize, usize>,
        comparison_counts: &HashMap<(usize, usize), usize>,
        variants: &[VariantInternal],
    ) {
        let top_candidates = Self::top_n_counts(candidate_counts, 10);
        if !top_candidates.is_empty() {
            log::debug!("Top 10 most frequent merge candidates:");
            for (i, (candidate_idx, count)) in top_candidates.iter().enumerate() {
                log::debug!(
                    "top={} {}: {} times candidate",
                    i + 1,
                    &variants[*candidate_idx],
                    count
                );
            }
        }

        let total_comparisons: usize = comparison_counts.values().sum();
        let redundant_comparisons: usize = comparison_counts
            .values()
            .filter(|&&count| count > 1)
            .map(|&count| count - 1)
            .sum();

        let top_comparisons = Self::top_n_counts(comparison_counts, 10);
        if !top_comparisons.is_empty() && redundant_comparisons > 0 {
            log::debug!("Top 10 most frequent variant comparisons:");
            for (i, (pair, count)) in top_comparisons.iter().enumerate() {
                log::debug!(
                    "  top={} pair=({}, {}): {} times",
                    i + 1,
                    pair.0,
                    pair.1,
                    count
                );
            }
        }

        if total_comparisons > 0 {
            log::debug!(
                "Total comparisons made: {}",
                format_number_with_commas(total_comparisons)
            );

            if redundant_comparisons > 0 {
                let redundancy_percent =
                    (redundant_comparisons as f64 / total_comparisons as f64) * 100.0;
                log::debug!(
                    "Total redundant comparisons (pair compared > 1 time): {} ({:.2}%)",
                    format_number_with_commas(redundant_comparisons),
                    redundancy_percent
                );
            }
        }
    }

    pub fn merge_block(&mut self) {
        log::info!(
            "Merge: Starting {} block for contig {}",
            self.variant_type,
            self.contig
        );
        if self.n <= 1 {
            log::debug!(
                "Skipping merge for {}: only {} variants",
                self.variant_type,
                self.n
            );
            return;
        }

        let candidate_counts;
        let comparison_counts;

        {
            let mut merger =
                VariantMerger::new(&self.variants, &self.tree, &mut self.forest, &self.args);
            merger.execute();
            candidate_counts = merger.candidate_counts;
            comparison_counts = merger.comparison_counts;
        }

        if self.args.refinement_iterations > 0 {
            self.refine_clusters();
        }

        Self::log_merge_statistics(&candidate_counts, &comparison_counts, &self.variants);

        log::info!(
            "Merge: Finished {} block for contig {}",
            self.variant_type,
            self.contig
        );
    }

    pub fn get_groups(&mut self) -> VariantBlockResult {
        let mut groups: Vec<Vec<VariantInternal>> = vec![Vec::new(); self.n];
        for i in 0..self.n {
            if self.forest.parent[i] < 0 {
                groups[i].push(self.variants[i].clone());
            } else {
                let root = self.forest.find(i);
                groups[root].push(self.variants[i].clone());
            }
        }

        VariantBlockResult {
            groups,
            contig: self.contig.clone(),
            variant_type: self.variant_type,
            n: self.n,
        }
    }

    fn refine_clusters(&mut self) {
        if self.args.refinement_iterations == 0 {
            return;
        }
        log::debug!(
            "Refine: Starting iterative refinement for {} block on contig {} (max {} iterations).",
            self.variant_type,
            self.contig,
            self.args.refinement_iterations
        );

        for iteration in 0..self.args.refinement_iterations {
            let mut groups_map: HashMap<usize, Vec<usize>> = HashMap::new();
            for i in 0..self.n {
                let root = self.forest.find(i);
                groups_map.entry(root).or_default().push(i);
            }

            if groups_map.len() <= 1 {
                log::debug!(
                    "Refine iter {}: Skipping, {} clusters found.",
                    iteration + 1,
                    groups_map.len()
                );
                break;
            }

            let mut centroids: Vec<(usize, [f64; 2])> = Vec::new(); // (root_id, point)
            for (root_id, variant_indices) in &groups_map {
                if variant_indices.is_empty() {
                    continue;
                }
                let mut sum_point = [0.0, 0.0];
                for &variant_idx in variant_indices {
                    let p = self.variants[variant_idx].point();
                    sum_point[0] += p[0];
                    sum_point[1] += p[1];
                }
                let count = variant_indices.len() as f64;
                centroids.push((*root_id, [sum_point[0] / count, sum_point[1] / count]));
            }

            let mut new_assignments: Vec<usize> = vec![0; self.n];
            let mut changed_count = 0;

            for (i, assignment_slot) in new_assignments.iter_mut().enumerate() {
                let variant = &self.variants[i];
                let current_root_from_prev_state = self.forest.find(i);

                let mut best_overall_target_root = current_root_from_prev_state; // Default to not moving
                let mut min_overall_dist = f32::MAX;

                // Find the best possible cluster for this variant among all centroids
                for (candidate_root_id, candidate_centroid_point) in &centroids {
                    let dist_to_candidate =
                        VariantInternal::point_distance(&variant.point(), candidate_centroid_point);

                    let mut consider_this_candidate = false;

                    if dist_to_candidate < min_overall_dist {
                        consider_this_candidate = true;
                    } else if dist_to_candidate == min_overall_dist {
                        // Tie-breaking: prefer smaller root_id for determinism
                        if *candidate_root_id < best_overall_target_root {
                            consider_this_candidate = true;
                        }
                    }

                    if consider_this_candidate {
                        let mut compatible = true;
                        if !self.args.allow_intrasample {
                            if let Some(target_group_members) = groups_map.get(candidate_root_id) {
                                for &member_idx in target_group_members {
                                    // Ensure we are not checking the variant against itself for intrasample compatibility
                                    if member_idx == i {
                                        continue;
                                    }
                                    if self.variants[member_idx].sample_id == variant.sample_id {
                                        compatible = false;
                                        break;
                                    }
                                }
                            } else {
                                // This case should ideally not be reached if centroids are derived from groups_map
                                log::error!(
                                    "Refine iter {}: candidate_root_id {} not found in groups_map during compatibility check.",
                                    iteration + 1,
                                    candidate_root_id
                                );
                                compatible = false; // Behave safely
                            }
                        }

                        if compatible {
                            min_overall_dist = dist_to_candidate;
                            best_overall_target_root = *candidate_root_id;
                        }
                    }
                }

                if min_overall_dist == f32::MAX && !centroids.is_empty() {
                    log::warn!(
                        "Refine iter {}: Variant {} (current root {}) could not find any compatible cluster or all distances were MAX. Staying in original cluster (or assigned default).",
                        iteration + 1,
                        variant.id,
                        current_root_from_prev_state
                    );
                }
                *assignment_slot = best_overall_target_root;
                if best_overall_target_root != current_root_from_prev_state {
                    changed_count += 1;
                }
            }

            if changed_count == 0 {
                log::debug!(
                    "Refine iter {}: No variants changed clusters. Stopping refinement.",
                    iteration + 1
                );
                break;
            } else {
                log::debug!(
                    "Refine iter {}: {} variants reassigned to new clusters.",
                    iteration + 1,
                    changed_count
                );
            }

            self.forest = Forest::new(&self.variants, self.args.allow_intrasample); // Reset forest
            let mut representative_map: HashMap<usize, usize> = HashMap::new();

            for (i, &conceptual_group_id) in new_assignments.iter().enumerate() {
                if let Some(&representative_idx) = representative_map.get(&conceptual_group_id) {
                    self.forest.union(i, representative_idx);
                } else {
                    representative_map.insert(conceptual_group_id, i);
                }
            }
        }
        log::info!(
            "Refine: Finished iterative refinement for {} block on contig {}.",
            self.variant_type,
            self.contig
        );
    }
}

struct NeighborSearchState<'a> {
    neighbors: Vec<&'a VariantInternal>,
    next_neighbor_idx: usize,
}

struct VariantMerger<'a> {
    variants: &'a [VariantInternal],
    tree: &'a VariantKdTree,
    forest: &'a mut Forest,
    args: &'a MergeArgsInner,
    ed_aligner: WFAligner,
    search_state: Vec<NeighborSearchState<'a>>,
    heap: BinaryHeap<Edge>,
    candidate_counts: HashMap<usize, usize>, // DEBUG: Tracks how often a given variant is a candidate
    comparison_counts: HashMap<(usize, usize), usize>, // DEBUG: Track how often each variant pair is compared
    comparison_cache: HashMap<(usize, usize), bool>, // Cache for comparison results (true=pass, false=fail)
}

impl<'a> VariantMerger<'a> {
    fn new(
        variants: &'a [VariantInternal],
        tree: &'a VariantKdTree,
        forest: &'a mut Forest,
        args: &'a MergeArgsInner,
    ) -> Self {
        let n = variants.len();
        let ed_aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
            .edit()
            .build();

        let mut search_state: Vec<NeighborSearchState<'a>> = Vec::with_capacity(n);
        for variant in variants.iter() {
            let initial_neighbors =
                Self::find_nearest_neighbors_internal(tree, variants, variant, 10);
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
            ed_aligner,
            search_state,
            heap: BinaryHeap::new(),
            candidate_counts: HashMap::new(),
            comparison_counts: HashMap::new(),
            comparison_cache: HashMap::new(),
        };
        merger.initialize_heap();
        merger
    }

    fn find_nearest_neighbors_internal<'b>(
        tree: &VariantKdTree,
        variants: &'b [VariantInternal],
        variant: &VariantInternal,
        k: usize,
    ) -> Vec<&'b VariantInternal> {
        tree.knn(variant, k)
            .iter()
            .map(|neighbor| &variants[neighbor.item as usize])
            .collect()
    }

    fn initialize_heap(&mut self) {
        for i in 0..self.variants.len() {
            if let Some(edge) = self.find_next_edge_for_variant(i) {
                self.heap.push(edge);
            }
        }
    }

    // Finds the next valid edge originating from `from_idx` based on its search state
    fn find_next_edge_for_variant(&mut self, from_idx: usize) -> Option<Edge> {
        let from_variant = &self.variants[from_idx];
        let state = &mut self.search_state[from_idx];
        loop {
            if state.next_neighbor_idx >= state.neighbors.len() {
                let current_k = state.neighbors.len();
                if current_k == 0 || current_k >= self.variants.len() - 1 {
                    log::trace!("Stopping neighbor search for {}: no more potential neighbors (current_k={}, total_variants={}).", from_idx, current_k, self.variants.len());
                    return None;
                }
                let next_k = (current_k * 2).min(self.variants.len() - 1);
                log::trace!(
                    "Finding more neighbors for {}: current_k={}, next_k={}",
                    from_idx,
                    current_k,
                    next_k
                );
                state.neighbors = Self::find_nearest_neighbors_internal(
                    self.tree,
                    self.variants,
                    from_variant,
                    next_k,
                );

                // Check again after potentially fetching more neighbors
                if state.next_neighbor_idx >= state.neighbors.len() {
                    log::trace!("Stopping neighbor search for {}: doubling k ({}) did not yield new neighbors beyond index {}.", from_idx, next_k, state.next_neighbor_idx);
                    return None;
                }
            }

            let candidate_to = state.neighbors[state.next_neighbor_idx];
            state.next_neighbor_idx += 1; // Consume this neighbor *before* any continue or return statements for this iteration

            // Start validation checks
            let distance = from_variant.distance(candidate_to);

            // Condition: Stop searching further for `from_idx` if the candidate is too far (exceeds own max_dist threshold)
            // This check is NOT symmetric and determines if we stop searching *for from_idx*, so it comes before caching.
            if distance > from_variant.max_dist + DISTANCE_OFFSET {
                log::trace!(
                    "Stopping search for {}: candidate {} is too far ({} > {} + offset).",
                    from_idx,
                    candidate_to.index,
                    distance,
                    from_variant.max_dist
                );
                return None;
            }

            // Cache check
            let to_idx = candidate_to.index;
            let cache_key = (from_idx.min(to_idx), from_idx.max(to_idx));
            if let Some(&cached_result) = self.comparison_cache.get(&cache_key) {
                if cached_result {
                    log::trace!(
                        "Found valid edge {} -> {} (dist: {}, cached pass).",
                        from_idx,
                        to_idx,
                        distance
                    );
                    return Some(Edge::new(from_idx, to_idx, distance));
                } else {
                    log::trace!(
                        "Skipping edge {} -> {}: failed cached checks.",
                        from_idx,
                        to_idx
                    );
                    continue;
                }
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
                        from_idx,
                        to_idx,
                        from_variant.trid.as_ref().unwrap().id,
                        candidate_to.trid.as_ref().unwrap().id
                    );
                    self.comparison_cache.insert(cache_key, false);
                    continue;
                }
                apply_tr_rules = true;
            }

            // NOTE: These are both used for debugging
            // Realistically with the cache check all comparisons are now non-redundant
            *self.candidate_counts.entry(candidate_to.index).or_insert(0) += 1;
            *self
                .comparison_counts
                .entry((
                    from_idx.min(candidate_to.index),
                    from_idx.max(candidate_to.index),
                ))
                .or_insert(0) += 1;
            //

            // If not in the cache, perform symmetric checks and cache the result
            let effective_max_dist_for_pair = if apply_tr_rules {
                self.args.tr_max_dist as f32
            } else if self.args.require_mutual_distance {
                from_variant.max_dist.min(candidate_to.max_dist)
            } else {
                from_variant.max_dist.max(candidate_to.max_dist)
            };

            if distance > effective_max_dist_for_pair + DISTANCE_OFFSET {
                log::trace!(
                    "Skipping edge {} -> {}: exceeds effective max distance for pair ({} > {} + offset). TR lenient: {}. Caching fail.",
                    from_idx,
                    to_idx,
                    distance,
                    effective_max_dist_for_pair,
                    apply_tr_rules
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            // Condition: Skip if intrasample merging is disallowed and samples are the same, independent of TR status.
            if !self.args.allow_intrasample && from_variant.sample_id == candidate_to.sample_id {
                log::trace!(
                    "Skipping intrasample edge {} -> {} (sample_id: {}). Caching fail.",
                    from_idx,
                    to_idx,
                    from_variant.sample_id
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            // Condition: Skip if sequence similarity check fails
            let similarity_passed;
            let threshold_for_check;
            // if apply_tr_rules && from_variant.trid.as_ref().unwrap().motif_len < 8 {
            if apply_tr_rules {
                threshold_for_check = self.args.tr_min_sequence_similarity;
                // similarity_passed =
                //     from_variant.passes_kmer_jaccard_similarity(candidate_to, threshold_for_check);
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
            };

            if !similarity_passed {
                log::trace!(
                    "Skipping edge {} -> {}: fails sequence similarity (threshold: {}). TR lenient: {}. Caching fail.",
                    from_idx,
                    to_idx,
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
                    from_idx,
                    to_idx,
                    overlap_threshold,
                    apply_tr_rules
                );
                self.comparison_cache.insert(cache_key, false);
                continue;
            }

            log::trace!(
                "Found valid edge {} -> {} (dist: {}). Caching pass.",
                from_idx,
                to_idx,
                distance
            );
            self.comparison_cache.insert(cache_key, true);
            return Some(Edge::new(from_idx, to_idx, distance));
        }
    }

    fn execute(&mut self) {
        let n = self.variants.len();
        if n <= 1 {
            return;
        }

        while let Some(edge) = self.heap.pop() {
            // Attempt to union the sets containing the variants of the edge
            if self.forest.can_union(edge.from, edge.to) {
                log::trace!("Uniting {} and {}", edge.from, edge.to);
                self.forest.union(edge.from, edge.to);
            } else {
                log::trace!(
                    "Cannot union {} and {} (already connected or sample conflict)",
                    edge.from,
                    edge.to
                );
            }

            // Find the *next* valid edge originating from `edge.from` and add it to the heap
            if let Some(next_edge) = self.find_next_edge_for_variant(edge.from) {
                log::trace!("Pushing next edge for {}: {}", edge.from, next_edge);
                self.heap.push(next_edge);
            } else {
                log::trace!("No more valid edges found for {}", edge.from);
            }
        }
        log::debug!("Merging heap empty, process finished.");
    }
}

#[derive(Debug)]
struct Edge {
    from: usize,
    to: usize,
    dist: f32,
}

impl Edge {
    fn new(from: usize, to: usize, dist: f32) -> Self {
        Self { from, to, dist }
    }
}

impl fmt::Display for Edge {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Edge[from={}, to={}, dist={}",
            self.from, self.to, self.dist
        )
    }
}

impl Ord for Edge {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .dist
            .partial_cmp(&self.dist)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.from.cmp(&other.from))
            .then_with(|| self.to.cmp(&other.to))
    }
}

impl PartialOrd for Edge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Edge {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{init_config, utils::util::init_logger, SvxConfig};

    pub fn create_variants(svtype: SvType) -> Vec<VariantInternal> {
        vec![
            VariantInternal::from_parts(0, "var1".to_string(), svtype, 10.0, 5.0).unwrap(),
            VariantInternal::from_parts(0, "var2".to_string(), svtype, 1.0, 5.0).unwrap(),
            VariantInternal::from_parts(0, "var3".to_string(), svtype, 18.0, 5.0).unwrap(),
            VariantInternal::from_parts(1, "var4".to_string(), svtype, 12.0, 7.0).unwrap(),
            VariantInternal::from_parts(1, "var5".to_string(), svtype, 10.0, 5.0).unwrap(),
            VariantInternal::from_parts(1, "var6".to_string(), svtype, 30.0, 30.0).unwrap(),
            VariantInternal::from_parts(1, "var7".to_string(), svtype, 0.0, 0.0).unwrap(),
            VariantInternal::from_parts(2, "var8".to_string(), svtype, 12.0, 12.0).unwrap(),
            VariantInternal::from_parts(2, "var9".to_string(), svtype, 15.0, 15.0).unwrap(),
            VariantInternal::from_parts(2, "var10".to_string(), svtype, 20.0, 20.0).unwrap(),
            VariantInternal::from_parts(2, "var11".to_string(), svtype, 28.0, 28.0).unwrap(),
            VariantInternal::from_parts(3, "var12".to_string(), svtype, 25.0, 25.0).unwrap(),
            VariantInternal::from_parts(4, "var13".to_string(), svtype, 22.0, 22.0).unwrap(),
        ]
    }

    #[test]
    #[ignore]
    fn test_variant_block_basic() {
        init_logger();

        init_config(SvxConfig {
            kd_tree_norm: 2,
            dump: false,
        });

        let variants = create_variants(SvType::INSERTION);
        let args = MergeArgsInner {
            min_recip_overlap: 0.0,
            ..Default::default()
        };
        let variant_blob = VariantBlob {
            variants,
            contig: "chr1".to_string(),
            variant_type: SvType::INSERTION,
            args,
        };
        let mut variant_block = VariantBlock::new(variant_blob);
        variant_block.merge_block();
        let vbr = variant_block.get_groups();

        for group in vbr.groups {
            if group.len() > 1 {
                println!("Variant:");
                for v in group {
                    println!(
                        "id: {}, sample: {}, start: {}, svlen: {}",
                        v.id, v.sample_id, v.start, v.svlen
                    );
                }
                println!();
            }
        }
    }
}
