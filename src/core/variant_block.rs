mod edge;
mod label;
mod merger;
mod sharding;

#[cfg(test)]
mod tests;

use super::{
    containers::{forest::Forest, kd_tree::VariantKdTree},
    svtype::SvType,
    variant::VariantInternal,
};
use crate::{
    cli::MergeArgsInner, commands::merge::VariantBlob, utils::util::format_number_with_commas,
};
use merger::VariantMerger;

#[cfg(test)]
use edge::Edge;
#[cfg(test)]
use label::VariantLabel;

pub struct VariantBlockResult {
    pub blob_ordinal: u64,
    pub groups: Vec<Vec<VariantInternal>>,
    pub contig: String,
    pub variant_type: SvType,
    pub n: usize,
}

pub struct VariantBlock {
    variants: Vec<VariantInternal>,
    tree: Option<VariantKdTree>,
    forest: Option<Forest>,
    shard_groups: Option<Vec<Vec<usize>>>,
    pub contig: String,
    pub variant_type: SvType,
    n: usize,
    blob_ordinal: u64,
    args: MergeArgsInner,
}

impl VariantBlock {
    pub fn new(blob: VariantBlob) -> Self {
        let n = blob.variants.len();
        log::debug!(
            "Creating VariantBlock for {} {} variants on contig {}",
            format_number_with_commas(n),
            blob.variant_type,
            blob.contig
        );

        let mut variants_with_indices = blob.variants;
        for (i, variant) in variants_with_indices.iter_mut().enumerate() {
            variant.index = i;
        }

        Self {
            variants: variants_with_indices,
            tree: None,
            forest: None,
            shard_groups: None,
            variant_type: blob.variant_type,
            contig: blob.contig,
            n,
            blob_ordinal: blob.blob_ordinal,
            args: blob.args,
        }
    }

    fn ensure_full_block_merge_state(&mut self) {
        if self.tree.is_some() && self.forest.is_some() {
            return;
        }
        self.tree = Some(VariantKdTree::new(&self.variants));
        self.forest = Some(Forest::new(&self.variants, self.args.allow_intrasample));
    }

    pub fn merge_block(&mut self) {
        if self.n <= 1 {
            log::debug!(
                "Skipping merge for {}: only {} variants",
                self.variant_type,
                self.n
            );
            return;
        }
        log::debug!(
            "Merge: Starting {} block for contig {}",
            self.variant_type,
            self.contig
        );

        let independent_shards = sharding::independent_shards(&self.variants, &self.args);
        let independent_shard_count = independent_shards.len();
        log::debug!(
            "{}",
            sharding::shard_count_log_message(
                self.variant_type,
                &self.contig,
                independent_shard_count
            )
        );

        let shards = sharding::coalesce_small_shards(independent_shards, self.args.min_shard_size);
        let shard_count = shards.len();
        if shard_count != independent_shard_count {
            log::debug!(
                "Merge: Coalesced independent shards for {} block on contig {} with min-shard-size {} ({} -> {})",
                self.variant_type,
                self.contig,
                format_number_with_commas(self.args.min_shard_size),
                format_number_with_commas(independent_shard_count),
                format_number_with_commas(shard_count)
            );
        }

        if sharding::should_use_sharded_merge(shard_count, self.args.no_shard) {
            log::debug!(
                "Merge: Using independent-shard merge path for {} block on contig {}",
                self.variant_type,
                self.contig
            );
            self.shard_groups = Some(sharding::merge_independent_shards(
                &self.variants,
                &self.args,
                shards,
            ));
            self.tree = None;
            self.forest = None;
        } else {
            if shard_count > 1 && self.args.no_shard {
                log::debug!(
                    "Merge: --no-shard set; using single-block merge path for {} block on contig {}",
                    self.variant_type,
                    self.contig
                );
            }
            self.ensure_full_block_merge_state();
            self.shard_groups = None;
            let tree = self.tree.as_ref().expect(
                "VariantBlock invariant violated: tree must be initialized for serial merge",
            );
            let forest = self.forest.as_mut().expect(
                "VariantBlock invariant violated: forest must be initialized for serial merge",
            );
            let mut merger = VariantMerger::new(&self.variants, tree, forest, &self.args);
            merger.execute();
        }

        log::debug!(
            "Merge: Finished {} block for contig {}",
            self.variant_type,
            self.contig
        );
    }

    pub fn get_groups(&mut self) -> VariantBlockResult {
        let variants = std::mem::take(&mut self.variants);
        let groups = if let Some(shard_groups) = self.shard_groups.take() {
            Self::groups_from_shard_indices(self.n, variants, shard_groups)
        } else if let Some(forest) = self.forest.as_mut() {
            Self::groups_from_forest(self.n, variants, forest)
        } else {
            Self::singleton_groups(self.n, variants)
        };

        VariantBlockResult {
            blob_ordinal: self.blob_ordinal,
            groups,
            contig: self.contig.clone(),
            variant_type: self.variant_type,
            n: self.n,
        }
    }

    fn groups_from_forest(
        n: usize,
        variants: Vec<VariantInternal>,
        forest: &mut Forest,
    ) -> Vec<Vec<VariantInternal>> {
        let mut groups: Vec<Vec<VariantInternal>> = vec![Vec::new(); n];
        for (i, variant) in variants.into_iter().enumerate() {
            let root = if forest.parent[i] < 0 {
                i
            } else {
                forest.find(i)
            };
            groups[root].push(variant);
        }
        groups
    }

    fn groups_from_shard_indices(
        n: usize,
        variants: Vec<VariantInternal>,
        shard_groups: Vec<Vec<usize>>,
    ) -> Vec<Vec<VariantInternal>> {
        let mut groups: Vec<Vec<VariantInternal>> = vec![Vec::new(); n];
        let mut variants_by_index: Vec<Option<VariantInternal>> =
            variants.into_iter().map(Some).collect();

        for shard_group in shard_groups {
            let Some(&anchor) = shard_group.first() else {
                continue;
            };
            for idx in shard_group {
                if let Some(variant) = variants_by_index[idx].take() {
                    groups[anchor].push(variant);
                }
            }
        }

        for (idx, variant) in variants_by_index.into_iter().enumerate() {
            if let Some(variant) = variant {
                groups[idx].push(variant);
            }
        }

        groups
    }

    fn singleton_groups(n: usize, variants: Vec<VariantInternal>) -> Vec<Vec<VariantInternal>> {
        let mut groups: Vec<Vec<VariantInternal>> = vec![Vec::new(); n];
        for (idx, variant) in variants.into_iter().enumerate() {
            groups[idx].push(variant);
        }
        groups
    }
}
