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
    pub groups: Vec<Vec<VariantInternal>>,
    pub contig: String,
    pub variant_type: SvType,
    pub n: usize,
}

pub struct VariantBlock {
    variants: Vec<VariantInternal>,
    tree: VariantKdTree,
    forest: Forest,
    contig: String,
    variant_type: SvType,
    n: usize,
    args: MergeArgsInner,
}

impl VariantBlock {
    fn independent_shards(&self) -> Vec<Vec<usize>> {
        sharding::independent_shards(&self.variants, &self.args)
    }

    fn coalesce_small_shards(shards: Vec<Vec<usize>>, min_shard_size: usize) -> Vec<Vec<usize>> {
        sharding::coalesce_small_shards(shards, min_shard_size)
    }

    fn shard_size_stats_log_message(&self, label: &str, shards: &[Vec<usize>]) -> Option<String> {
        sharding::shard_size_stats_log_message(
            self.variant_type,
            &self.contig,
            self.args.min_shard_size,
            label,
            shards,
        )
    }

    fn merge_independent_shards(&mut self, shards: Vec<Vec<usize>>) {
        sharding::merge_independent_shards(&mut self.forest, &self.variants, &self.args, shards)
    }

    fn should_use_sharded_merge(&self, shard_count: usize) -> bool {
        sharding::should_use_sharded_merge(shard_count, self.args.no_shard)
    }

    fn shard_count_log_message(&self, shard_count: usize) -> String {
        sharding::shard_count_log_message(self.variant_type, &self.contig, shard_count)
    }

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

    pub fn contig(&self) -> &str {
        &self.contig
    }

    pub fn variant_type(&self) -> SvType {
        self.variant_type
    }

    pub fn merge_block(&mut self) {
        log::debug!(
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

        let independent_shards = self.independent_shards();
        let independent_shard_count = independent_shards.len();
        log::debug!("{}", self.shard_count_log_message(independent_shard_count));
        if let Some(message) =
            self.shard_size_stats_log_message("independent", independent_shards.as_slice())
        {
            log::debug!("{message}");
        }

        let shards = Self::coalesce_small_shards(independent_shards, self.args.min_shard_size);
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
            if let Some(message) = self.shard_size_stats_log_message("coalesced", shards.as_slice())
            {
                log::debug!("{message}");
            }
        }

        if self.should_use_sharded_merge(shard_count) {
            log::debug!(
                "Merge: Using independent-shard merge path for {} block on contig {}",
                self.variant_type,
                self.contig
            );
            self.merge_independent_shards(shards);
        } else {
            if shard_count > 1 && self.args.no_shard {
                log::debug!(
                    "Merge: --no-shard set; using single-block merge path for {} block on contig {}",
                    self.variant_type,
                    self.contig
                );
            }
            let mut merger =
                VariantMerger::new(&self.variants, &self.tree, &mut self.forest, &self.args);
            merger.execute();
        }

        log::debug!(
            "Merge: Finished {} block for contig {}",
            self.variant_type,
            self.contig
        );
    }

    pub fn get_groups(&mut self) -> VariantBlockResult {
        let mut groups: Vec<Vec<VariantInternal>> = vec![Vec::new(); self.n];
        let variants = std::mem::take(&mut self.variants);
        for (i, variant) in variants.into_iter().enumerate() {
            let root = if self.forest.parent[i] < 0 {
                i
            } else {
                self.forest.find(i)
            };
            groups[root].push(variant);
        }

        VariantBlockResult {
            groups,
            contig: self.contig.clone(),
            variant_type: self.variant_type,
            n: self.n,
        }
    }
}
