use super::merger::VariantMerger;
use crate::{
    DISTANCE_OFFSET,
    cli::MergeArgsInner,
    core::{
        containers::{forest::Forest, kd_tree::VariantKdTree},
        svtype::SvType,
        variant::VariantInternal,
    },
    utils::util::format_number_with_commas,
};
use rayon::prelude::*;

struct ShardMergeResult {
    groups: Vec<Vec<usize>>,
}

fn variant_shard_reach(variant: &VariantInternal, args: &MergeArgsInner) -> f64 {
    if variant.trid.is_some() {
        variant.max_dist.max(f64::from(args.tr_max_dist))
    } else {
        variant.max_dist
    }
}

pub fn independent_shards(variants: &[VariantInternal], args: &MergeArgsInner) -> Vec<Vec<usize>> {
    if variants.is_empty() {
        return Vec::new();
    }

    let mut bounds: Vec<(usize, f64, f64)> = variants
        .iter()
        .enumerate()
        .map(|(idx, variant)| {
            let reach = variant_shard_reach(variant, args) + DISTANCE_OFFSET;
            (idx, variant.start - reach, variant.start + reach)
        })
        .collect();

    bounds.sort_by(|a, b| {
        a.1.total_cmp(&b.1)
            .then_with(|| a.2.total_cmp(&b.2))
            .then_with(|| a.0.cmp(&b.0))
    });

    let mut shards: Vec<Vec<usize>> = Vec::new();
    let mut current_shard: Vec<usize> = Vec::new();
    let mut current_right = f64::NEG_INFINITY;

    for (idx, left, right) in bounds {
        if current_shard.is_empty() {
            current_shard.push(idx);
            current_right = right;
            continue;
        }

        if left > current_right {
            current_shard.sort_unstable();
            shards.push(current_shard);
            current_shard = vec![idx];
            current_right = right;
        } else {
            current_shard.push(idx);
            current_right = current_right.max(right);
        }
    }

    if !current_shard.is_empty() {
        current_shard.sort_unstable();
        shards.push(current_shard);
    }

    shards
}

pub fn coalesce_small_shards(
    mut shards: Vec<Vec<usize>>,
    min_shard_size: usize,
) -> Vec<Vec<usize>> {
    if min_shard_size <= 1 || shards.len() <= 1 {
        return shards;
    }

    let mut coalesced: Vec<Vec<usize>> = Vec::with_capacity(shards.len());
    let mut pending: Vec<usize> = Vec::new();
    for mut shard in shards.drain(..) {
        if pending.is_empty() && shard.len() >= min_shard_size {
            coalesced.push(shard);
            continue;
        }

        pending.append(&mut shard);
        if pending.len() >= min_shard_size {
            pending.sort_unstable();
            coalesced.push(std::mem::take(&mut pending));
        }
    }

    if !pending.is_empty() {
        if let Some(last) = coalesced.last_mut() {
            last.append(&mut pending);
            last.sort_unstable();
        } else {
            pending.sort_unstable();
            coalesced.push(pending);
        }
    }

    coalesced
}

fn collect_groups_for_shard(shard_indices: &[usize], forest: &mut Forest) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); shard_indices.len()];
    for (local_idx, &global_idx) in shard_indices.iter().enumerate() {
        let root = if forest.parent[local_idx] < 0 {
            local_idx
        } else {
            forest.find(local_idx)
        };
        groups[root].push(global_idx);
    }
    groups
        .into_iter()
        .filter(|group| !group.is_empty())
        .collect()
}

fn merge_single_shard(
    shard_indices: &[usize],
    variants: &[VariantInternal],
    args: &MergeArgsInner,
) -> ShardMergeResult {
    let shard_variants: Vec<&VariantInternal> = shard_indices
        .iter()
        .map(|&global_idx| &variants[global_idx])
        .collect();

    let tree = VariantKdTree::new_from_refs(&shard_variants);
    let mut forest = Forest::new_from_refs(&shard_variants, args.allow_intrasample);

    let mut merger = VariantMerger::new(&shard_variants, &tree, &mut forest, args);
    merger.execute();
    drop(merger);

    let groups = collect_groups_for_shard(shard_indices, &mut forest);
    ShardMergeResult { groups }
}

pub fn merge_independent_shards(
    variants: &[VariantInternal],
    args: &MergeArgsInner,
    shards: Vec<Vec<usize>>,
) -> Vec<Vec<usize>> {
    let shard_results: Vec<ShardMergeResult> = shards
        .par_iter()
        .map(|shard_indices| merge_single_shard(shard_indices, variants, args))
        .collect();

    shard_results
        .into_iter()
        .flat_map(|shard_result| shard_result.groups)
        .collect()
}

pub fn should_use_sharded_merge(shard_count: usize, no_shard: bool) -> bool {
    shard_count > 1 && !no_shard
}

pub fn shard_count_log_message(variant_type: SvType, contig: &str, shard_count: usize) -> String {
    format!(
        "Merge: {} block for contig {} split into {} independent shards",
        variant_type,
        contig,
        format_number_with_commas(shard_count),
    )
}
