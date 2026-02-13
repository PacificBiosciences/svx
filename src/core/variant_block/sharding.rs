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

pub(super) fn independent_shards(
    variants: &[VariantInternal],
    args: &MergeArgsInner,
) -> Vec<Vec<usize>> {
    if variants.is_empty() {
        return Vec::new();
    }

    let distance_offset = DISTANCE_OFFSET;
    let mut bounds: Vec<(usize, f64, f64)> = variants
        .iter()
        .enumerate()
        .map(|(idx, variant)| {
            let reach = variant_shard_reach(variant, args) + distance_offset;
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

pub(super) fn coalesce_small_shards(
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

pub(super) fn shard_size_stats_log_message(
    variant_type: SvType,
    contig: &str,
    min_shard_size: usize,
    label: &str,
    shards: &[Vec<usize>],
) -> Option<String> {
    if shards.is_empty() {
        return None;
    }

    let mut shard_sizes: Vec<usize> = shards.iter().map(Vec::len).collect();
    shard_sizes.sort_unstable();
    let n_shards = shard_sizes.len();
    let total_variants: usize = shard_sizes.iter().sum();
    let min_size = shard_sizes[0];
    let p50_size = shard_sizes[(n_shards - 1) / 2];
    let p90_size = shard_sizes[((n_shards - 1) * 9) / 10];
    let max_size = shard_sizes[n_shards - 1];
    let singletons = shard_sizes.iter().filter(|&&size| size == 1).count();

    let small_summary = if min_shard_size > 1 {
        let n_small = shard_sizes
            .iter()
            .filter(|&&size| size < min_shard_size)
            .count();
        format!(
            ", lt_{}={}",
            format_number_with_commas(min_shard_size),
            format_number_with_commas(n_small)
        )
    } else {
        String::new()
    };

    Some(format!(
        "Merge: {} shard size stats ({}) on contig {}: shards={}, variants={}, min={}, p50={}, p90={}, max={}, singletons={}{}",
        variant_type,
        label,
        contig,
        format_number_with_commas(n_shards),
        format_number_with_commas(total_variants),
        format_number_with_commas(min_size),
        format_number_with_commas(p50_size),
        format_number_with_commas(p90_size),
        format_number_with_commas(max_size),
        format_number_with_commas(singletons),
        small_summary
    ))
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
    let mut shard_variants: Vec<VariantInternal> = shard_indices
        .iter()
        .map(|&global_idx| variants[global_idx].clone())
        .collect();
    for (local_idx, variant) in shard_variants.iter_mut().enumerate() {
        variant.index = local_idx;
    }

    let tree = VariantKdTree::new(&shard_variants);
    let mut forest = Forest::new(&shard_variants, args.allow_intrasample);

    let mut merger = VariantMerger::new(&shard_variants, &tree, &mut forest, args);
    merger.execute();
    drop(merger);

    let groups = collect_groups_for_shard(shard_indices, &mut forest);
    ShardMergeResult { groups }
}

pub(super) fn merge_independent_shards(
    forest: &mut Forest,
    variants: &[VariantInternal],
    args: &MergeArgsInner,
    shards: Vec<Vec<usize>>,
) {
    let shard_results: Vec<ShardMergeResult> = shards
        .par_iter()
        .map(|shard_indices| merge_single_shard(shard_indices, variants, args))
        .collect();

    for shard_result in shard_results {
        for group in shard_result.groups {
            if let Some((&anchor, rest)) = group.split_first() {
                for &idx in rest {
                    forest.union_unchecked(anchor, idx);
                }
            }
        }
    }
}

pub(super) fn should_use_sharded_merge(shard_count: usize, no_shard: bool) -> bool {
    shard_count > 1 && !no_shard
}

pub(super) fn shard_count_log_message(
    variant_type: SvType,
    contig: &str,
    shard_count: usize,
) -> String {
    let shard_label = if shard_count == 1 { "shard" } else { "shards" };
    format!(
        "Merge: {} block for contig {} split into {} independent {}",
        variant_type,
        contig,
        format_number_with_commas(shard_count),
        shard_label
    )
}
