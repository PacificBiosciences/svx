use super::*;
use crate::cli::MergeConstraint;
use crate::commands::merge::VariantBlob;
use crate::core::variant::test_utils;
use crate::io::bed_reader::TrId;
use crate::utils::util::init_logger;
use rust_htslib::bcf::Read;
use std::collections::BinaryHeap;
use std::{fs, time::SystemTime};

fn init_test_env() {
    init_logger();
}

fn make_temp_vcf(contents: &str) -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    let nanos = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    path.push(format!("svx_test_variant_block_bnd_{nanos}.vcf"));
    fs::write(&path, contents).unwrap();
    path
}

fn canonical_group_ids(groups: Vec<Vec<VariantInternal>>) -> Vec<Vec<String>> {
    let mut canonical: Vec<Vec<String>> = groups
        .into_iter()
        .filter(|group| !group.is_empty())
        .map(|group| {
            let mut ids: Vec<String> = group.into_iter().map(|v| v.id).collect();
            ids.sort();
            ids
        })
        .collect();
    canonical.sort();
    canonical
}

fn merge_groups_with_serial_variant_merger(
    mut variants: Vec<VariantInternal>,
    args: &MergeArgsInner,
) -> Vec<Vec<String>> {
    for (idx, variant) in variants.iter_mut().enumerate() {
        variant.index = idx;
    }
    let tree = VariantKdTree::new(&variants);
    let mut forest = Forest::new(&variants, args.allow_intrasample);
    let mut merger = VariantMerger::new(&variants, &tree, &mut forest, args);
    merger.execute();

    let mut groups: Vec<Vec<VariantInternal>> = vec![Vec::new(); variants.len()];
    for (i, variant) in variants.iter().enumerate() {
        let root = if forest.parent[i] < 0 {
            i
        } else {
            forest.find(i)
        };
        groups[root].push(variant.clone());
    }

    canonical_group_ids(groups)
}

pub fn create_variants(svtype: SvType) -> Vec<VariantInternal> {
    test_utils::insertion_fixture_variants(svtype)
}

#[test]
fn test_variant_label_formats_index_and_id() {
    init_test_env();

    let label = VariantLabel {
        idx: 631,
        id: "some_name",
    };
    assert_eq!(label.to_string(), "631 (some_name)");
}

#[test]
fn test_init_schedules_self_edge_per_variant() {
    init_test_env();

    let mut variants = vec![
        test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
        test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 10.0, 10.0).unwrap(),
        test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 20.0, 20.0).unwrap(),
    ];
    for (i, v) in variants.iter_mut().enumerate() {
        v.index = i;
    }

    let tree = VariantKdTree::new(&variants);
    let mut forest = Forest::new(&variants, false);
    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
    });

    let merger = VariantMerger::new(&variants, &tree, &mut forest, &args);

    assert_eq!(merger.heap.len(), variants.len());
    for edge in merger.heap.iter() {
        assert_eq!(edge.from, edge.to);
        assert!(edge.dist.abs() <= 1e-12);
    }
}

#[test]
fn test_edge_tie_break_prefers_lexicographic_id_over_index() {
    init_test_env();

    let mut heap = BinaryHeap::new();

    // Pretend:
    // - variant 0 has ID "a" (rank 0)
    // - variant 1 has ID "b" (rank 1)
    // - variant 2 has ID "c" (rank 2)
    // and all INFO hashes are equal (so ID ordering decides the tie).
    heap.push(Edge::new(0, 2, 1.0, 0, 0, 0, 2));
    heap.push(Edge::new(1, 2, 1.0, 0, 0, 1, 2));

    let first = heap.pop().expect("expected an edge");

    assert_eq!(first.from, 0);
}

#[test]
fn test_edge_tie_break_prefers_smaller_info_hash() {
    init_test_env();

    let mut heap = BinaryHeap::new();

    // Same distance. INFO-hash tie-break should prefer the smaller hash, regardless of indices.
    heap.push(Edge::new(0, 2, 1.0, 10, 0, 0, 2));
    heap.push(Edge::new(1, 2, 1.0, 5, 0, 1, 2));

    let first = heap.pop().expect("expected an edge");
    assert_eq!(first.from, 1);
}

#[test]
fn test_edge_display_has_closing_bracket() {
    init_test_env();

    let edge = Edge::new(0, 1, 2.5, 0, 0, 0, 1);
    assert_eq!(edge.to_string(), "Edge[from=0, to=1, dist=2.5]");
}

#[test]
fn test_next_edge_emits_intrasample_candidate() {
    init_test_env();

    let mut variants = vec![
        test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
        test_utils::from_parts(0, "v1".to_string(), SvType::INSERTION, 1.0, 0.0).unwrap(),
        test_utils::from_parts(1, "v2".to_string(), SvType::INSERTION, 50.0, 0.0).unwrap(),
    ];
    for (i, v) in variants.iter_mut().enumerate() {
        v.index = i;
    }

    let tree = VariantKdTree::new(&variants);
    let mut forest = Forest::new(&variants, false);
    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
    });

    let mut merger = VariantMerger::new(&variants, &tree, &mut forest, &args);

    // Simulate "after the sentinel/self edge was popped": the next candidate should be v1.
    merger.search_state[0].next_neighbor_idx = 1;
    let edge = merger.find_next_edge_for_variant(0).unwrap();
    assert_eq!(edge.from, 0);
    assert_eq!(edge.to, 1);
}

#[test]
fn merge_constraint_modes_handle_bridging_chain() {
    init_test_env();

    // Construct a simple "bridge" / chain:
    // v0 ~ v1 and v1 ~ v2 are within max_dist, but v0 ~ v2 is not.
    //
    // Single-linkage clustering will merge all three via v1 unless we add an extra constraint.
    fn run_chain(constraint: MergeConstraint) -> usize {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 5.0, 5.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 10.0, 10.0).unwrap(),
        ];

        for v in &mut variants {
            v.max_dist = 8.0;
        }

        let args = MergeArgsInner::default().with(|a| {
            a.allow_intrasample = false;
            a.knn_search_k = 4;
            a.min_sequence_similarity = 0.0;
            a.require_mutual_distance = true;
            a.merge_constraint = constraint;
        });

        let blob = VariantBlob {
            variants,
            contig: "chr1".to_string(),
            variant_type: SvType::INSERTION,
            args,
            dump_writer: None,
        };

        let mut block = VariantBlock::new(blob);
        block.merge_block();

        let result = block.get_groups();
        result.groups.iter().filter(|g| !g.is_empty()).count()
    }

    assert_eq!(run_chain(MergeConstraint::None), 1);
    assert_eq!(run_chain(MergeConstraint::BboxDiameter), 2);
    assert_eq!(run_chain(MergeConstraint::Clique), 2);
    assert_eq!(run_chain(MergeConstraint::Centroid), 1);
}

#[test]
fn merge_constraint_centroid_rejects_outlier_component() {
    init_test_env();

    // Single-linkage can merge these via v1, but the centroid of the merged component is too far
    // from v0 given v0's max_dist.
    fn run(constraint: MergeConstraint) -> usize {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 1.0, 1.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 9.0, 9.0).unwrap(),
        ];

        variants[0].max_dist = 2.0;
        variants[1].max_dist = 20.0;
        variants[2].max_dist = 20.0;

        let args = MergeArgsInner::default().with(|a| {
            a.allow_intrasample = false;
            a.knn_search_k = 4;
            a.min_sequence_similarity = 0.0;
            a.require_mutual_distance = true;
            a.merge_constraint = constraint;
        });

        let blob = VariantBlob {
            variants,
            contig: "chr1".to_string(),
            variant_type: SvType::INSERTION,
            args,
            dump_writer: None,
        };

        let mut block = VariantBlock::new(blob);
        block.merge_block();

        let result = block.get_groups();
        result.groups.iter().filter(|g| !g.is_empty()).count()
    }

    assert_eq!(run(MergeConstraint::None), 1);
    assert_eq!(run(MergeConstraint::Centroid), 2);
}

#[test]
fn tr_max_dist_can_merge_variants_beyond_per_variant_max_dist() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 2_000.0, 0.0).unwrap(),
        ];
        for v in &mut variants {
            v.max_dist = 100.0;
            v.trid = Some(TrId {
                id: "TR1".to_string(),
                motif_len: 2,
            });
        }
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
        a.tr_max_dist = 4_000;
        a.tr_min_sequence_similarity = 0.0;
    });

    let blob = VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    };

    let mut block = VariantBlock::new(blob);
    block.merge_block();

    let result = block.get_groups();
    assert_eq!(result.groups.iter().filter(|g| !g.is_empty()).count(), 1);
}

#[test]
fn independent_shards_split_variants_with_large_start_gap() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 3.0, 0.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 100.0, 0.0).unwrap(),
            test_utils::from_parts(3, "v3".to_string(), SvType::INSERTION, 105.0, 0.0).unwrap(),
        ];
        for v in &mut variants {
            v.max_dist = 10.0;
        }
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
    });

    let block = VariantBlock::new(VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    });

    assert_eq!(block.independent_shards(), vec![vec![0, 1], vec![2, 3]]);
}

#[test]
fn independent_shards_use_tr_max_dist_for_tr_contained_variants() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 1_500.0, 0.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 10_000.0, 0.0).unwrap(),
        ];
        variants[0].max_dist = 10.0;
        variants[1].max_dist = 10.0;
        variants[2].max_dist = 10.0;
        variants[0].trid = Some(TrId {
            id: "TR1".to_string(),
            motif_len: 2,
        });
        variants[1].trid = Some(TrId {
            id: "TR1".to_string(),
            motif_len: 2,
        });
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
        a.tr_max_dist = 2_000;
    });

    let block = VariantBlock::new(VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    });

    assert_eq!(block.independent_shards(), vec![vec![0, 1], vec![2]]);
}

#[test]
fn min_size_similarity_blocks_symbolic_length_mismatch_merges() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 1_000.0, 10.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 1_000.0, 100.0).unwrap(),
        ];
        variants[0].max_dist = 100.0;
        variants[1].max_dist = 100.0;
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.8;
        a.min_size_similarity = 0.7;
    });

    let blob = VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    };

    let mut block = VariantBlock::new(blob);
    block.merge_block();

    let result = block.get_groups();
    assert_eq!(result.groups.iter().filter(|g| !g.is_empty()).count(), 2);
}

#[test]
fn coalesce_small_shards_merges_adjacent_shards_to_minimum_size() {
    let shards = vec![vec![0], vec![1], vec![2, 3], vec![4], vec![5]];
    let coalesced = VariantBlock::coalesce_small_shards(shards, 3);
    assert_eq!(coalesced, vec![vec![0, 1, 2, 3, 4, 5]]);
}

#[test]
fn coalesce_small_shards_keeps_single_small_shard_when_it_is_the_only_shard() {
    let shards = vec![vec![0, 1]];
    let coalesced = VariantBlock::coalesce_small_shards(shards, 8);
    assert_eq!(coalesced, vec![vec![0, 1]]);
}

#[test]
fn coalesced_sharded_merge_matches_serial_variant_merger_groups() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 3.0, 0.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 100.0, 0.0).unwrap(),
            test_utils::from_parts(3, "v3".to_string(), SvType::INSERTION, 105.0, 0.0).unwrap(),
            test_utils::from_parts(4, "v4".to_string(), SvType::INSERTION, 1_000.0, 0.0).unwrap(),
            test_utils::from_parts(5, "v5".to_string(), SvType::INSERTION, 1_015.0, 0.0).unwrap(),
        ];
        for v in &mut variants {
            v.max_dist = 10.0;
        }
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 6;
        a.min_sequence_similarity = 0.0;
        a.require_mutual_distance = true;
        a.min_shard_size = 4;
    });

    let expected = merge_groups_with_serial_variant_merger(variants.clone(), &args);

    let mut block = VariantBlock::new(VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    });

    assert_eq!(block.independent_shards().len(), 3);
    assert_eq!(
        VariantBlock::coalesce_small_shards(block.independent_shards(), 4).len(),
        1
    );

    block.merge_block();
    let actual = canonical_group_ids(block.get_groups().groups);
    assert_eq!(actual, expected);
}

#[test]
fn sharded_merge_matches_serial_variant_merger_groups() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 3.0, 0.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 100.0, 0.0).unwrap(),
            test_utils::from_parts(3, "v3".to_string(), SvType::INSERTION, 105.0, 0.0).unwrap(),
            test_utils::from_parts(4, "v4".to_string(), SvType::INSERTION, 1_000.0, 0.0).unwrap(),
            test_utils::from_parts(5, "v5".to_string(), SvType::INSERTION, 1_015.0, 0.0).unwrap(),
        ];
        for v in &mut variants {
            v.max_dist = 10.0;
        }
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 6;
        a.min_sequence_similarity = 0.0;
        a.require_mutual_distance = true;
    });

    let expected = merge_groups_with_serial_variant_merger(variants.clone(), &args);

    let mut block = VariantBlock::new(VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    });

    assert_eq!(block.independent_shards().len(), 3);
    block.merge_block();
    let actual = canonical_group_ids(block.get_groups().groups);

    assert_eq!(actual, expected);
}

#[test]
fn no_shard_flag_disables_sharded_merge_path() {
    init_test_env();

    let variants = {
        let mut variants = vec![
            test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
            test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 3.0, 0.0).unwrap(),
            test_utils::from_parts(2, "v2".to_string(), SvType::INSERTION, 100.0, 0.0).unwrap(),
            test_utils::from_parts(3, "v3".to_string(), SvType::INSERTION, 105.0, 0.0).unwrap(),
        ];
        for v in &mut variants {
            v.max_dist = 10.0;
        }
        variants
    };

    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
        a.no_shard = true;
    });

    let block = VariantBlock::new(VariantBlob {
        variants,
        contig: "chr1".to_string(),
        variant_type: SvType::INSERTION,
        args,
        dump_writer: None,
    });

    assert_eq!(block.independent_shards().len(), 2);
    assert!(!block.should_use_sharded_merge(2));
}

#[test]
fn shard_count_log_message_reports_contig_and_shard_count() {
    init_test_env();

    let variants = vec![
        test_utils::from_parts(0, "v0".to_string(), SvType::INSERTION, 0.0, 0.0).unwrap(),
        test_utils::from_parts(1, "v1".to_string(), SvType::INSERTION, 50.0, 0.0).unwrap(),
    ];

    let block = VariantBlock::new(VariantBlob {
        variants,
        contig: "chr22".to_string(),
        variant_type: SvType::INSERTION,
        args: MergeArgsInner::default(),
        dump_writer: None,
    });

    let message = block.shard_count_log_message(3);
    assert!(message.contains("contig chr22"));
    assert!(message.contains("3 independent shards"));
}

#[test]
fn bnd_variants_from_both_contigs_merge_in_same_graph() {
    init_test_env();

    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tbnd_a\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
chr1\t200\tbnd_b\tA\t]chr9:100]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
";
    let path = make_temp_vcf(vcf);
    let mut reader = rust_htslib::bcf::Reader::from_path(&path).unwrap();
    let args = MergeArgsInner::default().with(|a| {
        a.allow_intrasample = false;
        a.knn_search_k = 4;
        a.min_sequence_similarity = 0.0;
    });

    let mut it = reader.records();
    let r0 = it.next().unwrap().unwrap();
    let r1 = it.next().unwrap().unwrap();

    let b0_a = VariantInternal::from_vcf_record(&r0, 0, "chr9", &args, &None).unwrap();
    let b0_b = VariantInternal::from_vcf_record(&r1, 0, "chr1", &args, &None).unwrap();
    let event0 = VariantInternal::from_bnd_pair(b0_a, b0_b, &args).unwrap();

    let b1_a = VariantInternal::from_vcf_record(&r0, 1, "chr9", &args, &None).unwrap();
    let b1_b = VariantInternal::from_vcf_record(&r1, 1, "chr1", &args, &None).unwrap();
    let event1 = VariantInternal::from_bnd_pair(b1_a, b1_b, &args).unwrap();

    let variant_blob = VariantBlob {
        variants: vec![event0, event1],
        contig: "chr1_chr9_TRA".to_string(),
        variant_type: SvType::BND,
        args,
        dump_writer: None,
    };
    let mut variant_block = VariantBlock::new(variant_blob);
    variant_block.merge_block();
    let res = variant_block.get_groups();

    let mut non_empty_groups: Vec<Vec<VariantInternal>> =
        res.groups.into_iter().filter(|g| !g.is_empty()).collect();
    non_empty_groups.sort_by_key(|g| g.len());
    assert_eq!(non_empty_groups.len(), 1);
    assert_eq!(non_empty_groups[0].len(), 2);
}

#[test]
#[ignore]
fn test_variant_block_basic() {
    init_logger();

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
        dump_writer: None,
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
