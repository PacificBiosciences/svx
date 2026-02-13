use crate::core::variant::VariantInternal;

const K: usize = 2;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Neighbor {
    pub distance: f64,
    pub item: usize,
}

mod kdtree_impl {
    use super::*;
    use std::cmp::Ordering;
    use std::collections::BinaryHeap;

    type NodeRef = u32;

    const NONE: NodeRef = 0;
    const LEAF_TAG: NodeRef = 1 << 31;
    const LEAF_SIZE: usize = 32;

    #[inline]
    fn stem_ref(stem_idx: u32) -> NodeRef {
        debug_assert!(stem_idx < LEAF_TAG - 1);
        stem_idx + 1
    }

    #[inline]
    fn leaf_ref(leaf_idx: u32) -> NodeRef {
        debug_assert!(leaf_idx < LEAF_TAG - 1);
        LEAF_TAG | (leaf_idx + 1)
    }

    #[inline]
    fn is_leaf(node: NodeRef) -> bool {
        (node & LEAF_TAG) != 0
    }

    #[inline]
    fn stem_idx(node: NodeRef) -> usize {
        debug_assert!(node != NONE);
        debug_assert!(!is_leaf(node));
        (node - 1) as usize
    }

    #[inline]
    fn leaf_idx(node: NodeRef) -> usize {
        debug_assert!(node != NONE);
        debug_assert!(is_leaf(node));
        ((node & !LEAF_TAG) - 1) as usize
    }

    #[repr(C)]
    #[derive(Clone, Copy, Debug)]
    struct Node {
        x: f64,
        y: f64,
        item: u32,
        left: NodeRef,
        right: NodeRef,
        axis: u8,
        _pad: [u8; 3],
    }

    #[repr(C)]
    #[derive(Clone, Copy, Debug)]
    struct Leaf {
        items: [u32; LEAF_SIZE],
        len: u16,
        _pad: [u8; 2],
    }

    impl Leaf {
        fn new() -> Self {
            Self {
                items: [0; LEAF_SIZE],
                len: 0,
                _pad: [0; 2],
            }
        }
    }

    #[derive(Clone, Copy, Debug)]
    struct HeapItem {
        distance: f64,
        item: u32,
    }

    impl PartialEq for HeapItem {
        fn eq(&self, other: &Self) -> bool {
            self.distance == other.distance && self.item == other.item
        }
    }

    impl Eq for HeapItem {}

    impl PartialOrd for HeapItem {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for HeapItem {
        fn cmp(&self, other: &Self) -> Ordering {
            self.distance
                .total_cmp(&other.distance)
                .then_with(|| self.item.cmp(&other.item))
        }
    }

    pub struct KnnScratch {
        heap: BinaryHeap<HeapItem>,
        out: Vec<Neighbor>,
    }

    impl KnnScratch {
        pub fn new(max_k: usize) -> Self {
            Self {
                heap: BinaryHeap::with_capacity(max_k),
                out: Vec::with_capacity(max_k),
            }
        }

        fn reset(&mut self) {
            self.heap.clear();
            self.out.clear();
        }
    }

    pub struct VariantKdTree {
        xs: Vec<f64>,
        ys: Vec<f64>,
        nodes: Vec<Node>,
        leaves: Vec<Leaf>,
        root: NodeRef,
        num_variants: u32,
    }

    impl VariantKdTree {
        pub fn new(variants: &[VariantInternal]) -> Self {
            assert!(
                variants.len() <= (LEAF_TAG - 1) as usize,
                "Variant count {} exceeds SVX KD-tree index capacity (u32, leaf-tagged)",
                variants.len()
            );

            let (xs, ys) = Self::coords_from_variants(variants);
            let num_variants = xs.len() as u32;
            let n = xs.len();
            if n == 0 {
                return Self {
                    xs,
                    ys,
                    nodes: Vec::new(),
                    leaves: Vec::new(),
                    root: NONE,
                    num_variants,
                };
            }

            let est_stems = n / (LEAF_SIZE + 1);
            let est_leaves = est_stems + 1;

            let mut nodes: Vec<Node> = Vec::with_capacity(est_stems);
            let mut leaves: Vec<Leaf> = Vec::with_capacity(est_leaves);
            let mut order: Vec<u32> = (0..num_variants).collect();
            let root = Self::build_recursive(&mut order, 0, &mut nodes, &mut leaves, &xs, &ys);

            Self {
                xs,
                ys,
                nodes,
                leaves,
                root,
                num_variants,
            }
        }

        fn coords_from_variants(variants: &[VariantInternal]) -> (Vec<f64>, Vec<f64>) {
            let mut xs: Vec<f64> = Vec::with_capacity(variants.len());
            let mut ys: Vec<f64> = Vec::with_capacity(variants.len());
            for v in variants {
                let p = v.point();
                assert!(
                    p[0].is_finite() && p[1].is_finite(),
                    "KD-tree coordinates must be finite (got [{}, {}])",
                    p[0],
                    p[1]
                );
                xs.push(p[0]);
                ys.push(p[1]);
            }
            (xs, ys)
        }

        fn build_recursive(
            order: &mut [u32],
            axis: u8,
            nodes: &mut Vec<Node>,
            leaves: &mut Vec<Leaf>,
            xs: &[f64],
            ys: &[f64],
        ) -> NodeRef {
            debug_assert_eq!(xs.len(), ys.len());
            if order.is_empty() {
                return NONE;
            }

            if order.len() <= LEAF_SIZE {
                let leaf_storage_idx = leaves.len() as u32;
                let mut leaf = Leaf::new();
                leaf.items[..order.len()].copy_from_slice(order);
                leaf.len = order.len() as u16;
                leaves.push(leaf);
                return leaf_ref(leaf_storage_idx);
            }

            let mid = order.len() / 2;
            if axis == 0 {
                order.select_nth_unstable_by(mid, |a, b| {
                    xs[*a as usize]
                        .total_cmp(&xs[*b as usize])
                        .then_with(|| a.cmp(b))
                });
            } else {
                order.select_nth_unstable_by(mid, |a, b| {
                    ys[*a as usize]
                        .total_cmp(&ys[*b as usize])
                        .then_with(|| a.cmp(b))
                });
            }

            let pivot = order[mid];
            let node_storage_idx = nodes.len() as u32;
            nodes.push(Node {
                x: xs[pivot as usize],
                y: ys[pivot as usize],
                item: pivot,
                left: NONE,
                right: NONE,
                axis,
                _pad: [0; 3],
            });

            let (left, right_with_mid) = order.split_at_mut(mid);
            let (_, right) = right_with_mid.split_first_mut().unwrap();

            let next_axis = axis ^ 1;
            let l = Self::build_recursive(left, next_axis, nodes, leaves, xs, ys);
            let r = Self::build_recursive(right, next_axis, nodes, leaves, xs, ys);
            let node = &mut nodes[node_storage_idx as usize];
            node.left = l;
            node.right = r;

            stem_ref(node_storage_idx)
        }

        pub fn knn_into<'s>(
            &self,
            query: &VariantInternal,
            k: usize,
            scratch: &'s mut KnnScratch,
        ) -> &'s [Neighbor] {
            scratch.reset();

            let k = k.min(self.num_variants as usize);
            if k == 0 {
                return &scratch.out;
            }
            if scratch.heap.capacity() < k {
                scratch.heap.reserve(k - scratch.heap.capacity());
            }
            if scratch.out.capacity() < k {
                scratch.out.reserve(k - scratch.out.capacity());
            }

            let query_point = query.point();
            if self.root != NONE {
                self.search_recursive(self.root, &query_point, k, &mut scratch.heap);
            }

            scratch.out.extend(scratch.heap.iter().map(|h| Neighbor {
                distance: h.distance,
                item: h.item as usize,
            }));
            scratch.out.sort_unstable_by(|a, b| {
                a.distance
                    .total_cmp(&b.distance)
                    .then_with(|| a.item.cmp(&b.item))
            });
            &scratch.out
        }

        pub fn knn(&self, query: &VariantInternal, k: usize) -> Vec<Neighbor> {
            let k = k.min(self.num_variants as usize);
            if k == 0 {
                return Vec::new();
            }

            let mut heap: BinaryHeap<HeapItem> = BinaryHeap::with_capacity(k);
            let query_point = query.point();
            if self.root != NONE {
                self.search_recursive(self.root, &query_point, k, &mut heap);
            }

            let mut out: Vec<HeapItem> = heap.into_vec();
            out.sort_unstable();
            out.into_iter()
                .map(|h| Neighbor {
                    distance: h.distance,
                    item: h.item as usize,
                })
                .collect()
        }

        fn search_recursive(
            &self,
            node: NodeRef,
            query: &[f64; K],
            k: usize,
            heap: &mut BinaryHeap<HeapItem>,
        ) {
            debug_assert_ne!(node, NONE);

            if is_leaf(node) {
                let leaf = &self.leaves[leaf_idx(node)];
                for &item in &leaf.items[..leaf.len as usize] {
                    let dx = self.xs[item as usize] - query[0];
                    let dy = self.ys[item as usize] - query[1];
                    let dist = dx * dx + dy * dy;
                    let new_item = HeapItem {
                        distance: dist,
                        item,
                    };
                    if heap.len() < k {
                        heap.push(new_item);
                    } else if let Some(mut worst) = heap.peek_mut() {
                        if new_item < *worst {
                            *worst = new_item;
                        }
                    }
                }
                return;
            }

            let stem = &self.nodes[stem_idx(node)];
            let dx = stem.x - query[0];
            let dy = stem.y - query[1];
            let dist = dx * dx + dy * dy;

            let new_item = HeapItem {
                distance: dist,
                item: stem.item,
            };
            if heap.len() < k {
                heap.push(new_item);
            } else if let Some(mut worst) = heap.peek_mut() {
                if new_item < *worst {
                    *worst = new_item;
                }
            }

            let diff = if stem.axis == 0 {
                query[0] - stem.x
            } else {
                query[1] - stem.y
            };
            let (near, far) = if diff < 0.0 {
                (stem.left, stem.right)
            } else {
                (stem.right, stem.left)
            };

            if near != NONE {
                self.search_recursive(near, query, k, heap);
            }

            let worst_dist = heap.peek().map_or(f64::INFINITY, |w| w.distance);
            if (heap.len() < k || diff * diff <= worst_dist) && far != NONE {
                self.search_recursive(far, query, k, heap);
            }
        }

        pub fn size(&self) -> usize {
            self.num_variants as usize
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::core::svtype::SvType;
        use crate::core::variant::test_utils;

        #[test]
        fn test_bucketed_build_uses_fewer_stems_than_points() {
            let n = LEAF_SIZE * 4;
            let variants: Vec<VariantInternal> = (0..n)
                .map(|i| {
                    test_utils::from_parts(
                        0,
                        i.to_string(),
                        SvType::INSERTION,
                        (i * 10) as f64,
                        (i * 10 + 1) as f64,
                    )
                    .unwrap()
                })
                .collect();
            let tree = VariantKdTree::new(&variants);
            assert!(
                tree.nodes.len() < variants.len(),
                "expected fewer stems than points with bucketed leaves (stems={}, points={})",
                tree.nodes.len(),
                variants.len()
            );
        }

        #[test]
        fn test_node_is_compact() {
            assert!(
                std::mem::size_of::<Node>() <= 40,
                "Node should stay compact for cache locality (got {} bytes)",
                std::mem::size_of::<Node>()
            );
        }
    }
}

pub use kdtree_impl::KnnScratch;
pub use kdtree_impl::VariantKdTree;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::svtype::SvType;
    use crate::core::variant::test_utils;
    use proptest::prelude::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn brute_knn(variants: &[VariantInternal], query: &VariantInternal, k: usize) -> Vec<Neighbor> {
        let k = k.min(variants.len());
        if k == 0 {
            return Vec::new();
        }

        let qp = query.point();
        let mut all: Vec<Neighbor> = variants
            .iter()
            .enumerate()
            .map(|(idx, v)| {
                let p = v.point();
                let dx = p[0] - qp[0];
                let dy = p[1] - qp[1];
                Neighbor {
                    distance: dx * dx + dy * dy,
                    item: idx,
                }
            })
            .collect();

        all.sort_by(|a, b| {
            a.distance
                .total_cmp(&b.distance)
                .then_with(|| a.item.cmp(&b.item))
        });
        all.truncate(k);
        all
    }

    fn has_unique_items(neighbors: &[Neighbor]) -> bool {
        let mut items: Vec<usize> = neighbors.iter().map(|n| n.item).collect();
        items.sort_unstable();
        items.dedup();
        items.len() == neighbors.len()
    }

    fn build_variants(points: &[(u32, u32)]) -> Vec<VariantInternal> {
        points
            .iter()
            .enumerate()
            .map(|(i, (start, end))| {
                test_utils::from_parts(
                    0,
                    i.to_string(),
                    SvType::INSERTION,
                    f64::from(*start),
                    f64::from(*end),
                )
                .unwrap()
            })
            .collect()
    }

    #[test]
    fn knn_k_zero_returns_empty() {
        let variants = vec![
            test_utils::from_parts(0, "0".to_string(), SvType::INSERTION, 10.0, 5.0).unwrap(),
            test_utils::from_parts(0, "1".to_string(), SvType::INSERTION, 20.0, 5.0).unwrap(),
        ];
        let tree = VariantKdTree::new(&variants);
        let query =
            test_utils::from_parts(0, "q".to_string(), SvType::INSERTION, 12.0, 5.0).unwrap();
        assert!(tree.knn(&query, 0).is_empty());
    }

    #[test]
    fn test_variant_kd_tree() {
        let variants = vec![
            test_utils::from_parts(2, "0".to_string(), SvType::INSERTION, 30.0, 36.0).unwrap(),
            test_utils::from_parts(1, "1".to_string(), SvType::INSERTION, 20.0, 24.0).unwrap(),
            test_utils::from_parts(0, "2".to_string(), SvType::INSERTION, 3.0, 20.0).unwrap(),
            test_utils::from_parts(3, "3".to_string(), SvType::INSERTION, 4.0, 5.0).unwrap(),
        ];

        let kdtree = VariantKdTree::new(&variants);
        let query =
            test_utils::from_parts(3, "0".to_string(), SvType::INSERTION, 4.0, 5.0).unwrap();
        let nearest = kdtree.knn(&query, 3);
        assert_eq!(nearest.len(), 3);
        assert_eq!(nearest[0].item, 3);
        assert_eq!(nearest[1].item, 2);
        assert_eq!(nearest[2].item, 1);
    }

    #[test]
    fn test_knn_ordering() {
        let mut rng = StdRng::seed_from_u64(42);
        let num_points = 1000;
        let k = 50;
        let variants: Vec<VariantInternal> = (0..num_points)
            .map(|i| {
                let start: f64 = rng.random_range(0.0..1000.0);
                let svlen: f64 = rng.random_range(10.0..500.0);
                test_utils::from_parts(0, i.to_string(), SvType::INSERTION, start, svlen).unwrap()
            })
            .collect();

        let kdtree = VariantKdTree::new(&variants);

        let query_start: f64 = rng.random_range(0.0..1000.0);
        let query_svlen: f64 = rng.random_range(10.0..500.0);
        let query_variant = test_utils::from_parts(
            0,
            "query".to_string(),
            SvType::INSERTION,
            query_start,
            query_svlen,
        )
        .unwrap();

        let nearest = kdtree.knn(&query_variant, k);

        assert!(nearest.len() <= k);
        if num_points > 0 {
            assert!(!nearest.is_empty());
        }

        for win in nearest.windows(2) {
            assert!(
                win[0].distance <= win[1].distance,
                "KNN results are not sorted by distance: {} > {}",
                win[0].distance,
                win[1].distance
            );
        }
    }

    #[test]
    fn test_knn_matches_bruteforce() {
        let mut rng = StdRng::seed_from_u64(7);
        let num_points = 400;
        let variants: Vec<VariantInternal> = (0..num_points)
            .map(|i| {
                let start: f64 = rng.random_range(0.0..10_000.0);
                let svlen: f64 = rng.random_range(1.0..2_000.0);
                test_utils::from_parts(0, i.to_string(), SvType::INSERTION, start, svlen).unwrap()
            })
            .collect();
        let tree = VariantKdTree::new(&variants);

        for _ in 0..20 {
            let q_idx = rng.random_range(0..num_points);
            let query = &variants[q_idx];
            let k = rng.random_range(1..=32);
            let kd = tree.knn(query, k);
            let brute = brute_knn(&variants, query, k);
            assert_eq!(
                kd.iter().map(|n| n.item).collect::<Vec<_>>(),
                brute.iter().map(|n| n.item).collect::<Vec<_>>(),
            );
        }
    }

    #[test]
    fn test_knn_into_matches_knn_and_resets() {
        let variants = build_variants(&[(10, 12), (20, 21), (25, 29), (30, 35), (100, 101)]);
        let tree = VariantKdTree::new(&variants);
        let query = variants[2].clone();

        let mut scratch = KnnScratch::new(4);
        let first = tree.knn_into(&query, 4, &mut scratch).to_vec();
        let baseline = tree.knn(&query, 4);
        assert_eq!(first, baseline);

        let second = tree.knn_into(&query, 2, &mut scratch).to_vec();
        let baseline2 = tree.knn(&query, 2);
        assert_eq!(second, baseline2);
    }

    proptest! {
        #![proptest_config(ProptestConfig {
            cases: 64,
            .. ProptestConfig::default()
        })]

        #[test]
        fn prop_knn_basic_invariants_hold(
            points in prop::collection::vec((0u32..1_000_000, 1u32..1_000_000), 0..256),
            query_start in 0u32..1_000_000,
            query_end in 1u32..1_000_000,
            k in 0usize..300,
        ) {
            let variants = build_variants(&points);
            let tree = VariantKdTree::new(&variants);
            prop_assert_eq!(tree.size(), variants.len());

            let query = test_utils::from_parts(
                0,
                "query".to_string(),
                SvType::INSERTION,
                f64::from(query_start),
                f64::from(query_end),
            ).unwrap();

            let nearest = tree.knn(&query, k);
            let expected_len = variants.len().min(k);
            prop_assert_eq!(nearest.len(), expected_len);
            prop_assert!(has_unique_items(&nearest));
            prop_assert!(nearest.iter().all(|n| n.item < variants.len()));

            prop_assert!(nearest.windows(2).all(|w| w[0].distance <= w[1].distance));

            if expected_len > 0 {
                let brute = brute_knn(&variants, &query, k);
                let cutoff = brute[expected_len - 1].distance;

                prop_assert!(nearest.iter().all(|n| n.distance <= cutoff));

                for cand in brute_knn(&variants, &query, variants.len()) {
                    if cand.distance < cutoff {
                        prop_assert!(nearest.iter().any(|n| n.item == cand.item));
                    }
                }
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig {
            cases: 64,
            .. ProptestConfig::default()
        })]

        #[test]
        fn prop_knn_matches_bruteforce_exactly(
            points in prop::collection::vec((0u32..1_000_000, 1u32..1_000_000), 0..256),
            query_start in 0u32..1_000_000,
            query_end in 1u32..1_000_000,
            k in 0usize..300,
        ) {
            let variants = build_variants(&points);
            let tree = VariantKdTree::new(&variants);

            let query = test_utils::from_parts(
                0,
                "query".to_string(),
                SvType::INSERTION,
                f64::from(query_start),
                f64::from(query_end),
            ).unwrap();

            let kd = tree.knn(&query, k);
            let brute = brute_knn(&variants, &query, k);
            prop_assert_eq!(
                kd.iter().map(|n| (n.distance, n.item)).collect::<Vec<_>>(),
                brute.iter().map(|n| (n.distance, n.item)).collect::<Vec<_>>()
            );

            let k_small = k.min(variants.len()).min(16);
            let k_large = (k_small + 1).min(variants.len());
            if k_large > 0 {
                let small = tree.knn(&query, k_small);
                let large = tree.knn(&query, k_large);
                let expected = large.into_iter().take(small.len()).collect::<Vec<_>>();
                prop_assert_eq!(small, expected);
            }
        }
    }
}
