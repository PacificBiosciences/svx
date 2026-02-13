use crate::{DISTANCE_OFFSET, core::variant::VariantInternal};
use std::mem;

// Each component may require multiple 64-bit integers to hold its bitset of sample IDs if there are many samples
// parent[i] - is negative if root, more negative means bigger set; if nonnegative, indicates the parent
// sample_masks - For each node (only guaranteed current at roots), a bitmask of which samples are present in its component.
// Stored node-major: sample_masks[node * mask_words + word].
const SAMPLES_PER_MASK: usize = u64::BITS as usize;

#[derive(Debug)]
pub struct Forest {
    pub parent: Vec<i32>,
    pub sample_masks: Vec<u64>,
    pub mask_words: usize,
    min_start: Vec<f64>,
    max_start: Vec<f64>,
    min_end: Vec<f64>,
    max_end: Vec<f64>,
    allow_intrasample: bool,
}

impl Forest {
    pub fn new(variants: &[VariantInternal], allow_intrasample: bool) -> Self {
        let n = variants.len();
        let max_sample = variants.iter().map(|v| v.sample_id).max().unwrap_or(0);
        let mask_words = if n == 0 {
            0
        } else {
            max_sample / SAMPLES_PER_MASK + 1
        };
        let parent = vec![-1; n];
        let mut sample_masks = vec![0u64; n * mask_words];
        let mut min_start = vec![0.0; n];
        let mut max_start = vec![0.0; n];
        let mut min_end = vec![0.0; n];
        let mut max_end = vec![0.0; n];
        for (i, variant) in variants.iter().enumerate() {
            let mask_id = variant.sample_id / SAMPLES_PER_MASK;
            let bit_position = variant.sample_id % SAMPLES_PER_MASK;
            sample_masks[i * mask_words + mask_id] |= 1u64 << bit_position;

            min_start[i] = variant.start;
            max_start[i] = variant.start;
            min_end[i] = variant.end;
            max_end[i] = variant.end;
        }
        Forest {
            parent,
            sample_masks,
            mask_words,
            min_start,
            max_start,
            min_end,
            max_end,
            allow_intrasample,
        }
    }

    #[inline(always)]
    pub fn find(&mut self, x: usize) -> usize {
        let mut current = x;
        while self.parent[current] >= 0 {
            current = self.parent[current] as usize;
        }
        let root = current;
        let mut current = x;
        while self.parent[current] >= 0 {
            let parent = self.parent[current] as usize;
            self.parent[current] = root as i32;
            current = parent;
        }
        root
    }

    #[inline(always)]
    fn union_roots_unchecked(&mut self, mut root_a: usize, mut root_b: usize) {
        debug_assert!(root_a != root_b);

        // Union by size: attach smaller component to larger (more negative) root.
        if self.parent[root_a] < self.parent[root_b] {
            mem::swap(&mut root_a, &mut root_b);
        }

        self.parent[root_b] += self.parent[root_a];
        self.parent[root_a] = root_b as i32;

        self.min_start[root_b] = self.min_start[root_b].min(self.min_start[root_a]);
        self.max_start[root_b] = self.max_start[root_b].max(self.max_start[root_a]);
        self.min_end[root_b] = self.min_end[root_b].min(self.min_end[root_a]);
        self.max_end[root_b] = self.max_end[root_b].max(self.max_end[root_a]);

        let mask_words = self.mask_words;
        let a_off = root_a * mask_words;
        let b_off = root_b * mask_words;
        for word in 0..mask_words {
            self.sample_masks[b_off + word] |= self.sample_masks[a_off + word];
        }
    }

    #[inline(always)]
    pub fn union_unchecked(&mut self, a: usize, b: usize) {
        let root_a = self.find(a);
        let root_b = self.find(b);
        if root_a != root_b {
            self.union_roots_unchecked(root_a, root_b);
        }
    }

    fn merged_bbox_diameter(&self, root_a: usize, root_b: usize) -> f64 {
        let min_start = self.min_start[root_a].min(self.min_start[root_b]);
        let max_start = self.max_start[root_a].max(self.max_start[root_b]);
        let min_end = self.min_end[root_a].min(self.min_end[root_b]);
        let max_end = self.max_end[root_a].max(self.max_end[root_b]);

        let d_start = max_start - min_start;
        let d_end = max_end - min_end;
        (d_start * d_start + d_end * d_end).sqrt()
    }

    #[inline(always)]
    pub fn try_union(&mut self, a: usize, b: usize) -> bool {
        self.try_union_with_diameter(a, b, None)
    }

    #[inline(always)]
    pub fn try_union_with_diameter(
        &mut self,
        a: usize,
        b: usize,
        max_diameter: Option<f64>,
    ) -> bool {
        let root_a = self.find(a);
        let root_b = self.find(b);
        if root_a == root_b {
            return false;
        }
        if !self.valid_edge(root_a, root_b) {
            return false;
        }
        if let Some(max_diameter) = max_diameter {
            if self.merged_bbox_diameter(root_a, root_b) > max_diameter + DISTANCE_OFFSET {
                return false;
            }
        }
        self.union_roots_unchecked(root_a, root_b);
        true
    }

    #[inline(always)]
    pub fn can_union(&mut self, a: usize, b: usize) -> bool {
        let root_a = self.find(a);
        let root_b = self.find(b);
        root_a != root_b && self.valid_edge(root_a, root_b)
    }

    #[inline(always)]
    fn valid_edge(&self, root_a: usize, root_b: usize) -> bool {
        if self.allow_intrasample {
            return true;
        }
        let mask_words = self.mask_words;
        let a_off = root_a * mask_words;
        let b_off = root_b * mask_words;
        for word in 0..mask_words {
            if (self.sample_masks[a_off + word] & self.sample_masks[b_off + word]) != 0 {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::svtype::SvType;
    use crate::core::variant::test_utils;
    use proptest::prelude::*;
    use rand::seq::SliceRandom;
    use std::collections::HashSet;

    pub fn create_variants() -> Vec<VariantInternal> {
        test_utils::insertion_fixture_variants(SvType::INSERTION)
    }

    fn create_variants_disjoint() -> Vec<VariantInternal> {
        test_utils::insertion_fixture_variants_disjoint(SvType::INSERTION)
    }

    fn create_large_variant_set(num_variants: usize, max_sample_id: usize) -> Vec<VariantInternal> {
        let mut variants: Vec<VariantInternal> = (0..num_variants)
            .map(|i| {
                test_utils::from_parts(
                    i % max_sample_id,
                    format!("var{}", i),
                    SvType::INSERTION,
                    1.0,
                    1.0,
                )
                .unwrap()
            })
            .collect();

        let mut rng = rand::rng();
        variants.shuffle(&mut rng);

        variants
    }

    fn variants_from_sample_ids(sample_ids: &[usize]) -> Vec<VariantInternal> {
        sample_ids
            .iter()
            .enumerate()
            .map(|(i, &sample_id)| {
                test_utils::from_parts(sample_id, format!("var{}", i), SvType::INSERTION, 1.0, 1.0)
                    .unwrap()
            })
            .collect()
    }

    fn boundary_sample_id_strategy() -> impl Strategy<Value = usize> {
        prop_oneof![
            Just(0usize),
            Just(SAMPLES_PER_MASK.saturating_sub(1)),
            Just(SAMPLES_PER_MASK),
            Just(SAMPLES_PER_MASK * 2 - 1),
            Just(SAMPLES_PER_MASK * 2),
            Just(SAMPLES_PER_MASK * 64),
        ]
    }

    fn sample_id_strategy() -> impl Strategy<Value = usize> {
        prop_oneof![
            8 => 0usize..(SAMPLES_PER_MASK * 10),
            2 => boundary_sample_id_strategy(),
        ]
    }

    fn nonempty_sample_ids_strategy() -> impl Strategy<Value = Vec<usize>> {
        (
            proptest::collection::vec(sample_id_strategy(), 0..64),
            boundary_sample_id_strategy(),
        )
            .prop_map(|(mut sample_ids, boundary)| {
                sample_ids.push(boundary);
                sample_ids
            })
    }

    fn union_case_strategy() -> impl Strategy<Value = (Vec<usize>, Vec<(usize, usize)>)> {
        nonempty_sample_ids_strategy().prop_flat_map(|sample_ids| {
            let n = sample_ids.len();
            let ops = proptest::collection::vec((0..n, 0..n), 0..256);
            (Just(sample_ids), ops)
        })
    }

    #[derive(Debug)]
    struct ModelForest {
        parent: Vec<usize>,
        size: Vec<usize>,
        samples: Vec<HashSet<usize>>,
        allow_intrasample: bool,
    }

    impl ModelForest {
        fn new(sample_ids: &[usize], allow_intrasample: bool) -> Self {
            let n = sample_ids.len();
            let parent = (0..n).collect();
            let size = vec![1usize; n];
            let samples = sample_ids
                .iter()
                .map(|&sample_id| HashSet::from([sample_id]))
                .collect();
            Self {
                parent,
                size,
                samples,
                allow_intrasample,
            }
        }

        fn find(&mut self, x: usize) -> usize {
            let mut root = x;
            while self.parent[root] != root {
                root = self.parent[root];
            }

            let mut current = x;
            while self.parent[current] != current {
                let next = self.parent[current];
                self.parent[current] = root;
                current = next;
            }
            root
        }

        fn try_union(&mut self, a: usize, b: usize) -> bool {
            let mut root_a = self.find(a);
            let mut root_b = self.find(b);
            if root_a == root_b {
                return false;
            }

            if !self.allow_intrasample && !self.samples[root_a].is_disjoint(&self.samples[root_b]) {
                return false;
            }

            // Union by size.
            if self.size[root_a] < self.size[root_b] {
                mem::swap(&mut root_a, &mut root_b);
            }

            self.parent[root_b] = root_a;
            self.size[root_a] += self.size[root_b];

            let merged = mem::take(&mut self.samples[root_b]);
            self.samples[root_a].extend(merged);
            true
        }
    }

    proptest! {
        #[test]
        fn prop_forest_init_sets_expected_bits(sample_ids in nonempty_sample_ids_strategy()) {
            let variants = variants_from_sample_ids(&sample_ids);
            let f = Forest::new(&variants, false);

            let n = sample_ids.len();
            prop_assert_eq!(f.parent, vec![-1; n]);

            let max_sample_id = sample_ids.iter().copied().max().unwrap_or(0);
            let expected_mask_words = max_sample_id / SAMPLES_PER_MASK + 1;
            prop_assert_eq!(f.mask_words, expected_mask_words);
            prop_assert_eq!(f.sample_masks.len(), expected_mask_words * n);

            for (i, &sample_id) in sample_ids.iter().enumerate() {
                let mask_id = sample_id / SAMPLES_PER_MASK;
                let bit_position = sample_id % SAMPLES_PER_MASK;
                for m in 0..expected_mask_words {
                    let expected = if m == mask_id {
                        1u64 << bit_position
                    } else {
                        0u64
                    };
                    prop_assert_eq!(f.sample_masks[i * expected_mask_words + m], expected);
                }
            }
        }

        #[test]
        fn prop_forest_union_find_invariants(
            allow_intrasample in any::<bool>(),
            (sample_ids, ops) in union_case_strategy(),
        ) {
            let variants = variants_from_sample_ids(&sample_ids);
            let mut f = Forest::new(&variants, allow_intrasample);

            for (a, b) in ops {
                if a == b {
                    continue;
                }
                let _ = f.try_union(a, b);
            }

            let n = sample_ids.len();
            let mask_words = f.mask_words;

            // Force path compression and collect each node's root.
            let mut roots = vec![0usize; n];
            for (i, root) in roots.iter_mut().enumerate() {
                *root = f.find(i);
            }

            // Compute expected component sizes and sample masks from the final partition.
            let mut counts = vec![0usize; n];
            let mut expected_masks = vec![0u64; n * mask_words];
            for (i, &sample_id) in sample_ids.iter().enumerate() {
                let root = roots[i];
                counts[root] += 1;
                let mask_id = sample_id / SAMPLES_PER_MASK;
                let bit_position = sample_id % SAMPLES_PER_MASK;
                expected_masks[root * mask_words + mask_id] |= 1u64 << bit_position;
            }

            // Roots store negative size; sample masks are only guaranteed to be current at roots.
            for (i, parent) in f.parent.iter().enumerate() {
                if *parent < 0 {
                    prop_assert_eq!((-*parent) as usize, counts[i]);
                    let off = i * mask_words;
                    for word in 0..mask_words {
                        prop_assert_eq!(
                            f.sample_masks[off + word],
                            expected_masks[off + word]
                        );
                    }

                    if !allow_intrasample {
                        let distinct_samples: usize = expected_masks[off..off + mask_words]
                            .iter()
                            .map(|m| m.count_ones() as usize)
                            .sum();
                        prop_assert_eq!(distinct_samples, counts[i]);
                    }
                } else {
                    prop_assert!((*parent as usize) < n);
                }
            }

            // After one `find` per node, all non-roots should point directly at their root.
            for (i, &root) in roots.iter().enumerate() {
                if i != root {
                    prop_assert_eq!(f.parent[i], root as i32);
                }
                prop_assert!(f.parent[root] < 0);
            }

            // `can_union` should match the intended semantics for all pairs.
            for a in 0..n {
                for b in 0..n {
                    if a == b {
                        continue;
                    }
                    let ra = roots[a];
                    let rb = roots[b];
                    let expected = if ra == rb {
                        false
                    } else if allow_intrasample {
                        true
                    } else {
                        let a_off = ra * mask_words;
                        let b_off = rb * mask_words;
                        (0..mask_words).all(|word| {
                            (expected_masks[a_off + word] & expected_masks[b_off + word]) == 0
                        })
                    };

                    prop_assert_eq!(f.can_union(a, b), expected);
                    prop_assert_eq!(f.can_union(b, a), expected);
                }
            }
        }

        #[test]
        fn prop_forest_try_union_return_matches_model(
            allow_intrasample in any::<bool>(),
            (sample_ids, ops) in union_case_strategy(),
        ) {
            let variants = variants_from_sample_ids(&sample_ids);
            let mut f = Forest::new(&variants, allow_intrasample);
            let mut model = ModelForest::new(&sample_ids, allow_intrasample);

            for (a, b) in ops {
                let expected = model.try_union(a, b);
                let got = f.try_union(a, b);
                prop_assert_eq!(got, expected);
            }

            let n = sample_ids.len();
            let mut f_roots = vec![0usize; n];
            for (i, root) in f_roots.iter_mut().enumerate() {
                *root = f.find(i);
            }

            let mut model_roots = vec![0usize; n];
            for (i, root) in model_roots.iter_mut().enumerate() {
                *root = model.find(i);
            }

            for i in 0..n {
                for j in 0..n {
                    prop_assert_eq!(f_roots[i] == f_roots[j], model_roots[i] == model_roots[j]);
                }
            }
        }
    }

    #[test]
    fn test_forest_init_basic() {
        let variants = create_variants();
        let f = Forest::new(&variants, false);

        assert_eq!(f.parent.len(), variants.len());
        assert_eq!(f.parent, vec![-1; variants.len()]);

        let max_sample_id = 4;
        let expected_mask_words = max_sample_id / SAMPLES_PER_MASK + 1;
        assert_eq!(f.mask_words, expected_mask_words);

        assert_eq!(f.sample_masks[0], 0b1);
        assert_eq!(f.sample_masks[expected_mask_words], 0b1);
        assert_eq!(f.sample_masks[3 * expected_mask_words], 0b10);
        assert_eq!(f.sample_masks[7 * expected_mask_words], 0b100);
        assert_eq!(f.sample_masks[11 * expected_mask_words], 0b1000);
        assert_eq!(f.sample_masks[12 * expected_mask_words], 0b10000);
    }

    #[test]
    fn test_forest_init_empty() {
        let variants = vec![];
        let f = Forest::new(&variants, false);
        assert_eq!(f.parent.len(), 0);
        assert_eq!(f.sample_masks.len(), 0);
    }

    #[test]
    fn test_find_basic() {
        let variants = create_variants();
        let mut f = Forest::new(&variants, false);
        assert_eq!(f.find(0), 0);
        assert_eq!(f.find(1), 1);
        assert_eq!(f.find(2), 2);
        assert_eq!(f.find(3), 3);
    }

    #[test]
    #[should_panic]
    fn test_find_with_empty_data() {
        let variants: Vec<VariantInternal> = vec![];
        let mut f = Forest::new(&variants, false);
        f.find(0);
    }

    #[test]
    fn test_union_two_distinct_components() {
        let variants = vec![
            test_utils::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0).unwrap(),
            test_utils::from_parts(1, "var2".to_string(), SvType::INSERTION, 15.0, 10.0).unwrap(),
        ];
        let mut f = Forest::new(&variants, false);
        assert!(f.can_union(0, 1));
        assert!(f.can_union(1, 0));
        f.union_unchecked(0, 1);
        assert!(!f.can_union(0, 1));
        assert!(!f.can_union(1, 0));
        assert_eq!(f.find(0), 1);
        assert_eq!(f.find(1), 1);

        assert_eq!(f.sample_masks, vec![1, 3]);
        assert_eq!(f.parent, vec![1, -2]);
    }

    #[test]
    fn test_forest_union_find() {
        let variants = create_variants();
        let mut f = Forest::new(&variants, false);

        f.union_unchecked(0, 1);
        assert_eq!(f.find(0), f.find(1));

        f.union_unchecked(1, 2);
        assert_eq!(f.find(0), f.find(2));
        assert_eq!(f.find(1), f.find(2));

        f.union_unchecked(3, 4);
        assert_eq!(f.find(3), f.find(4));

        f.union_unchecked(0, 3);
        assert_eq!(f.find(0), f.find(4));

        assert_eq!(f.sample_masks, vec![1, 3, 1, 2, 2, 2, 2, 4, 4, 4, 4, 8, 16]);

        assert_eq!(
            f.parent,
            vec![1, -5, 1, 4, 1, -1, -1, -1, -1, -1, -1, -1, -1]
        );

        // Connect everything where allowed by sample ids
        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                let _ = f.try_union(i, j);
            }
        }

        assert_eq!(
            f.sample_masks,
            vec![1, 31, 1, 2, 2, 2, 2, 4, 6, 6, 4, 8, 16]
        );

        assert_eq!(f.parent, vec![1, -8, 1, 1, 1, 8, 9, 1, -2, -2, -1, 1, 1]);

        // Connect everything regardless of sample_id
        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                f.union_unchecked(i, j);
            }
        }

        assert_eq!(f.parent, vec![1, -13, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
    }

    #[test]
    fn test_forest_union_find_large() {
        let variants = create_large_variant_set(1000, 200);
        let mut f = Forest::new(&variants, false);

        f.union_unchecked(0, 1);
        assert_eq!(f.find(0), f.find(1));

        f.union_unchecked(1, 2);
        assert_eq!(f.find(0), f.find(2));
        assert_eq!(f.find(1), f.find(2));

        f.union_unchecked(3, 4);
        assert_eq!(f.find(3), f.find(4));

        f.union_unchecked(0, 900);
        assert_eq!(f.find(0), f.find(900));

        f.union_unchecked(500, 501);
        assert_eq!(f.find(500), f.find(501));
    }

    #[test]
    fn test_forest_can_union_edge_cases() {
        // Disjoint is used here since can_union will check valid edges (reject intra-sample union)
        let variants = create_variants_disjoint();
        let mut f = Forest::new(&variants, false);

        // Any two variants can be unioned
        assert!(f.can_union(0, 1));
        f.union_unchecked(0, 1);
        assert!(!f.can_union(0, 1)); // Already unioned, cannot union again
        f.union_unchecked(1, 2);
        assert!(!f.can_union(0, 2)); // All three are connected
        f.union_unchecked(3, 4);
        assert!(f.can_union(0, 3)); // 0 and 3 are still separate, can be unioned
        f.union_unchecked(0, 3);
        assert!(!f.can_union(0, 4)); // All are now connected
    }

    #[test]
    fn test_forest_invalid_union_cases() {
        let variants = create_variants_disjoint();
        let mut f = Forest::new(&variants, false);

        f.union_unchecked(0, 1);
        f.union_unchecked(2, 3);
        assert!(f.can_union(0, 2));

        // Modify the sample mask directly to simulate an invalid state
        let i = f.find(2);
        let j = f.find(0);
        let mask_words = f.mask_words;
        f.sample_masks[i * mask_words] |= f.sample_masks[j * mask_words];
        // Invalid union due to overlapping masks
        assert!(!f.can_union(0, 2));
    }

    #[test]
    fn test_forest_path_compression() {
        let variants = create_large_variant_set(100, 10);
        let mut f = Forest::new(&variants, false);

        for i in 1..variants.len() {
            f.union_unchecked(0, i);
        }

        for i in 1..variants.len() {
            assert_eq!(f.find(i), f.find(0));
        }
    }

    #[test]
    fn test_forest_sample_mask_integrity() {
        let variants = create_large_variant_set(63, 63);
        let mut f = Forest::new(&variants, false);

        f.union_unchecked(0, 1);
        f.union_unchecked(2, 3);
        f.union_unchecked(0, 2);

        let root = f.find(0);
        let mask_words = f.mask_words;
        let expected_mask = f.sample_masks[root * mask_words];

        for i in 0..4 {
            let x = f.find(i);
            assert_eq!(f.sample_masks[x * mask_words], expected_mask);
        }
    }

    #[test]
    fn test_forest_allow_intrasample_true() {
        let variants = vec![
            test_utils::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0).unwrap(),
            test_utils::from_parts(0, "var2".to_string(), SvType::INSERTION, 20.0, 10.0).unwrap(),
        ];

        let mut f_allow_true = Forest::new(&variants, true);
        assert!(
            f_allow_true.can_union(0, 1),
            "Should be able to union variants from the same sample when allow_intrasample is true"
        );

        let mut f_allow_false = Forest::new(&variants, false);
        assert!(
            !f_allow_false.can_union(0, 1),
            "Should NOT be able to union variants from the same sample when allow_intrasample is false"
        );
    }
}
