use crate::core::variant_internal::VariantInternal;
use std::mem;

// Each component may require multiple 64-bit integers to hold its bitset of sample IDs if there are many samples
// parent[i] - is negative if root, more negative means bigger set; if nonnegative, indicates the parent
// sample_masks - For each root node, a bitmask of which samples are present in its component
const SAMPLES_PER_MASK: usize = 63;

#[derive(Debug)]
pub struct Forest {
    pub parent: Vec<i32>,
    pub sample_masks: Vec<Vec<u64>>,
    allow_intrasample: bool,
}

impl Forest {
    pub fn new(variants: &[VariantInternal], allow_intrasample: bool) -> Self {
        let n = variants.len();
        let max_sample = variants.iter().map(|v| v.sample_id).max().unwrap_or(0);
        let masks_needed = if n == 0 {
            0
        } else {
            max_sample / SAMPLES_PER_MASK + 1
        };
        let parent = vec![-1; n];
        let mut sample_masks = vec![vec![0; n]; masks_needed];
        for (i, variant) in variants.iter().enumerate() {
            let mask_id = variant.sample_id / SAMPLES_PER_MASK;
            let bit_position = variant.sample_id % SAMPLES_PER_MASK;
            sample_masks[mask_id][i] |= 1u64 << bit_position;
        }
        Forest {
            parent,
            sample_masks,
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
    pub fn union(&mut self, a: usize, b: usize) {
        let mut root_a = self.find(a);
        let mut root_b = self.find(b);
        if root_a != root_b {
            if self.parent[root_a] < self.parent[root_b] {
                mem::swap(&mut root_a, &mut root_b);
            }
            self.parent[root_b] += self.parent[root_a];
            self.parent[root_a] = root_b as i32;
            for masks in &mut self.sample_masks {
                masks[root_b] |= masks[root_a];
            }
        }
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
        self.sample_masks
            .iter()
            .all(|mask| (mask[root_a] & mask[root_b]) == 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::svtype::SvType;
    use rand::{rngs::StdRng, seq::SliceRandom, Rng, SeedableRng};
    use std::time::Instant;

    pub fn create_variants() -> Vec<VariantInternal> {
        vec![
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(0, "var2".to_string(), SvType::INSERTION, 1.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(0, "var3".to_string(), SvType::INSERTION, 18.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(1, "var4".to_string(), SvType::INSERTION, 12.0, 7.0)
                .unwrap(),
            VariantInternal::from_parts(1, "var5".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(1, "var6".to_string(), SvType::INSERTION, 30.0, 30.0)
                .unwrap(),
            VariantInternal::from_parts(1, "var7".to_string(), SvType::INSERTION, 0.0, 0.0)
                .unwrap(),
            VariantInternal::from_parts(2, "var8".to_string(), SvType::INSERTION, 12.0, 12.0)
                .unwrap(),
            VariantInternal::from_parts(2, "var9".to_string(), SvType::INSERTION, 15.0, 15.0)
                .unwrap(),
            VariantInternal::from_parts(2, "var10".to_string(), SvType::INSERTION, 20.0, 20.0)
                .unwrap(),
            VariantInternal::from_parts(2, "var11".to_string(), SvType::INSERTION, 28.0, 28.0)
                .unwrap(),
            VariantInternal::from_parts(3, "var12".to_string(), SvType::INSERTION, 25.0, 25.0)
                .unwrap(),
            VariantInternal::from_parts(4, "var13".to_string(), SvType::INSERTION, 22.0, 22.0)
                .unwrap(),
        ]
    }

    fn create_variants_disjoint() -> Vec<VariantInternal> {
        vec![
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(1, "var2".to_string(), SvType::INSERTION, 1.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(2, "var3".to_string(), SvType::INSERTION, 18.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(3, "var4".to_string(), SvType::INSERTION, 12.0, 7.0)
                .unwrap(),
            VariantInternal::from_parts(4, "var5".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(5, "var6".to_string(), SvType::INSERTION, 30.0, 30.0)
                .unwrap(),
            VariantInternal::from_parts(6, "var7".to_string(), SvType::INSERTION, 0.0, 0.0)
                .unwrap(),
            VariantInternal::from_parts(7, "var8".to_string(), SvType::INSERTION, 12.0, 12.0)
                .unwrap(),
            VariantInternal::from_parts(8, "var9".to_string(), SvType::INSERTION, 15.0, 15.0)
                .unwrap(),
            VariantInternal::from_parts(9, "var10".to_string(), SvType::INSERTION, 20.0, 20.0)
                .unwrap(),
            VariantInternal::from_parts(10, "var11".to_string(), SvType::INSERTION, 28.0, 28.0)
                .unwrap(),
            VariantInternal::from_parts(11, "var12".to_string(), SvType::INSERTION, 25.0, 25.0)
                .unwrap(),
            VariantInternal::from_parts(12, "var13".to_string(), SvType::INSERTION, 22.0, 22.0)
                .unwrap(),
        ]
    }

    fn create_large_variant_set(num_variants: usize, max_sample_id: usize) -> Vec<VariantInternal> {
        let mut variants: Vec<VariantInternal> = (0..num_variants)
            .map(|i| {
                VariantInternal::from_parts(
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

    fn create_variant_set(
        num_variants: usize,
        max_sample_id: usize,
        seed: u64,
    ) -> Vec<VariantInternal> {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut variants: Vec<VariantInternal> = (0..num_variants)
            .map(|i| {
                VariantInternal::from_parts(
                    rng.random_range(0..max_sample_id),
                    format!("var{}", i),
                    SvType::INSERTION,
                    1.0,
                    1.0,
                )
                .unwrap()
            })
            .collect();

        variants.shuffle(&mut rng);

        variants
    }

    #[test]
    fn test_forest_init_basic() {
        let variants = create_variants();
        let f = Forest::new(&variants, false);

        assert_eq!(f.parent.len(), variants.len());
        assert_eq!(f.parent, vec![-1; variants.len()]);

        let max_sample_id = 4;
        let masks_needed = max_sample_id / SAMPLES_PER_MASK + 1;
        assert_eq!(f.sample_masks.len(), masks_needed);

        assert_eq!(f.sample_masks[0][0], 0b1);
        assert_eq!(f.sample_masks[0][1], 0b1);
        assert_eq!(f.sample_masks[0][3], 0b10);
        assert_eq!(f.sample_masks[0][7], 0b100);
        assert_eq!(f.sample_masks[0][11], 0b1000);
        assert_eq!(f.sample_masks[0][12], 0b10000);
    }

    #[test]
    fn test_forest_init_empty() {
        let variants = vec![];
        let f = Forest::new(&variants, false);
        assert_eq!(f.parent.len(), 0);
        assert_eq!(f.sample_masks.len(), 0);
    }

    #[test]
    fn test_forest_init_single_variant() {
        let variant =
            vec![
                VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                    .unwrap(),
            ];
        let f = Forest::new(&variant, false);

        assert_eq!(f.parent.len(), 1);
        assert_eq!(f.parent, vec![-1; 1]);
        assert_eq!(f.sample_masks.len(), 1);
        assert_eq!(f.sample_masks[0][0], 0b1);
    }

    #[test]
    fn test_forest_init_same_sample_id() {
        let variants = vec![
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(0, "var2".to_string(), SvType::INSERTION, 20.0, 10.0)
                .unwrap(),
        ];
        let f = Forest::new(&variants, false);

        assert_eq!(f.parent.len(), 2);
        assert_eq!(f.parent, vec![-1; 2]);
        assert_eq!(f.sample_masks.len(), 1);
        assert_eq!(f.sample_masks[0][0], 0b1);
        assert_eq!(f.sample_masks[0][1], 0b1);
    }

    #[test]
    fn test_forest_init_max_sample_id() {
        let sample_id = SAMPLES_PER_MASK * 64;
        let variant = vec![VariantInternal::from_parts(
            sample_id,
            "var1".to_string(),
            SvType::INSERTION,
            10.0,
            5.0,
        )
        .unwrap()];
        let f = Forest::new(&variant, false);

        assert_eq!(f.parent.len(), 1);
        assert_eq!(f.parent, vec![-1; 1]);
        assert_eq!(f.sample_masks.len(), 65);
        assert_eq!(f.sample_masks[64][0], 0b1);
    }

    #[test]
    fn test_forest_init_boundary_condition() {
        let sample_id = SAMPLES_PER_MASK - 1;
        let variant = vec![VariantInternal::from_parts(
            sample_id,
            "var1".to_string(),
            SvType::INSERTION,
            10.0,
            5.0,
        )
        .unwrap()];
        let f = Forest::new(&variant, false);

        assert_eq!(f.parent.len(), 1);
        assert_eq!(f.parent, vec![-1; 1]);
        assert_eq!(f.sample_masks.len(), 1);
        assert_eq!(f.sample_masks[0][0], 1 << (SAMPLES_PER_MASK - 1));
    }

    #[test]
    fn test_forest_large_init() {
        let num_variants = 1000;
        let max_sample_id = 512;
        let variants = create_large_variant_set(num_variants, max_sample_id);
        let f = Forest::new(&variants, false);
        assert_eq!(f.parent.len(), num_variants);
        let masks_needed = max_sample_id / SAMPLES_PER_MASK + 1;
        assert_eq!(f.sample_masks.len(), masks_needed);
        for (i, variant) in variants.iter().enumerate() {
            let mask_id = variant.sample_id / SAMPLES_PER_MASK;
            let bit_position = variant.sample_id % SAMPLES_PER_MASK;
            assert_eq!(f.sample_masks[mask_id][i], 1 << bit_position);
        }
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
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(1, "var2".to_string(), SvType::INSERTION, 15.0, 10.0)
                .unwrap(),
        ];
        let mut f = Forest::new(&variants, false);
        assert!(f.can_union(0, 1));
        assert!(f.can_union(1, 0));
        f.union(0, 1);
        assert!(!f.can_union(0, 1));
        assert!(!f.can_union(1, 0));
        assert_eq!(f.find(0), 1);
        assert_eq!(f.find(1), 1);

        assert_eq!(f.sample_masks, vec![vec![1, 3]]);
        assert_eq!(f.parent, vec![1, -2]);
    }

    #[test]
    fn test_union_boundary_condition() {
        let sample_id = SAMPLES_PER_MASK - 1;
        let variant = vec![
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(
                sample_id,
                "var2".to_string(),
                SvType::INSERTION,
                10.0,
                5.0,
            )
            .unwrap(),
        ];
        let mut f = Forest::new(&variant, false);
        assert!(f.can_union(0, 1));
        assert!(f.can_union(1, 0));
        f.union(0, 1);
        assert!(!f.can_union(0, 1));
        assert!(!f.can_union(1, 0));
        assert_eq!(f.find(0), 1);
        assert_eq!(f.find(1), 1);
    }

    #[test]
    fn test_union_max_samples() {
        let sample_id = SAMPLES_PER_MASK * 64;
        let variant = vec![
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(
                sample_id,
                "var2".to_string(),
                SvType::INSERTION,
                10.0,
                5.0,
            )
            .unwrap(),
        ];
        let mut f = Forest::new(&variant, false);
        assert!(f.can_union(1, 0));
        assert!(f.can_union(0, 1));
        f.union(1, 0);
        assert!(!f.can_union(1, 0));
        assert!(!f.can_union(0, 1));
        assert_eq!(f.find(0), 0);
        assert_eq!(f.find(1), 0);
    }

    #[test]
    fn test_forest_union_find() {
        let variants = create_variants();
        let mut f = Forest::new(&variants, false);

        f.union(0, 1);
        assert_eq!(f.find(0), f.find(1));

        f.union(1, 2);
        assert_eq!(f.find(0), f.find(2));
        assert_eq!(f.find(1), f.find(2));

        f.union(3, 4);
        assert_eq!(f.find(3), f.find(4));

        f.union(0, 3);
        assert_eq!(f.find(0), f.find(4));

        assert_eq!(
            f.sample_masks,
            vec![vec![1, 3, 1, 2, 2, 2, 2, 4, 4, 4, 4, 8, 16]]
        );

        assert_eq!(
            f.parent,
            vec![1, -5, 1, 4, 1, -1, -1, -1, -1, -1, -1, -1, -1]
        );

        // Connect everything where allowed by sample ids
        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                if f.can_union(i, j) {
                    f.union(i, j);
                }
            }
        }

        assert_eq!(
            f.sample_masks,
            vec![vec![1, 31, 1, 2, 2, 2, 2, 4, 6, 6, 4, 8, 16]]
        );

        assert_eq!(f.parent, vec![1, -8, 1, 1, 1, 8, 9, 1, -2, -2, -1, 1, 1]);

        // Connect everything regardless of sample_id
        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                f.union(i, j);
            }
        }

        assert_eq!(f.parent, vec![1, -13, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
    }

    #[test]
    fn test_forest_union_find_random_sample_ids() {
        let variants = create_variant_set(100, 10, 42);

        let mut f = Forest::new(&variants, false);

        assert!(f.can_union(0, 1));
        f.union(0, 1);
        assert_eq!(f.find(0), f.find(1));

        for i in (10..90).skip(2) {
            if f.can_union(i, i + 2) {
                f.union(i, i + 2);
            }
        }

        assert_eq!(
            f.sample_masks,
            vec![vec![
                64, 65, 64, 256, 256, 8, 16, 1, 4, 32, 16, 8, 32, 16, 48, 528, 16, 512, 146, 513,
                2, 1, 2, 901, 2, 256, 1015, 4, 512, 512, 64, 512, 128, 866, 256, 2, 1, 64, 32, 32,
                16, 2, 64, 262, 594, 256, 16, 256, 512, 257, 16, 256, 560, 377, 512, 1, 512, 8,
                820, 64, 256, 16, 32, 1, 16, 11, 256, 8, 352, 8, 32, 153, 256, 128, 352, 1, 32, 1,
                32, 27, 96, 16, 32, 2, 313, 2, 16, 262, 8, 4, 1, 256, 16, 64, 4, 512, 32, 256, 8,
                16,
            ]]
        );

        assert_eq!(
            f.parent,
            vec![
                1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 14, 15, -2, -2, 18, 19, -3, -2, 18,
                23, -1, -5, 26, 23, -9, 23, 26, 23, 26, 33, 26, -5, 26, 33, 26, 33, 26, 33, 26, 43,
                44, -3, -4, 43, 44, 49, 44, -2, 52, 53, -3, -6, 52, 53, 58, 53, -5, 53, 58, 53, 58,
                65, 58, -3, 68, 65, -3, 71, 68, -4, 74, 71, -3, 71, 74, 79, 80, -4, -2, 79, 84, 79,
                -5, 87, 84, -3, 84, 87, 84, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            ]
        );

        for i in (1..97).skip(3) {
            if f.can_union(i, i + 3) {
                f.union(i, i + 3);
            }
        }

        assert_eq!(
            f.sample_masks,
            vec![vec![
                64, 65, 64, 256, 256, 8, 16, 273, 12, 48, 16, 8, 32, 16, 571, 528, 16, 512, 146,
                513, 2, 1, 2, 901, 2, 256, 1015, 4, 512, 512, 64, 512, 128, 866, 256, 2, 1, 64, 32,
                32, 16, 2, 64, 262, 851, 256, 16, 256, 512, 257, 16, 256, 560, 377, 512, 1, 512, 8,
                831, 64, 256, 16, 32, 1, 16, 11, 256, 8, 352, 8, 32, 505, 256, 128, 352, 1, 32, 1,
                32, 379, 96, 16, 32, 2, 377, 2, 16, 798, 8, 4, 1, 256, 16, 64, 260, 512, 32, 256,
                8, 48,
            ]]
        );

        assert_eq!(
            f.parent,
            vec![
                1, -2, -1, -1, 7, 8, 9, -3, -2, -2, 7, 14, 14, 15, -6, -2, 18, 14, -3, 14, 18, 23,
                14, -5, 26, 23, -9, 23, 26, 23, 26, 33, 26, -5, 26, 33, 26, 33, 26, 33, 26, 43, 44,
                -3, -6, 43, 44, 44, 44, 44, 52, 53, -3, -6, 52, 53, 58, 53, -8, 53, 58, 53, 58, 58,
                58, 58, 68, 58, 71, 71, 71, -7, 74, 71, 79, 71, 79, 79, 80, -7, -2, 79, 84, 79, -6,
                87, 84, -6, 84, 87, 84, 94, 87, 84, -2, 87, 99, -1, 87, -2,
            ]
        );

        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                if f.can_union(i, j) {
                    f.union(i, j);
                }
            }
        }

        assert_eq!(
            f.sample_masks,
            vec![vec![
                64, 381, 64, 256, 256, 8, 16, 337, 12, 48, 16, 8, 32, 16, 831, 528, 16, 512, 402,
                513, 2, 1, 2, 949, 2, 256, 1015, 4, 512, 512, 64, 512, 128, 866, 256, 2, 1, 64, 32,
                32, 16, 2, 64, 886, 851, 256, 16, 256, 512, 257, 16, 256, 560, 377, 512, 1, 512, 8,
                831, 64, 256, 16, 32, 1, 16, 11, 256, 8, 352, 8, 32, 505, 256, 128, 352, 1, 32, 1,
                32, 379, 96, 16, 32, 2, 377, 2, 16, 798, 8, 4, 1, 256, 16, 64, 260, 512, 32, 256,
                8, 48,
            ]]
        );

        assert_eq!(
            f.parent,
            vec![
                1, -7, 7, 1, 7, 1, 1, -4, 1, 1, 7, 14, 14, 43, -8, 43, 18, 14, -4, 14, 18, 23, 14,
                -7, 26, 23, -9, 23, 26, 23, 26, 33, 26, -5, 26, 33, 26, 33, 26, 33, 26, 43, 44, -7,
                -6, 43, 44, 44, 44, 44, 52, 53, -3, -6, 52, 53, 58, 53, -8, 53, 58, 53, 58, 58, 58,
                58, 71, 58, 71, 71, 71, -7, 79, 71, 79, 71, 79, 79, 43, -7, 43, 79, 84, 79, -6, 87,
                84, -6, 84, 87, 84, 14, 87, 84, 14, 87, 23, 18, 87, 23,
            ]
        );

        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                f.union(i, j);
            }
        }

        assert_eq!(
            f.parent,
            vec![
                1, -100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            ]
        );
    }

    #[test]
    fn test_forest_union_find_random_sample_ids_large_sample_ids() {
        let variants = create_variant_set(100, SAMPLES_PER_MASK * 2, 42);
        let mut f = Forest::new(&variants, false);

        for i in 0..variants.len() {
            for j in (i + 1)..variants.len() {
                if f.can_union(i, j) {
                    f.union(i, j);
                }
            }
        }

        assert_eq!(
            f.sample_masks,
            vec![
                vec![
                    0,
                    8592754029918624219,
                    0,
                    0,
                    0,
                    140737488355328,
                    4611686018427387904,
                    2048,
                    2147483648,
                    0,
                    288230376151711744,
                    17592186044416,
                    0,
                    144115188075855872,
                    4503599627370496,
                    0,
                    72057594037927936,
                    0,
                    0,
                    256,
                    131072,
                    64,
                    8388608,
                    0,
                    16777216,
                    0,
                    2147483648,
                    33554432,
                    0,
                    4914007737883623441,
                    0,
                    0,
                    0,
                    0,
                    0,
                    131072,
                    8,
                    0,
                    0,
                    0,
                    2251799813685248,
                    1048576,
                    0,
                    8589934592,
                    524288,
                    0,
                    1152921504606846976,
                    0,
                    0,
                    2,
                    9007199254740992,
                    0,
                    0,
                    0,
                    0,
                    1,
                    0,
                    562949953421312,
                    4294967296,
                    0,
                    0,
                    9007199254740992,
                    0,
                    128,
                    1125899906842624,
                    1048576,
                    0,
                    274877906944,
                    0,
                    17592186044416,
                    0,
                    4503599627370496,
                    0,
                    0,
                    9007199254740992,
                    4096,
                    0,
                    16,
                    0,
                    562949953421312,
                    0,
                    288230376151711744,
                    0,
                    65536,
                    0,
                    65536,
                    2305843009213693952,
                    0,
                    8796093022208,
                    134217728,
                    16,
                    0,
                    4611686018427387904,
                    0,
                    134217728,
                    0,
                    0,
                    0,
                    281474976710656,
                    9007199254740992,
                ],
                vec![
                    131072,
                    4134825489604769325,
                    262144,
                    70368744177664,
                    8796093022208,
                    0,
                    0,
                    0,
                    0,
                    8,
                    0,
                    0,
                    4096,
                    0,
                    0,
                    2305843009213693952,
                    0,
                    576460752303423488,
                    137438953472,
                    0,
                    0,
                    0,
                    0,
                    536870912,
                    0,
                    1099511627776,
                    0,
                    0,
                    18014398509481984,
                    2900689107762417701,
                    65536,
                    2305843009213693952,
                    67108864,
                    70368744177664,
                    274877906944,
                    0,
                    0,
                    16384,
                    512,
                    1,
                    0,
                    0,
                    8192,
                    0,
                    0,
                    17592186044416,
                    0,
                    281474976710656,
                    9007199254740992,
                    0,
                    0,
                    274877906944,
                    4,
                    1024,
                    1152921504606846976,
                    0,
                    18014398509481984,
                    0,
                    0,
                    131072,
                    70368744177664,
                    0,
                    4096,
                    0,
                    0,
                    0,
                    140737488355328,
                    0,
                    262144,
                    0,
                    4,
                    0,
                    17592186044416,
                    137438953472,
                    70368744439812,
                    0,
                    4,
                    0,
                    32,
                    0,
                    1048576,
                    0,
                    1,
                    0,
                    281474976710656,
                    0,
                    0,
                    549755813888,
                    0,
                    0,
                    0,
                    1099511627776,
                    0,
                    524288,
                    0,
                    72057594037927936,
                    32,
                    70368744177664,
                    0,
                    0,
                ],
            ]
        )
    }

    #[test]
    fn test_forest_union_find_large() {
        let variants = create_large_variant_set(1000, 200);
        let mut f = Forest::new(&variants, false);

        f.union(0, 1);
        assert_eq!(f.find(0), f.find(1));

        f.union(1, 2);
        assert_eq!(f.find(0), f.find(2));
        assert_eq!(f.find(1), f.find(2));

        f.union(3, 4);
        assert_eq!(f.find(3), f.find(4));

        f.union(0, 900);
        assert_eq!(f.find(0), f.find(900));

        f.union(500, 501);
        assert_eq!(f.find(500), f.find(501));
    }

    #[test]
    fn test_forest_can_union_edge_cases() {
        // Disjoint is used here since can_union will check valid edges (reject intra-sample union)
        let variants = create_variants_disjoint();
        let mut f = Forest::new(&variants, false);

        // Any two variants can be unioned
        assert!(f.can_union(0, 1));
        f.union(0, 1);
        assert!(!f.can_union(0, 1)); // Already unioned, cannot union again
        f.union(1, 2);
        assert!(!f.can_union(0, 2)); // All three are connected
        f.union(3, 4);
        assert!(f.can_union(0, 3)); // 0 and 3 are still separate, can be unioned
        f.union(0, 3);
        assert!(!f.can_union(0, 4)); // All are now connected
    }

    #[test]
    fn test_forest_invalid_union_cases() {
        let variants = create_variants_disjoint();
        let mut f = Forest::new(&variants, false);

        f.union(0, 1);
        f.union(2, 3);
        assert!(f.can_union(0, 2));

        // Modify the sample mask directly to simulate an invalid state
        let i = f.find(2);
        let j = f.find(0);
        f.sample_masks[0][i] |= f.sample_masks[0][j];
        // Invalid union due to overlapping masks
        assert!(!f.can_union(0, 2));
    }

    #[test]
    fn test_forest_path_compression() {
        let variants = create_large_variant_set(100, 10);
        let mut f = Forest::new(&variants, false);

        for i in 1..variants.len() {
            f.union(0, i);
        }

        for i in 1..variants.len() {
            assert_eq!(f.find(i), f.find(0));
        }
    }

    #[test]
    fn test_forest_sample_mask_integrity() {
        let variants = create_large_variant_set(63, 63);
        let mut f = Forest::new(&variants, false);

        f.union(0, 1);
        f.union(2, 3);
        f.union(0, 2);

        let root = f.find(0);
        let expected_mask = f.sample_masks[0][root];

        for i in 0..4 {
            let x = f.find(i);
            assert_eq!(f.sample_masks[0][x], expected_mask);
        }
    }

    #[cfg(debug_assertions)]
    #[test]
    #[ignore = "Benchmark test"]
    fn test_union_find_timing() {
        let num_variants = 5000;
        let max_sample_id = 1000;
        let variants = create_variant_set(num_variants, max_sample_id, 42);
        let mut f = Forest::new(&variants, false);

        let t0 = Instant::now();
        for i in 0..num_variants / 2 {
            if f.can_union(i, num_variants - 1 - i) {
                f.union(i, num_variants - 1 - i);
            }
        }
        let duration = t0.elapsed();
        println!("Duration: {:?}", duration);
    }

    #[test]
    fn test_forest_allow_intrasample_true() {
        let variants = vec![
            VariantInternal::from_parts(0, "var1".to_string(), SvType::INSERTION, 10.0, 5.0)
                .unwrap(),
            VariantInternal::from_parts(0, "var2".to_string(), SvType::INSERTION, 20.0, 10.0)
                .unwrap(),
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
