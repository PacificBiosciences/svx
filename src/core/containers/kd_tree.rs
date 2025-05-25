use kiddo::{traits::DistanceMetric, NearestNeighbour, SquaredEuclidean};

#[cfg(not(any(feature = "use_mutable_kdtree", feature = "use_immutable_kdtree")))]
compile_error!(
    "Enable at least one of 'use_mutable_kdtree' or 'use_immutable_kdtree'. \
Check Cargo.toml default features or build command. \
If both are enabled, mutable will be preferred."
);

/*
TODO: Possibly normalize the spatial coordinates
TODO: Implement a immutable KD-tree from scratch. The immutable KD-tree of kiddo has some degenerate behaviour when the leaf size is small (I suspect this is the same thing that would otherwise cause a crash in the mutable tree)
*/

// NOTE: The DistanceMetric trait does not allow passing weights... This is a hacky way of doing it.
const WEIGHTS_2D: [f64; 2] = [1.0, 0.5];

pub struct WeightedManhattan<const K: usize> {}

impl<const K: usize> WeightedManhattan<K> {
    pub fn new() -> Self {
        WeightedManhattan {}
    }
}

impl<const K: usize> Default for WeightedManhattan<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl DistanceMetric<f64, 2> for WeightedManhattan<2> {
    fn dist(a: &[f64; 2], b: &[f64; 2]) -> f64 {
        a.iter()
            .zip(b.iter())
            .zip(WEIGHTS_2D.iter())
            .map(|((&a_val, &b_val), &weight)| (a_val - b_val).abs() * weight)
            .sum()
    }

    fn dist1(a: f64, b: f64) -> f64 {
        (a - b).abs()
    }
}

const K: usize = 2;

#[cfg(all(feature = "use_immutable_kdtree", not(feature = "use_mutable_kdtree")))]
mod kdtree_impl {
    use super::*;
    use crate::core::variant_internal::VariantInternal;
    use kiddo::immutable::float::kdtree::ImmutableKdTree;
    use std::num::NonZero;

    const BUCKET_SIZE: usize = 512;
    pub type KdTreeItemType = u64;

    pub struct VariantKdTree {
        kdtree: ImmutableKdTree<f64, KdTreeItemType, K, BUCKET_SIZE>,
        num_variants: usize,
    }

    impl VariantKdTree {
        pub fn new(variants: &[VariantInternal]) -> Self {
            let points: Vec<[f64; K]> = variants.iter().map(|v| v.point()).collect();
            let kdtree = ImmutableKdTree::new_from_slice(&points);
            let num_variants = variants.len();
            assert_eq!(num_variants, kdtree.size());
            Self {
                kdtree,
                num_variants,
            }
        }

        pub fn knn(&self, query: &VariantInternal, k: usize) -> Vec<NearestNeighbour<f64, u64>> {
            self.kdtree
                .nearest_n::<SquaredEuclidean>(&query.point(), NonZero::new(k).unwrap())
        }

        pub fn dump_coordinates(&self, variants: &[VariantInternal]) {
            assert_eq!(
                variants.len(),
                self.num_variants,
                "Provided variants slice length does not match the one used to build the KdTree"
            );
            for (i, variant) in variants.iter().enumerate() {
                let p = variant.point();
                println!(
                    "{},{},{},{:?},{},{},{}",
                    variant.svtype, variant.vcf_id, variant.id, variant.tr_contained, i, p[0], p[1]
                );
            }
        }

        pub fn size(&self) -> usize {
            self.num_variants
        }
    }
}

#[cfg(feature = "use_mutable_kdtree")]
mod kdtree_impl {
    use super::*;
    use crate::core::variant_internal::VariantInternal;
    use kiddo::float::kdtree::KdTree;

    const BUCKET_SIZE: usize = 512;
    pub type KdTreeItemType = u32;

    pub struct VariantKdTree {
        kdtree: KdTree<f64, KdTreeItemType, K, BUCKET_SIZE, u32>,
        num_variants: usize,
    }

    impl VariantKdTree {
        pub fn new(variants: &[VariantInternal]) -> Self {
            let mut kdtree = KdTree::with_capacity(variants.len());
            for (idx, v) in variants.iter().enumerate() {
                let point = v.point();
                kdtree.add(&point, idx as KdTreeItemType);
            }
            let num_variants = variants.len();
            assert_eq!(num_variants, kdtree.size() as usize);
            Self {
                kdtree,
                num_variants,
            }
        }

        pub fn knn(
            &self,
            query: &VariantInternal,
            k: usize,
        ) -> Vec<NearestNeighbour<f64, KdTreeItemType>> {
            self.kdtree.nearest_n::<SquaredEuclidean>(&query.point(), k)
        }

        pub fn dump_coordinates(&self, variants: &[VariantInternal]) {
            assert_eq!(
                variants.len(),
                self.num_variants,
                "Provided variants slice length does not match the one used to build the KdTree"
            );
            for (i, variant) in variants.iter().enumerate() {
                let p = variant.point();
                println!(
                    "{},{},{},{:?},{},{},{}",
                    variant.svtype, variant.vcf_id, variant.id, variant.trid, i, p[0], p[1]
                );
            }
        }

        pub fn size(&self) -> usize {
            self.num_variants
        }
    }
}

pub use kdtree_impl::{KdTreeItemType, VariantKdTree};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::svtype::SvType;
    use crate::core::variant_internal::VariantInternal;
    use approx::assert_relative_eq;
    use kiddo::Manhattan;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn test_weighted_manhattan_2d() {
        let point1 = [0.0, 0.0];
        let point2 = [3.0, 4.0];

        let d1 = Manhattan::dist(&point1, &point2);
        let d2 = WeightedManhattan::dist(&point1, &point2);

        assert_relative_eq!(d1, 7.0, epsilon = 1e-10);
        assert_relative_eq!(d2, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_variant_kd_tree() {
        let variants = vec![
            VariantInternal::from_parts(2, "0".to_string(), SvType::INSERTION, 30.0, 36.0).unwrap(),
            VariantInternal::from_parts(1, "1".to_string(), SvType::INSERTION, 20.0, 24.0).unwrap(),
            VariantInternal::from_parts(0, "2".to_string(), SvType::INSERTION, 3.0, 20.0).unwrap(),
            VariantInternal::from_parts(3, "3".to_string(), SvType::INSERTION, 4.0, 5.0).unwrap(),
        ];

        let kdtree = VariantKdTree::new(&variants);
        let query =
            VariantInternal::from_parts(3, "0".to_string(), SvType::INSERTION, 4.0, 5.0).unwrap();
        let nearest = kdtree.knn(&query, 3);
        assert_eq!(nearest.len(), 3);
        assert_eq!(nearest[0].item, 3 as KdTreeItemType);
        assert_eq!(nearest[1].item, 2 as KdTreeItemType);
        assert_eq!(nearest[2].item, 1 as KdTreeItemType);
    }

    #[test]
    fn test_knn_ordering() {
        let mut rng = StdRng::seed_from_u64(42);
        let num_points = 1000;
        let k = 50;
        let variants: Vec<VariantInternal> = (0..num_points)
            .map(|i| {
                let start: f32 = rng.random_range(0.0..1000.0);
                let svlen: f32 = rng.random_range(10.0..500.0);
                VariantInternal::from_parts(
                    0,             // sample_id
                    i.to_string(), // id
                    SvType::INSERTION,
                    start,
                    svlen,
                )
                .unwrap()
            })
            .collect();

        let kdtree = VariantKdTree::new(&variants);

        let query_start: f32 = rng.random_range(0.0..1000.0);
        let query_svlen: f32 = rng.random_range(10.0..500.0);
        let query_variant = VariantInternal::from_parts(
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

        for i in 0..(nearest.len() - 1) {
            assert!(
                nearest[i].distance <= nearest[i + 1].distance,
                "KNN results are not sorted by distance: item {} distance {} > item {} distance {}",
                i,
                nearest[i].distance,
                i + 1,
                nearest[i + 1].distance
            );
        }
    }

    #[test]
    fn test_knn() {
        dbg!(SquaredEuclidean::dist(
            &[3299543.0, 732.0],
            &[4065783.0, 570.0]
        ));
    }
}
