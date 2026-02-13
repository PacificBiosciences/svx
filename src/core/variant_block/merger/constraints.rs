use super::VariantMerger;
use crate::{DISTANCE_OFFSET, cli::MergeConstraint, core::variant::VariantInternal};

pub(super) struct CentroidState {
    sum_start: Vec<f64>,
    sum_end: Vec<f64>,
    count: Vec<usize>,
}

impl CentroidState {
    pub(super) fn new(variants: &[VariantInternal]) -> Self {
        let mut sum_start = Vec::with_capacity(variants.len());
        let mut sum_end = Vec::with_capacity(variants.len());
        let mut count = Vec::with_capacity(variants.len());
        for v in variants {
            sum_start.push(v.start);
            sum_end.push(v.end);
            count.push(1);
        }
        Self {
            sum_start,
            sum_end,
            count,
        }
    }
}

impl VariantMerger<'_> {
    pub(super) fn effective_max_dist_for_pair(
        &self,
        v1: &VariantInternal,
        v2: &VariantInternal,
    ) -> f64 {
        let v1_is_tr = v1.trid.is_some();
        let v2_is_tr = v2.trid.is_some();
        if v1_is_tr && v2_is_tr && v1.trid.as_ref().unwrap().id == v2.trid.as_ref().unwrap().id {
            f64::from(self.args.tr_max_dist)
        } else if self.args.require_mutual_distance {
            v1.max_dist.min(v2.max_dist)
        } else {
            v1.max_dist.max(v2.max_dist)
        }
    }

    fn pair_is_distance_mergeable(&self, v1: &VariantInternal, v2: &VariantInternal) -> bool {
        if let (Some(tr1), Some(tr2)) = (&v1.trid, &v2.trid) {
            if tr1.id != tr2.id {
                return false;
            }
        }

        v1.distance(v2) <= self.effective_max_dist_for_pair(v1, v2) + DISTANCE_OFFSET
    }

    fn union_root_order(&self, mut root_a: usize, mut root_b: usize) -> (usize, usize) {
        if self.forest.parent[root_a] < self.forest.parent[root_b] {
            std::mem::swap(&mut root_a, &mut root_b);
        }

        // After the swap, root_b is the "winner" (larger component or tie-break), and root_a is absorbed.
        (root_a, root_b)
    }

    fn merge_members(&mut self, absorbed_root: usize, new_root: usize) {
        let members = self
            .members
            .as_mut()
            .expect("merge-constraint requires membership tracking");
        let mut absorbed = std::mem::take(&mut members[absorbed_root]);
        members[new_root].append(&mut absorbed);
    }

    fn try_union_clique(&mut self, from: usize, to: usize) -> bool {
        if !self.forest.can_union(from, to) {
            return false;
        }

        let root_a = self.forest.find(from);
        let root_b = self.forest.find(to);
        if root_a == root_b {
            return false;
        }

        {
            let members = self
                .members
                .as_ref()
                .expect("clique merge-constraint requires membership tracking");
            let (small_root, large_root) = if members[root_a].len() <= members[root_b].len() {
                (root_a, root_b)
            } else {
                (root_b, root_a)
            };

            for &i in &members[small_root] {
                let vi = &self.variants[i];
                for &j in &members[large_root] {
                    if !self.pair_is_distance_mergeable(vi, &self.variants[j]) {
                        return false;
                    }
                }
            }
        }

        let (absorbed_root, new_root) = self.union_root_order(root_a, root_b);
        let unioned = self.forest.try_union(from, to);
        if unioned {
            self.merge_members(absorbed_root, new_root);
        }
        unioned
    }

    fn try_union_centroid(&mut self, from: usize, to: usize) -> bool {
        if !self.forest.can_union(from, to) {
            return false;
        }

        let root_a = self.forest.find(from);
        let root_b = self.forest.find(to);
        if root_a == root_b {
            return false;
        }

        let (centroid_start, centroid_end) = {
            let centroid = self
                .centroid_state
                .as_ref()
                .expect("centroid merge-constraint requires centroid state");
            let sum_start = centroid.sum_start[root_a] + centroid.sum_start[root_b];
            let sum_end = centroid.sum_end[root_a] + centroid.sum_end[root_b];
            let count = centroid.count[root_a] + centroid.count[root_b];
            debug_assert!(count > 0);
            (sum_start / count as f64, sum_end / count as f64)
        };

        {
            let members = self
                .members
                .as_ref()
                .expect("centroid merge-constraint requires membership tracking");
            for &idx in members[root_a].iter().chain(members[root_b].iter()) {
                let v = &self.variants[idx];
                let d_start = v.start - centroid_start;
                let d_end = v.end - centroid_end;
                let dist = (d_start * d_start + d_end * d_end).sqrt();
                if dist > v.max_dist + DISTANCE_OFFSET {
                    return false;
                }
            }
        }

        let (absorbed_root, new_root) = self.union_root_order(root_a, root_b);
        let unioned = self.forest.try_union(from, to);
        if unioned {
            self.merge_members(absorbed_root, new_root);

            let centroid = self
                .centroid_state
                .as_mut()
                .expect("centroid merge-constraint requires centroid state");
            centroid.sum_start[new_root] += centroid.sum_start[absorbed_root];
            centroid.sum_end[new_root] += centroid.sum_end[absorbed_root];
            centroid.count[new_root] += centroid.count[absorbed_root];
            centroid.sum_start[absorbed_root] = 0.0;
            centroid.sum_end[absorbed_root] = 0.0;
            centroid.count[absorbed_root] = 0;
        }
        unioned
    }

    pub(super) fn try_union_with_constraint(&mut self, from: usize, to: usize) -> bool {
        match self.merge_constraint {
            MergeConstraint::None => self.forest.try_union(from, to),
            MergeConstraint::BboxDiameter => self.forest.try_union_with_diameter(
                from,
                to,
                Some(self.effective_max_dist_for_pair(&self.variants[from], &self.variants[to])),
            ),
            MergeConstraint::Clique => self.try_union_clique(from, to),
            MergeConstraint::Centroid => self.try_union_centroid(from, to),
        }
    }
}
