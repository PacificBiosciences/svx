use std::{cmp::Ordering, fmt};

#[derive(Debug)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub dist: f64,
    dist_key: i64,
    from_info_hash: i32,
    to_info_hash: i32,
    from_id_rank: u32,
    to_id_rank: u32,
}

impl Edge {
    fn dist_key(dist: f64) -> i64 {
        if !dist.is_finite() {
            return i64::MAX;
        }
        let scaled = (dist * 1e9_f64).round();
        if scaled >= i64::MAX as f64 {
            i64::MAX
        } else if scaled <= i64::MIN as f64 {
            i64::MIN
        } else {
            scaled as i64
        }
    }

    pub fn new(
        from: usize,
        to: usize,
        dist: f64,
        from_info_hash: i32,
        to_info_hash: i32,
        from_id_rank: u32,
        to_id_rank: u32,
    ) -> Self {
        Self {
            from,
            to,
            dist,
            dist_key: Self::dist_key(dist),
            from_info_hash,
            to_info_hash,
            from_id_rank,
            to_id_rank,
        }
    }
}

impl fmt::Display for Edge {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Edge[from={}, to={}, dist={}]",
            self.from, self.to, self.dist
        )
    }
}

impl Ord for Edge {
    fn cmp(&self, other: &Self) -> Ordering {
        // Primary: smaller distance first (with 1e-9 tie tolerance)
        // Ties: from.INFO_hash, to.INFO_hash, from.ID, to.ID
        let dist_cmp = if self.dist_key.abs_diff(other.dist_key) > 1 {
            other.dist_key.cmp(&self.dist_key)
        } else {
            Ordering::Equal
        };

        dist_cmp
            .then_with(|| other.from_info_hash.cmp(&self.from_info_hash))
            .then_with(|| other.to_info_hash.cmp(&self.to_info_hash))
            .then_with(|| other.from_id_rank.cmp(&self.from_id_rank))
            .then_with(|| other.to_id_rank.cmp(&self.to_id_rank))
            .then_with(|| other.from.cmp(&self.from))
            .then_with(|| other.to.cmp(&self.to))
    }
}

impl PartialOrd for Edge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Edge {}
