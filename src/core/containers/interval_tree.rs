use num_traits::{One, Zero};
use std::{
    cmp::Ordering,
    fmt,
    ops::{Add, Div},
};

#[derive(Debug, Clone, PartialEq)]
pub struct Interval<Scalar, Value> {
    pub start: Scalar,
    pub stop: Scalar,
    pub value: Value,
}

impl<Scalar, Value> Interval<Scalar, Value>
where
    Scalar: PartialOrd + Copy,
{
    pub fn new(s: Scalar, e: Scalar, v: Value) -> Self {
        let (start, stop) = if s <= e { (s, e) } else { (e, s) };
        Self {
            start,
            stop,
            value: v,
        }
    }
}

impl<Scalar, Value> fmt::Display for Interval<Scalar, Value>
where
    Scalar: fmt::Display,
    Value: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Interval({}, {}): {}", self.start, self.stop, self.value)
    }
}

#[derive(Default)]
pub struct IntervalTree<Scalar, Value> {
    intervals: Vec<Interval<Scalar, Value>>,
    left: Option<Box<IntervalTree<Scalar, Value>>>,
    right: Option<Box<IntervalTree<Scalar, Value>>>,
    center: Scalar,
    // Cached bounds of `self.intervals` (the "center" bucket) for fast pruning in queries.
    // IMPORTANT: `self.intervals` are sorted by `start`, so `intervals.last().stop` is NOT
    // necessarily the maximum stop; we must track it explicitly.
    min_start: Scalar,
    max_stop: Scalar,
}

impl<Scalar, Value> IntervalTree<Scalar, Value>
where
    Scalar:
        PartialOrd + Copy + Add<Output = Scalar> + Div<Output = Scalar> + Zero + One + fmt::Debug,
{
    #[inline]
    fn default() -> Self {
        Self {
            intervals: Vec::new(),
            left: None,
            right: None,
            center: Scalar::zero(),
            min_start: Scalar::zero(),
            max_stop: Scalar::zero(),
        }
    }

    pub fn new(mut intervals: Vec<Interval<Scalar, Value>>) -> Self {
        let n = intervals.len();
        if n == 0 {
            return Self::default();
        }
        intervals.sort_unstable_by(|a, b| a.start.partial_cmp(&b.start).unwrap_or(Ordering::Equal));
        Self::build_tree(intervals, 16, 64, 512, None, None)
    }

    fn build_tree(
        intervals: Vec<Interval<Scalar, Value>>,
        depth: usize,
        minbucket: usize,
        maxbucket: usize,
        left_extent: Option<Scalar>,
        right_extent: Option<Scalar>,
    ) -> Self {
        let left_extent = left_extent.unwrap_or_else(|| intervals.first().unwrap().start);
        let right_extent = right_extent.unwrap_or_else(|| {
            intervals
                .iter()
                .map(|i| i.stop)
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap()
        });

        let center = (left_extent + right_extent) / (Scalar::one() + Scalar::one());

        if depth == 0 || (intervals.len() < minbucket && intervals.len() < maxbucket) {
            let (min_start, max_stop) = if intervals.is_empty() {
                (center, center)
            } else {
                let min_start = intervals.first().unwrap().start;
                let max_stop = intervals
                    .iter()
                    .map(|i| i.stop)
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap();
                (min_start, max_stop)
            };
            return Self {
                intervals,
                left: None,
                right: None,
                center,
                min_start,
                max_stop,
            };
        }

        let (lefts, centers): (Vec<_>, Vec<_>) =
            intervals.into_iter().partition(|i| i.stop < center);
        let (centers, rights): (Vec<_>, Vec<_>) =
            centers.into_iter().partition(|i| i.start <= center);

        let (min_start, max_stop) = if centers.is_empty() {
            (center, center)
        } else {
            let min_start = centers.first().unwrap().start;
            let max_stop = centers
                .iter()
                .map(|i| i.stop)
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap();
            (min_start, max_stop)
        };

        let left = (!lefts.is_empty()).then(|| {
            Box::new(Self::build_tree(
                lefts,
                depth - 1,
                minbucket,
                maxbucket,
                Some(left_extent),
                Some(center),
            ))
        });
        let right = (!rights.is_empty()).then(|| {
            Box::new(Self::build_tree(
                rights,
                depth - 1,
                minbucket,
                maxbucket,
                Some(center),
                Some(right_extent),
            ))
        });

        Self {
            intervals: centers,
            left,
            right,
            center,
            min_start,
            max_stop,
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.left.as_ref().is_none_or(|l| l.is_empty())
            && self.intervals.is_empty()
            && self.right.as_ref().is_none_or(|r| r.is_empty())
    }

    fn visit_near_while<F>(&self, start: Scalar, stop: Scalar, f: &mut F) -> bool
    where
        F: FnMut(&Interval<Scalar, Value>) -> bool,
    {
        if !(self.intervals.is_empty() || stop < self.min_start || start > self.max_stop) {
            for interval in &self.intervals {
                if interval.stop >= start && interval.start <= stop && !f(interval) {
                    return false;
                }
            }
        }

        if start <= self.center
            && self
                .left
                .as_ref()
                .is_some_and(|left| !left.visit_near_while(start, stop, f))
        {
            return false;
        }

        if stop >= self.center
            && self
                .right
                .as_ref()
                .is_some_and(|right| !right.visit_near_while(start, stop, f))
        {
            return false;
        }

        true
    }

    pub fn visit_near<F>(&self, start: Scalar, stop: Scalar, f: &mut F)
    where
        F: FnMut(&Interval<Scalar, Value>),
    {
        self.visit_near_while(start, stop, &mut |i| {
            f(i);
            true
        });
    }

    #[inline]
    pub fn visit_overlapping<F>(&self, start: Scalar, stop: Scalar, f: &mut F)
    where
        F: FnMut(&Interval<Scalar, Value>),
    {
        self.visit_near(start, stop, &mut |i| {
            if i.stop >= start && i.start <= stop {
                f(i);
            }
        });
    }

    #[inline]
    pub fn visit_contained<F>(&self, start: Scalar, stop: Scalar, f: &mut F)
    where
        F: FnMut(&Interval<Scalar, Value>),
    {
        self.visit_near(start, stop, &mut |i| {
            if i.start >= start && i.stop <= stop {
                f(i);
            }
        });
    }

    #[inline]
    pub fn visit_containing_while<F>(&self, start: Scalar, stop: Scalar, f: &mut F) -> bool
    where
        F: FnMut(&Interval<Scalar, Value>) -> bool,
    {
        self.visit_near_while(start, stop, &mut |i| {
            if i.start <= start && i.stop >= stop {
                return f(i);
            }
            true
        })
    }

    #[inline]
    pub fn any_containing(&self, start: Scalar, stop: Scalar) -> bool {
        let mut found = false;
        self.visit_containing_while(start, stop, &mut |_| {
            found = true;
            false
        });
        found
    }

    pub fn find_overlapping(&self, start: Scalar, stop: Scalar) -> Vec<Interval<Scalar, Value>>
    where
        Value: Clone,
    {
        let mut result = Vec::new();
        self.visit_overlapping(start, stop, &mut |i| result.push(i.clone()));
        result
    }

    pub fn find_contained(&self, start: Scalar, stop: Scalar) -> Vec<Interval<Scalar, Value>>
    where
        Value: Clone,
    {
        let mut result = Vec::new();
        self.visit_contained(start, stop, &mut |i| result.push(i.clone()));
        result
    }

    pub fn find_containing(&self, start: Scalar, stop: Scalar) -> Vec<Interval<Scalar, Value>>
    where
        Value: Clone,
    {
        let mut result = Vec::new();
        self.visit_containing_while(start, stop, &mut |i| {
            result.push(i.clone());
            true
        });
        result
    }

    pub fn visit_all<F>(&self, f: &mut F)
    where
        F: FnMut(&Interval<Scalar, Value>),
    {
        if let Some(ref left) = self.left {
            left.visit_all(f);
        }
        self.intervals.iter().for_each(&mut *f);
        if let Some(ref right) = self.right {
            right.visit_all(f);
        }
    }
}

impl<Scalar: fmt::Display, Value: fmt::Display> fmt::Display for IntervalTree<Scalar, Value>
where
    Scalar: PartialOrd + Copy + Add<Output = Scalar> + Div<Output = Scalar> + Zero + One,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_out(f, 0)
    }
}

impl<Scalar, Value> IntervalTree<Scalar, Value>
where
    Scalar:
        fmt::Display + PartialOrd + Copy + Add<Output = Scalar> + Div<Output = Scalar> + Zero + One,
    Value: fmt::Display,
{
    fn write_out(&self, f: &mut fmt::Formatter<'_>, depth: usize) -> fmt::Result {
        let pad = |f: &mut fmt::Formatter<'_>, depth: usize| {
            for _ in 0..depth {
                write!(f, " ")?;
            }
            Ok(())
        };
        pad(f, depth)?;
        writeln!(f, "center: {}", self.center)?;
        for interval in &self.intervals {
            pad(f, depth)?;
            writeln!(f, "{}", interval)?;
        }
        if let Some(ref left) = self.left {
            pad(f, depth)?;
            writeln!(f, "left:")?;
            left.write_out(f, depth + 1)?;
        } else {
            pad(f, depth)?;
            writeln!(f, "left: None")?;
        }
        if let Some(ref right) = self.right {
            pad(f, depth)?;
            writeln!(f, "right:")?;
            right.write_out(f, depth + 1)?;
        } else {
            pad(f, depth)?;
            writeln!(f, "right: None")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const MIN: i64 = -1_000_000;
    const MAX: i64 = 1_000_000;

    #[test]
    fn test_interval_display() {
        let interval = Interval::new(1, 5, "value");
        let output = format!("{}", interval);
        assert_eq!(output, "Interval(1, 5): value");
    }

    #[test]
    fn overlap_not_pruned_by_last_stop() {
        let tree = IntervalTree::new(vec![
            Interval::new(0i32, 100i32, 1usize),
            Interval::new(40i32, 60i32, 2usize),
        ]);

        let hits = tree.find_overlapping(90, 95);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0], Interval::new(0, 100, 1));
    }

    #[test]
    fn interval_new_normalizes_endpoints() {
        let interval = Interval::new(5, 2, 42);
        assert_eq!(interval.start, 2);
        assert_eq!(interval.stop, 5);
        assert_eq!(interval.value, 42);

        let interval = Interval::new(3, 3, 42);
        assert_eq!(interval.start, 3);
        assert_eq!(interval.stop, 3);

        let interval = Interval::new(-5, -2, 42);
        assert_eq!(interval.start, -5);
        assert_eq!(interval.stop, -2);
    }

    #[test]
    fn empty_tree_is_empty() {
        let t: IntervalTree<i32, i32> = IntervalTree::new(Vec::new());
        assert!(t.is_empty());

        assert!(t.find_overlapping(0, 0).is_empty());
        assert!(t.find_contained(0, 0).is_empty());
        assert!(t.find_containing(0, 0).is_empty());
    }

    #[test]
    fn any_containing_matches_find_containing_non_empty() {
        let t = IntervalTree::new(vec![
            Interval::new(0, 10, 1usize),
            Interval::new(20, 30, 2usize),
            Interval::new(25, 40, 3usize),
        ]);

        let queries = vec![(0, 1), (5, 6), (10, 11), (20, 21), (31, 32), (24, 29)];
        for (qs, qe) in queries {
            let expected = !t.find_containing(qs, qe).is_empty();
            assert_eq!(t.any_containing(qs, qe), expected);
        }
    }

    #[test]
    fn visit_containing_while_stops_after_first_match() {
        let t = IntervalTree::new(vec![
            Interval::new(0, 100, 1usize),
            Interval::new(10, 90, 2usize),
            Interval::new(20, 80, 3usize),
        ]);

        let mut visited = 0usize;
        let completed = t.visit_containing_while(50, 51, &mut |_| {
            visited += 1;
            false
        });

        assert!(!completed);
        assert_eq!(visited, 1);
    }

    #[test]
    fn test_custom_scalar_type() {
        use std::ops::Mul;

        #[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
        struct CustomScalar(i32);

        impl Mul for CustomScalar {
            type Output = Self;
            fn mul(self, other: Self) -> Self {
                CustomScalar(self.0 * other.0)
            }
        }

        impl Add for CustomScalar {
            type Output = Self;
            fn add(self, other: Self) -> Self {
                CustomScalar(self.0 + other.0)
            }
        }

        impl Div for CustomScalar {
            type Output = Self;
            fn div(self, other: Self) -> Self {
                CustomScalar(self.0 / other.0)
            }
        }

        impl Zero for CustomScalar {
            fn zero() -> Self {
                CustomScalar(0)
            }

            fn is_zero(&self) -> bool {
                self.0 == 0
            }
        }

        impl One for CustomScalar {
            fn one() -> Self {
                CustomScalar(1)
            }
        }

        impl fmt::Display for CustomScalar {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}", self.0)
            }
        }

        let intervals = vec![
            Interval::new(CustomScalar(1), CustomScalar(5), "a"),
            Interval::new(CustomScalar(3), CustomScalar(7), "b"),
        ];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(CustomScalar(4), CustomScalar(4));
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_display_interval_tree() {
        let intervals = vec![Interval::new(1, 5, "a"), Interval::new(3, 7, "b")];
        let t = IntervalTree::new(intervals);
        let output = format!("{}", t);
        assert!(output.contains(
            "center: 4\nInterval(1, 5): a\nInterval(3, 7): b\nleft: None\nright: None\n"
        ));
    }

    #[test]
    fn test_interval_with_nan() {
        let intervals = vec![Interval::new(f64::NAN, f64::NAN, 42)];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(f64::NAN, f64::NAN);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_interval_with_infinity() {
        let intervals = vec![Interval::new(f64::NEG_INFINITY, f64::INFINITY, 42)];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(0.0, 1.0);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].value, 42);
    }

    #[test]
    fn test_large_number_of_intervals() {
        let mut intervals = Vec::new();
        for i in 0..1000 {
            intervals.push(Interval::new(i, i + 10, i));
        }
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(500, 505);
        assert!(!result.is_empty());
        for interval in result {
            assert!(interval.start <= 505 && interval.stop >= 500);
        }
    }

    #[test]
    fn wild_search_values() {
        let intervals = vec![Interval::new(1.0f64, 2.0f64, 0)];
        let t = IntervalTree::new(intervals);
        let inf = f64::INFINITY;
        let nan = f64::NAN;

        let sanity_results = t.find_overlapping(inf, inf);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(-inf, inf);
        assert_eq!(sanity_results.len(), 1);

        let sanity_results = t.find_overlapping(0.0, inf);
        assert_eq!(sanity_results.len(), 1);

        let sanity_results = t.find_overlapping(0.5, inf);
        assert_eq!(sanity_results.len(), 1);

        let sanity_results = t.find_overlapping(1.1, inf);
        assert_eq!(sanity_results.len(), 1);

        let sanity_results = t.find_overlapping(-inf, 1.0);
        assert_eq!(sanity_results.len(), 1);

        let sanity_results = t.find_overlapping(-inf, 0.5);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(-inf, 0.0);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(-inf, -0.1);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(nan, nan);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(-nan, nan);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(nan, 1.0);
        assert_eq!(sanity_results.len(), 0);

        let sanity_results = t.find_overlapping(0.0, nan);
        assert_eq!(sanity_results.len(), 0);
    }

    #[inline]
    fn norm(a: i64, b: i64) -> (i64, i64) {
        if a <= b { (a, b) } else { (b, a) }
    }

    #[inline]
    fn intervals_from_raw(raw: Vec<(i64, i64)>) -> Vec<Interval<i64, usize>> {
        raw.into_iter()
            .enumerate()
            .map(|(id, (s, e))| Interval::new(s, e, id))
            .collect()
    }

    #[inline]
    fn sorted(mut ids: Vec<usize>) -> Vec<usize> {
        ids.sort_unstable();
        ids
    }

    #[inline]
    fn naive_overlapping(intervals: &[Interval<i64, usize>], qs: i64, qe: i64) -> Vec<usize> {
        sorted(
            intervals
                .iter()
                .filter(|i| i.stop >= qs && i.start <= qe)
                .map(|i| i.value)
                .collect(),
        )
    }

    #[inline]
    fn naive_contained(intervals: &[Interval<i64, usize>], qs: i64, qe: i64) -> Vec<usize> {
        sorted(
            intervals
                .iter()
                .filter(|i| i.start >= qs && i.stop <= qe)
                .map(|i| i.value)
                .collect(),
        )
    }

    #[inline]
    fn naive_containing(intervals: &[Interval<i64, usize>], qs: i64, qe: i64) -> Vec<usize> {
        sorted(
            intervals
                .iter()
                .filter(|i| i.start <= qs && i.stop >= qe)
                .map(|i| i.value)
                .collect(),
        )
    }

    proptest! {
        #[test]
        fn prop_overlapping_matches_naive(
            raw in proptest::collection::vec((MIN..=MAX, MIN..=MAX), 0..200),
            q in (MIN..=MAX, MIN..=MAX),
        ) {
            let intervals = intervals_from_raw(raw);
            let tree = IntervalTree::new(intervals.clone());
            let (qs, qe) = norm(q.0, q.1);

            let naive = naive_overlapping(&intervals, qs, qe);
            let got = sorted(
                tree.find_overlapping(qs, qe)
                    .into_iter()
                    .map(|i| i.value)
                    .collect(),
            );
            prop_assert_eq!(got, naive);
        }

        #[test]
        fn prop_contained_matches_naive(
            raw in proptest::collection::vec((MIN..=MAX, MIN..=MAX), 0..200),
            q in (MIN..=MAX, MIN..=MAX),
        ) {
            let intervals = intervals_from_raw(raw);
            let tree = IntervalTree::new(intervals.clone());
            let (qs, qe) = norm(q.0, q.1);

            let naive = naive_contained(&intervals, qs, qe);
            let got = sorted(
                tree.find_contained(qs, qe)
                    .into_iter()
                    .map(|i| i.value)
                    .collect(),
            );
            prop_assert_eq!(got, naive);
        }

        #[test]
        fn prop_containing_matches_naive(
            raw in proptest::collection::vec((MIN..=MAX, MIN..=MAX), 0..200),
            q in (MIN..=MAX, MIN..=MAX),
        ) {
            let intervals = intervals_from_raw(raw);
            let tree = IntervalTree::new(intervals.clone());
            let (qs, qe) = norm(q.0, q.1);

            let naive = naive_containing(&intervals, qs, qe);
            let got = sorted(
                tree.find_containing(qs, qe)
                    .into_iter()
                    .map(|i| i.value)
                    .collect(),
            );
            prop_assert_eq!(got, naive);
        }

        #[test]
        fn prop_visit_near_matches_naive_overlap(
            raw in proptest::collection::vec((MIN..=MAX, MIN..=MAX), 0..200),
            q in (MIN..=MAX, MIN..=MAX),
        ) {
            let intervals = intervals_from_raw(raw);
            let tree = IntervalTree::new(intervals.clone());
            let (qs, qe) = norm(q.0, q.1);

            let naive = naive_overlapping(&intervals, qs, qe);
            let mut got = Vec::new();
            tree.visit_near(qs, qe, &mut |i| got.push(i.value));
            prop_assert_eq!(sorted(got), naive);
        }

        #[test]
        fn prop_visit_all_returns_every_interval_exactly_once(
            raw in proptest::collection::vec((MIN..=MAX, MIN..=MAX), 0..200),
        ) {
            let intervals = intervals_from_raw(raw);
            let n = intervals.len();
            let tree = IntervalTree::new(intervals);

            let mut got = Vec::with_capacity(n);
            tree.visit_all(&mut |i| got.push(i.value));
            prop_assert_eq!(sorted(got), (0..n).collect::<Vec<_>>());
        }

        #[test]
        fn prop_is_empty_matches_input_len(
            raw in proptest::collection::vec((MIN..=MAX, MIN..=MAX), 0..200),
        ) {
            let intervals = intervals_from_raw(raw);
            let tree = IntervalTree::new(intervals.clone());
            prop_assert_eq!(tree.is_empty(), intervals.is_empty());
        }
    }
}
