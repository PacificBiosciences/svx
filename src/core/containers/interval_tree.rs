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
            return Self {
                intervals,
                left: None,
                right: None,
                center,
            };
        }

        let (lefts, centers): (Vec<_>, Vec<_>) =
            intervals.into_iter().partition(|i| i.stop < center);
        let (centers, rights): (Vec<_>, Vec<_>) =
            centers.into_iter().partition(|i| i.start <= center);

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
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.left.as_ref().is_none_or(|l| l.is_empty())
            && self.intervals.is_empty()
            && self.right.as_ref().is_none_or(|r| r.is_empty())
    }

    pub fn visit_near<F>(&self, start: Scalar, stop: Scalar, f: &mut F)
    where
        F: FnMut(&Interval<Scalar, Value>),
    {
        if !(self.intervals.is_empty()
            || stop < self.intervals[0].start
            || start > self.intervals.last().unwrap().stop)
        {
            for interval in &self.intervals {
                if interval.stop >= start && interval.start <= stop {
                    f(interval);
                }
            }
        }

        if start <= self.center {
            if let Some(ref left) = self.left {
                left.visit_near(start, stop, f);
            }
        }

        if stop >= self.center {
            if let Some(ref right) = self.right {
                right.visit_near(start, stop, f);
            }
        }
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
        self.visit_overlapping(start, stop, &mut |i| {
            if i.start <= start && i.stop >= stop {
                result.push(i.clone());
            }
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
    use rand::Rng;
    use std::cmp::min;
    use std::collections::HashSet;
    use std::f64;
    use std::time::Instant;

    fn random_interval<Value>(
        min_start: i32,
        max_stop: i32,
        min_len: i32,
        max_len: i32,
        value: Value,
    ) -> Interval<i32, Value>
    where
        Value: Clone,
    {
        let mut rng = rand::rng();
        let len = rng.random_range(min_len..=max_len);
        let start = rng.random_range(min_start..=max_stop - len);
        let stop = min(start + len - 1, max_stop);
        Interval::new(start, stop, value)
    }

    #[test]
    fn test_interval_display() {
        let interval = Interval::new(1, 5, "value");
        let output = format!("{}", interval);
        assert_eq!(output, "Interval(1, 5): value");
    }

    #[test]
    fn test_interval_new_ordering() {
        let interval = Interval::new(5, 2, 42);
        assert_eq!(interval.start, 2);
        assert_eq!(interval.stop, 5);
        assert_eq!(interval.value, 42);
    }

    #[test]
    fn test_empty_tree_default_center() {
        let t: IntervalTree<i32, i32> = IntervalTree::new(Vec::new());
        assert!(t.is_empty());
        assert_eq!(t.center, 0);
    }

    #[test]
    fn test_default_impl() {
        let t: IntervalTree<i32, i32> = IntervalTree::default();
        assert!(t.is_empty());
        assert_eq!(t.center, 0);
        assert!(t.left.is_none());
        assert!(t.right.is_none());
        assert!(t.intervals.is_empty());
    }

    #[test]
    fn test_build_tree_with_empty_intervals() {
        let intervals: Vec<Interval<i32, i32>> = Vec::new();
        let t = IntervalTree::new(intervals);
        assert!(t.is_empty());
        assert_eq!(t.center, 0);
    }

    #[test]
    fn test_center_initialization_with_non_zero_scalar() {
        let t: IntervalTree<f64, i32> = IntervalTree::default();
        assert!(t.is_empty());
        assert_eq!(t.center, 0.0);
    }

    #[test]
    fn test_zero_length_intervals_in_search() {
        let intervals = vec![Interval::new(5, 5, "zero")];
        let t = IntervalTree::new(intervals);

        let result_overlapping = t.find_overlapping(5, 5);
        assert_eq!(result_overlapping.len(), 1);
        assert_eq!(result_overlapping[0].value, "zero");

        let result_contained = t.find_contained(5, 5);
        assert_eq!(result_contained.len(), 1);
        assert_eq!(result_contained[0].value, "zero");
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
    fn test_interval_new_zero_length() {
        let interval = Interval::new(3, 3, 42);
        assert_eq!(interval.start, 3);
        assert_eq!(interval.stop, 3);
        assert_eq!(interval.value, 42);
    }

    #[test]
    fn test_interval_new_negative_scalars() {
        let interval = Interval::new(-5, -2, 42);
        assert_eq!(interval.start, -5);
        assert_eq!(interval.stop, -2);
        assert_eq!(interval.value, 42);
    }

    #[test]
    fn test_is_empty_tree() {
        let t: IntervalTree<isize, i32> = IntervalTree::new(Vec::new());
        assert!(t.is_empty());

        let intervals = vec![Interval::new(1, 5, 42)];
        let t = IntervalTree::new(intervals);
        assert!(!t.is_empty());
    }

    #[test]
    fn test_zero_length_intervals() {
        let intervals = vec![Interval::new(5, 5, 42)];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(5, 5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].start, 5);
        assert_eq!(result[0].stop, 5);
        assert_eq!(result[0].value, 42);
    }

    #[test]
    fn test_zero_length_interval_no_overlap() {
        let intervals = vec![Interval::new(5, 5, 42)];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(4, 4);
        assert_eq!(result.len(), 0);
        let result = t.find_overlapping(6, 6);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_negative_intervals() {
        let intervals = vec![Interval::new(-10, -5, 42)];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(-7, -6);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].start, -10);
        assert_eq!(result[0].stop, -5);
        assert_eq!(result[0].value, 42);
    }

    #[test]
    fn test_overlapping_intervals() {
        let intervals = vec![
            Interval::new(1, 5, "a"),
            Interval::new(3, 7, "b"),
            Interval::new(4, 6, "c"),
        ];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(4, 4);
        assert_eq!(result.len(), 3);
        let values: Vec<_> = result.iter().map(|i| i.value).collect();
        assert!(values.contains(&"a"));
        assert!(values.contains(&"b"));
        assert!(values.contains(&"c"));
    }

    #[test]
    fn test_nested_intervals() {
        let intervals = vec![
            Interval::new(1, 10, "outer"),
            Interval::new(3, 7, "middle"),
            Interval::new(4, 5, "inner"),
        ];
        let t = IntervalTree::new(intervals);
        let result = t.find_contained(2, 8);
        assert_eq!(result.len(), 2);
        let values: Vec<_> = result.iter().map(|i| i.value).collect();
        assert!(values.contains(&"middle"));
        assert!(values.contains(&"inner"));
    }

    #[test]
    fn test_boundary_overlap() {
        let intervals = vec![Interval::new(1, 5, "a"), Interval::new(5, 10, "b")];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(5, 5);
        assert_eq!(result.len(), 2);
        let values: Vec<_> = result.iter().map(|i| i.value).collect();
        assert!(values.contains(&"a"));
        assert!(values.contains(&"b"));
    }

    #[test]
    fn test_adjacent_intervals_no_overlap() {
        let intervals = vec![Interval::new(1, 5, "a"), Interval::new(5, 10, "b")];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(0, 1);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].value, "a");

        let result = t.find_overlapping(11, 15);
        assert_eq!(result.len(), 0);
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
    fn test_center_computation() {
        let intervals = vec![Interval::new(1, 7, "a"), Interval::new(3, 5, "b")];
        let t = IntervalTree::new(intervals);

        let expected_center = 1 + (7 - 1) / 2;
        assert_eq!(t.center, expected_center);
    }

    #[test]
    fn test_center_computation_f64() {
        let intervals = vec![Interval::new(1.0, 7.0, "a"), Interval::new(3.0, 5.0, "b")];
        let t = IntervalTree::new(intervals);

        let expected_center = 1.0 + (7.0 - 1.0) / 2.0;
        assert_eq!(t.center, expected_center);
    }

    #[test]
    fn test_center_computation_multiple() {
        let intervals = vec![
            Interval::new(4, 7, "a"),
            Interval::new(90, 130, "b"),
            Interval::new(3, 5, "c"),
            Interval::new(12, 15, "d"),
        ];
        let t = IntervalTree::new(intervals);
        let expected_center = 3 + (130 - 3) / 2;
        assert_eq!(t.center, expected_center);
    }

    #[test]
    fn test_visit_all() {
        let intervals = vec![
            Interval::new(1, 5, "a"),
            Interval::new(3, 7, "b"),
            Interval::new(6, 10, "c"),
        ];
        let t = IntervalTree::new(intervals.clone());

        let mut visited_intervals = Vec::new();
        t.visit_all(&mut |interval| {
            visited_intervals.push(interval.clone());
        });

        assert_eq!(visited_intervals.len(), intervals.len());
        for interval in intervals {
            assert!(visited_intervals.contains(&interval));
        }
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
    fn test_find_overlapping_no_overlap() {
        let intervals = vec![Interval::new(1, 5, "a"), Interval::new(6, 10, "b")];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(11, 15);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_find_contained_no_containment() {
        let intervals = vec![Interval::new(1, 5, "a"), Interval::new(6, 10, "b")];
        let t = IntervalTree::new(intervals);

        let result = t.find_contained(0, 7);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].value, "a");

        let result = t.find_contained(2, 12);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].value, "b");

        let result = t.find_contained(0, 15);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_duplicate_intervals() {
        let intervals = vec![Interval::new(1, 5, "a"), Interval::new(1, 5, "b")];
        let t = IntervalTree::new(intervals);

        let result = t.find_overlapping(3, 3);
        assert_eq!(result.len(), 2);
        let values: Vec<_> = result.iter().map(|i| i.value).collect();
        assert!(values.contains(&"a"));
        assert!(values.contains(&"b"));
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
    fn test_visit_near() {
        let intervals = vec![
            Interval::new(1, 5, "a"),
            Interval::new(6, 10, "b"),
            Interval::new(11, 15, "c"),
        ];
        let t = IntervalTree::new(intervals);

        let mut visited_values = Vec::new();
        t.visit_near(4, 12, &mut |interval| {
            visited_values.push(interval.value);
        });

        assert_eq!(visited_values.len(), 3);
        assert!(visited_values.contains(&"a"));
        assert!(visited_values.contains(&"b"));
        assert!(visited_values.contains(&"c"));
    }

    #[test]
    fn test_visit_overlapping_edge_cases() {
        let intervals = vec![
            Interval::new(1, 5, "a"),
            Interval::new(5, 10, "b"),
            Interval::new(10, 15, "c"),
        ];
        let t = IntervalTree::new(intervals);

        let mut collected = Vec::new();
        t.visit_overlapping(5, 10, &mut |interval| {
            collected.push(interval.value);
        });
        assert_eq!(collected.len(), 3);
        assert!(collected.contains(&"a"));
        assert!(collected.contains(&"b"));
        assert!(collected.contains(&"c"));
    }

    #[test]
    fn test_visit_near_outside_range() {
        let intervals = vec![Interval::new(1, 3, "a"), Interval::new(7, 9, "b")];
        let t = IntervalTree::new(intervals);

        let mut visited = Vec::new();
        t.visit_near(4, 6, &mut |interval| {
            visited.push(interval.value);
        });
        assert!(visited.is_empty());
    }

    #[test]
    fn test_singleton_tree() {
        let intervals = vec![Interval::new(1usize, 3usize, 5.5)];
        let t = IntervalTree::new(intervals);

        // Point query on left
        {
            let v = t.find_overlapping(1, 1);
            assert_eq!(v.len(), 1);
            assert_eq!(v[0].start, 1);
            assert_eq!(v[0].stop, 3);
            assert_eq!(v[0].value, 5.5);
        }

        // Point query in middle
        {
            let v = t.find_overlapping(2, 2);
            assert_eq!(v.len(), 1);
            assert_eq!(v[0].start, 1);
            assert_eq!(v[0].stop, 3);
            assert_eq!(v[0].value, 5.5);
        }

        // Point query on right
        {
            let v = t.find_overlapping(3, 3);
            assert_eq!(v.len(), 1);
            assert_eq!(v[0].start, 1);
            assert_eq!(v[0].stop, 3);
            assert_eq!(v[0].value, 5.5);
        }

        // Non-overlapping queries
        let v = t.find_overlapping(4, 4);
        assert_eq!(v.len(), 0);

        let v = t.find_overlapping(0, 0);
        assert_eq!(v.len(), 0);
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

    #[test]
    fn test_two_identical_intervals_with_different_contents() {
        let intervals = vec![
            Interval::new(5usize, 10usize, 10.5),
            Interval::new(5usize, 10usize, 5.5),
        ];
        let t = IntervalTree::new(intervals);

        let v = t.find_overlapping(6, 6);
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].start, 5);
        assert_eq!(v[0].stop, 10);
        assert_eq!(v[1].start, 5);
        assert_eq!(v[1].stop, 10);

        let expected = vec![5.5, 10.5];
        let mut actual: Vec<_> = v.iter().map(|i| i.value).collect();
        actual.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut expected_sorted = expected.clone();
        expected_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(actual, expected_sorted);
    }

    #[test]
    fn test_find_containing() {
        let intervals = vec![
            Interval::new(0, 10, "A"),
            Interval::new(5, 15, "B"),
            Interval::new(10, 20, "C"),
            Interval::new(15, 25, "D"),
            Interval::new(20, 30, "E"),
            Interval::new(0, 30, "F"),
        ];

        let tree = IntervalTree::new(intervals);

        let result = tree.find_containing(10, 10);
        assert_eq!(result.len(), 4);
        let expected: HashSet<_> = ["A", "B", "C", "F"].into();
        assert!(result.iter().all(|i| expected.contains(&i.value)));

        let result = tree.find_containing(12, 18);
        assert_eq!(result.len(), 2);
        let expected: HashSet<_> = ["C", "F"].into();
        assert!(result.iter().all(|i| expected.contains(&i.value)));

        let result = tree.find_containing(5, 25);
        assert_eq!(result.len(), 1);
        let expected: HashSet<_> = ["F"].into();
        assert!(result.iter().all(|i| expected.contains(&i.value)));

        let result = tree.find_containing(15, 25);
        assert_eq!(result.len(), 2);
        let expected: HashSet<_> = ["D", "F"].into();
        assert!(result.iter().all(|i| expected.contains(&i.value)));

        let result = tree.find_containing(12, 31);
        assert_eq!(result.len(), 0);

        let result = tree.find_containing(-1, 28);
        assert_eq!(result.len(), 0);

        let result = tree.find_containing(0, 30);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].value, "F");

        let result = tree.find_containing(40, 50);
        assert_eq!(result.len(), 0);
    }

    #[test]
    #[ignore]
    fn brute_force() {
        const N_INTERVALS: usize = 10_000;
        const N_QUERIES: usize = 5000;

        let mut intervals = Vec::new();
        for _ in 0..N_INTERVALS {
            intervals.push(random_interval(0, 1_000_000, 20, 2000, 1));
        }

        let mut queries = Vec::new();
        for _ in 0..N_QUERIES {
            queries.push(random_interval(0, 1_000_000, 20, 2000, 1));
        }

        let mut bruteforcecounts = Vec::new();
        let t0 = Instant::now();
        for q in &queries {
            let mut results = Vec::new();
            for i in &intervals {
                if i.start >= q.start && i.stop <= q.stop {
                    results.push(i.clone());
                }
            }
            bruteforcecounts.push(results.len());
        }
        let t1 = Instant::now();
        let duration = t1 - t0;
        println!("Brute force:\t{:?}", duration);

        let tree = IntervalTree::new(intervals.clone());
        let mut treecounts = Vec::new();
        let t0 = Instant::now();
        for q in &queries {
            let results = tree.find_contained(q.start, q.stop);
            treecounts.push(results.len());
        }
        let t1 = Instant::now();
        let duration = t1 - t0;
        println!("Interval tree:\t{:?}", duration);
        for (b, t) in bruteforcecounts.iter().zip(treecounts.iter()) {
            assert_eq!(b, t);
        }
        println!("All counts match!");
    }
}
