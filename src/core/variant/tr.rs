use crate::core::containers::interval_tree::{Interval, IntervalTree};
use crate::io::bed_reader::TrId;

pub const TR_CONTAINMENT_SCALE: u32 = 1_000_000;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrContainmentQuery {
    Deletion { start: u32, end: u32 },
    Insertion { pos: u32 },
}

#[derive(Clone, Copy, Debug)]
pub struct TrContainmentConfig {
    pub span_query_slop: u32,
    pub min_span_containment_scaled: u32,
    pub min_span_overlap_bp: u32,
    pub ins_max_dist: u32,
}

pub fn annotate_tr_containment(
    query: TrContainmentQuery,
    tr_it: &Option<&IntervalTree<u32, TrId>>,
    cfg: TrContainmentConfig,
) -> Option<TrId> {
    let Some(interval_tree) = tr_it else {
        return None;
    };

    match query {
        TrContainmentQuery::Deletion { start, end } => {
            annotate_deletion(start, end, interval_tree, cfg)
        }
        TrContainmentQuery::Insertion { pos } => annotate_insertion(pos, interval_tree, cfg),
    }
}

#[inline]
fn annotate_deletion(
    sv_start: u32,
    sv_end: u32,
    interval_tree: &IntervalTree<u32, TrId>,
    cfg: TrContainmentConfig,
) -> Option<TrId> {
    if sv_end <= sv_start {
        return None;
    }

    let sv_len = sv_end - sv_start;
    let q_start = sv_start.saturating_sub(cfg.span_query_slop);
    let q_end = sv_end.saturating_add(cfg.span_query_slop);

    let mut best_id: Option<TrId> = None;
    let mut best_dice: u32 = 0;
    let mut best_overlap: u32 = 0;
    let mut best_tr_len: u32 = u32::MAX;
    let mut best_tr_start: u32 = u32::MAX;

    interval_tree.visit_overlapping(q_start, q_end, &mut |tr_interval: &Interval<u32, TrId>| {
        if tr_interval.stop <= tr_interval.start {
            return;
        }

        let tr_len = tr_interval.stop - tr_interval.start;
        let overlap = overlap_len(sv_start, sv_end, tr_interval.start, tr_interval.stop);
        if overlap < cfg.min_span_overlap_bp {
            return;
        }
        if !passes_containment(overlap, sv_len, cfg.min_span_containment_scaled) {
            return;
        }

        let dice = dice_scaled(overlap, sv_len, tr_len);
        let better = best_id.is_none()
            || dice > best_dice
            || (dice == best_dice && overlap > best_overlap)
            || (dice == best_dice && overlap == best_overlap && tr_len < best_tr_len)
            || (dice == best_dice
                && overlap == best_overlap
                && tr_len == best_tr_len
                && tr_interval.start < best_tr_start);
        if better {
            best_id = Some(tr_interval.value.clone());
            best_dice = dice;
            best_overlap = overlap;
            best_tr_len = tr_len;
            best_tr_start = tr_interval.start;
        }
    });

    best_id
}

#[inline]
fn annotate_insertion(
    pos: u32,
    interval_tree: &IntervalTree<u32, TrId>,
    cfg: TrContainmentConfig,
) -> Option<TrId> {
    let pad = cfg.ins_max_dist.saturating_add(1);
    let q_start = pos.saturating_sub(pad);
    let q_end = pos.saturating_add(pad);

    let mut best_id: Option<TrId> = None;
    let mut best_inside = false;
    let mut best_dist = u32::MAX;
    let mut best_tr_len = u32::MAX;
    let mut best_tr_start = u32::MAX;

    interval_tree.visit_overlapping(q_start, q_end, &mut |tr_interval: &Interval<u32, TrId>| {
        if tr_interval.stop <= tr_interval.start {
            return;
        }

        let tr_len = tr_interval.stop - tr_interval.start;
        let d0 = dist_point_to_closed_interval(pos, tr_interval.start, tr_interval.stop);
        let d1 = if pos > 0 {
            dist_point_to_closed_interval(pos - 1, tr_interval.start, tr_interval.stop)
        } else {
            u32::MAX
        };
        let dist = d0.min(d1);
        if dist > cfg.ins_max_dist {
            return;
        }

        let inside = in_interval_by_flank(pos, tr_interval.start, tr_interval.stop);
        let better = best_id.is_none()
            || (inside && !best_inside)
            || (inside == best_inside && dist < best_dist)
            || (inside == best_inside && dist == best_dist && tr_len < best_tr_len)
            || (inside == best_inside
                && dist == best_dist
                && tr_len == best_tr_len
                && tr_interval.start < best_tr_start);
        if better {
            best_id = Some(tr_interval.value.clone());
            best_inside = inside;
            best_dist = dist;
            best_tr_len = tr_len;
            best_tr_start = tr_interval.start;
        }
    });

    best_id
}

#[inline(always)]
fn overlap_len(a_start: u32, a_end: u32, b_start: u32, b_end: u32) -> u32 {
    let start = a_start.max(b_start);
    let end = a_end.min(b_end);
    end.saturating_sub(start)
}

#[inline(always)]
fn passes_containment(overlap: u32, sv_len: u32, min_containment_scaled: u32) -> bool {
    if sv_len == 0 {
        return false;
    }
    u64::from(overlap) * u64::from(TR_CONTAINMENT_SCALE)
        >= u64::from(min_containment_scaled) * u64::from(sv_len)
}

#[inline(always)]
fn dice_scaled(overlap: u32, sv_len: u32, tr_len: u32) -> u32 {
    let denominator = u64::from(sv_len) + u64::from(tr_len);
    if denominator == 0 {
        return 0;
    }
    ((2_u64 * u64::from(overlap) * u64::from(TR_CONTAINMENT_SCALE)) / denominator) as u32
}

#[inline(always)]
fn dist_point_to_closed_interval(pos: u32, start: u32, end_exclusive: u32) -> u32 {
    if end_exclusive < start {
        return u32::MAX;
    }
    if pos < start {
        return start - pos;
    }
    if pos > end_exclusive {
        return pos - end_exclusive;
    }
    0
}

#[inline(always)]
fn in_interval_by_flank(pos: u32, start: u32, end_exclusive: u32) -> bool {
    let right_flank_inside = pos >= start && pos < end_exclusive;
    let left_flank_inside = pos > 0 && pos > start && pos - 1 < end_exclusive;
    right_flank_inside || left_flank_inside
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_tr_id() -> TrId {
        TrId {
            id: "TR1".to_string(),
            motif_len: 2,
        }
    }

    fn test_tr2_id() -> TrId {
        TrId {
            id: "TR2".to_string(),
            motif_len: 3,
        }
    }

    fn cfg_with_containment(min_span_containment_scaled: u32) -> TrContainmentConfig {
        TrContainmentConfig {
            span_query_slop: 0,
            min_span_containment_scaled,
            min_span_overlap_bp: 1,
            ins_max_dist: 2,
        }
    }

    #[test]
    fn ins_at_tr_end_boundary_is_contained() {
        let interval_tree = IntervalTree::new(vec![Interval::new(10, 20, test_tr_id())]);
        let tr_it = Some(&interval_tree);
        let cfg = TrContainmentConfig {
            ins_max_dist: 0,
            ..cfg_with_containment(800_000)
        };

        let inside =
            annotate_tr_containment(TrContainmentQuery::Insertion { pos: 19 }, &tr_it, cfg);
        assert_eq!(inside, Some(test_tr_id()));

        let at_exclusive_end =
            annotate_tr_containment(TrContainmentQuery::Insertion { pos: 20 }, &tr_it, cfg);
        assert_eq!(at_exclusive_end, Some(test_tr_id()));
    }

    #[test]
    fn del_majority_overlap_is_contained_even_if_start_is_outside_tr() {
        let interval_tree = IntervalTree::new(vec![Interval::new(100, 200, test_tr_id())]);
        let tr_it = Some(&interval_tree);
        let cfg = cfg_with_containment(950_000);

        // DEL-like span [99,151) overlaps TR [100,200) by 51bp out of 52bp.
        let contained = annotate_tr_containment(
            TrContainmentQuery::Deletion {
                start: 99,
                end: 151,
            },
            &tr_it,
            cfg,
        );
        assert_eq!(contained, Some(test_tr_id()));
    }

    #[test]
    fn del_spanning_two_trs_without_majority_in_single_tr_is_not_contained() {
        let interval_tree = IntervalTree::new(vec![
            Interval::new(100, 150, test_tr_id()),
            Interval::new(150, 200, test_tr2_id()),
        ]);
        let tr_it = Some(&interval_tree);
        let cfg = cfg_with_containment(600_000);

        // DEL-like span [120,180) overlaps each TR by 30bp and should not be considered contained.
        let contained = annotate_tr_containment(
            TrContainmentQuery::Deletion {
                start: 120,
                end: 180,
            },
            &tr_it,
            cfg,
        );
        assert!(contained.is_none());
    }

    #[test]
    fn nested_tr_prefers_innermost_for_deletion() {
        let interval_tree = IntervalTree::new(vec![
            Interval::new(100, 200, test_tr_id()),
            Interval::new(130, 160, test_tr2_id()),
        ]);
        let tr_it = Some(&interval_tree);
        let cfg = cfg_with_containment(800_000);

        let contained = annotate_tr_containment(
            TrContainmentQuery::Deletion {
                start: 140,
                end: 150,
            },
            &tr_it,
            cfg,
        );
        assert_eq!(contained, Some(test_tr2_id()));
    }

    #[test]
    fn insertion_proximity_respects_max_distance() {
        let interval_tree = IntervalTree::new(vec![Interval::new(100, 200, test_tr_id())]);
        let tr_it = Some(&interval_tree);
        let cfg = TrContainmentConfig {
            ins_max_dist: 2,
            ..cfg_with_containment(800_000)
        };

        let near_left =
            annotate_tr_containment(TrContainmentQuery::Insertion { pos: 99 }, &tr_it, cfg);
        assert_eq!(near_left, Some(test_tr_id()));

        let near_right =
            annotate_tr_containment(TrContainmentQuery::Insertion { pos: 202 }, &tr_it, cfg);
        assert_eq!(near_right, Some(test_tr_id()));

        let far = annotate_tr_containment(TrContainmentQuery::Insertion { pos: 204 }, &tr_it, cfg);
        assert!(far.is_none());
    }

    #[test]
    fn no_tree_returns_none() {
        let cfg = cfg_with_containment(800_000);
        let got = annotate_tr_containment(TrContainmentQuery::Insertion { pos: 10 }, &None, cfg);
        assert!(got.is_none());
    }
}
