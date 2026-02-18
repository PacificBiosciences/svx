use super::VariantInternal;
use crate::{
    DISTANCE_OFFSET,
    core::{aligner::WFAligner, svtype::SvType},
};

impl VariantInternal {
    pub fn point(&self) -> [f64; 2] {
        [self.start, self.end]
    }

    // pub fn point_distance(p1: &[f64; 2], p2: &[f64; 2]) -> f64 {
    //     let d_start = p1[0] - p2[0];
    //     let d_end = p1[1] - p2[1];
    //     (d_start * d_start + d_end * d_end).sqrt()
    // }

    #[inline]
    pub fn distance(&self, other: &VariantInternal) -> f64 {
        let d_start = self.start - other.start;
        let d_end = self.end - other.end;
        (d_start * d_start + d_end * d_end).sqrt()
    }

    pub fn passes_seq_similarity(
        &self,
        other: &VariantInternal,
        aligner: &mut WFAligner,
        required_similarity: f32,
    ) -> bool {
        let (Some(s), Some(t)) = (self.sequence.as_deref(), other.sequence.as_deref()) else {
            return true;
        };

        if seq_similarity_passes(s, t, aligner, required_similarity) {
            return true;
        }

        let Some((start_delta, motif_len)) = tr_unroll_params(self, other) else {
            return false;
        };

        seq_similarity_passes_with_unroll(
            s,
            t,
            aligner,
            required_similarity,
            start_delta,
            motif_len,
        )
    }

    #[inline]
    pub fn size_similarity(&self, other: &VariantInternal) -> f64 {
        let mut size_a = self.svlen.abs();
        let mut size_b = other.svlen.abs();
        if size_a <= DISTANCE_OFFSET || size_b <= DISTANCE_OFFSET {
            if (size_a - size_b).abs() <= DISTANCE_OFFSET {
                return 1.0;
            }
            if size_a <= DISTANCE_OFFSET {
                size_a = 1.0;
            }
            if size_b <= DISTANCE_OFFSET {
                size_b = 1.0;
            }
        }

        size_a.min(size_b) / size_a.max(size_b)
    }

    #[inline]
    pub fn passes_size_similarity(
        &self,
        other: &VariantInternal,
        min_size_similarity: f64,
    ) -> bool {
        if min_size_similarity <= DISTANCE_OFFSET {
            return true;
        }
        self.size_similarity(other) >= min_size_similarity - DISTANCE_OFFSET
    }

    #[inline]
    pub fn passes_overlap(&self, other: &VariantInternal, min_recip_overlap: f32) -> bool {
        if self.interval.is_none() || other.interval.is_none() {
            return true;
        }

        let self_interval = self.interval.unwrap();
        let other_interval = other.interval.unwrap();

        let max_start = self_interval[0].max(other_interval[0]);
        let min_end = self_interval[1].min(other_interval[1]);

        if min_end <= max_start + DISTANCE_OFFSET {
            return false;
        }

        let max_interval_size =
            (self_interval[1] - self_interval[0]).max(other_interval[1] - other_interval[0]);

        (min_end - max_start + DISTANCE_OFFSET) >= max_interval_size * f64::from(min_recip_overlap)
    }
}

fn seq_similarity_passes(
    s: &[u8],
    t: &[u8],
    aligner: &mut WFAligner,
    required_similarity: f32,
) -> bool {
    let distance_offset = DISTANCE_OFFSET as f32;
    let min_length = s.len().min(t.len());
    let max_length = s.len().max(t.len());

    if max_length == 0 {
        return required_similarity <= distance_offset;
    }

    // If the shorter sequence is too short relative to the longer one to ever meet the
    // required_similarity (even with zero edits), reject without aligning.
    // Max possible similarity = min_length / max_length
    if (min_length as f32 / max_length as f32) < required_similarity - distance_offset {
        return false;
    }

    let _status = aligner.align_end_to_end(s, t);
    let edit_distance = aligner.score() as f32;
    let dist = 1.0 - (edit_distance / max_length as f32);

    dist >= required_similarity - distance_offset
}

fn tr_unroll_params(this: &VariantInternal, other: &VariantInternal) -> Option<(i64, usize)> {
    if this.svtype != SvType::INSERTION || other.svtype != SvType::INSERTION {
        return None;
    }

    let (Some(this_tr), Some(other_tr)) = (&this.trid, &other.trid) else {
        return None;
    };
    if this_tr.id != other_tr.id {
        return None;
    }

    let motif_len = this_tr.motif_len.min(other_tr.motif_len);
    if motif_len == 0 {
        return None;
    }

    let start_delta = (this.start - other.start).round() as i64;
    Some((start_delta, motif_len))
}

fn seq_similarity_passes_with_unroll(
    s: &[u8],
    t: &[u8],
    aligner: &mut WFAligner,
    required_similarity: f32,
    start_delta: i64,
    motif_len: usize,
) -> bool {
    let mut shifts = vec![start_delta, -start_delta];
    let phase_delta = if motif_len > 1 {
        i64::try_from(motif_len)
            .ok()
            .map(|motif_len_i64| start_delta.rem_euclid(motif_len_i64))
    } else {
        None
    };
    if let Some(phase_delta) = phase_delta {
        shifts.push(phase_delta);
        shifts.push(-phase_delta);
    }
    shifts.sort_unstable();
    shifts.dedup();

    for shift in shifts {
        if shift == 0 {
            continue;
        }

        let rotated_t = rotate_right(t, shift);
        if seq_similarity_passes(s, &rotated_t, aligner, required_similarity) {
            return true;
        }

        let rotated_s = rotate_right(s, shift);
        if seq_similarity_passes(&rotated_s, t, aligner, required_similarity) {
            return true;
        }
    }

    false
}

fn rotate_right(seq: &[u8], shift: i64) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }

    let len = seq.len();
    let Ok(len_i64) = i64::try_from(len) else {
        return seq.to_vec();
    };
    let normalized_shift = shift.rem_euclid(len_i64) as usize;
    if normalized_shift == 0 {
        return seq.to_vec();
    }

    let split = len - normalized_shift;
    let mut rotated = Vec::with_capacity(len);
    rotated.extend_from_slice(&seq[split..]);
    rotated.extend_from_slice(&seq[..split]);
    rotated
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::aligner::{AlignmentScope, MemoryModel};
    use crate::core::svtype::SvType;
    use crate::core::variant::test_utils;
    use crate::io::bed_reader::TrId;

    fn new_ed_aligner() -> WFAligner {
        WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
            .edit()
            .build()
    }

    #[test]
    fn seq_similarity_empty_sequences() {
        let mut aligner = new_ed_aligner();
        assert!(seq_similarity_passes(b"", b"", &mut aligner, 0.0));
        assert!(!seq_similarity_passes(b"", b"", &mut aligner, 0.1));
    }

    #[test]
    fn seq_similarity_identical_sequences() {
        let mut aligner = new_ed_aligner();
        assert!(seq_similarity_passes(
            b"ACGTACGT",
            b"ACGTACGT",
            &mut aligner,
            1.0
        ));
    }

    #[test]
    fn seq_similarity_length_ratio_prefilter_rejects_impossible_pairs() {
        let mut aligner = new_ed_aligner();
        assert!(!seq_similarity_passes(
            b"AAAAAAAAAA",
            b"AAAAAAAAAAAAAAAAAAAA",
            &mut aligner,
            0.60
        ));
    }

    #[test]
    fn tr_shared_id_allows_unrolled_sequence_similarity_for_insertions() {
        let mut aligner = new_ed_aligner();

        let mut a = test_utils::from_parts(0, "a".to_string(), SvType::INSERTION, 100.0, 0.0)
            .expect("variant should build");
        let mut b = test_utils::from_parts(1, "b".to_string(), SvType::INSERTION, 101.0, 0.0)
            .expect("variant should build");

        a.sequence = Some(b"ACGTACGT".to_vec());
        b.sequence = Some(b"CGTACGTA".to_vec());
        a.trid = Some(TrId {
            id: "TR1".to_string(),
            motif_len: 4,
        });
        b.trid = Some(TrId {
            id: "TR1".to_string(),
            motif_len: 4,
        });

        assert!(a.passes_seq_similarity(&b, &mut aligner, 0.95));
    }

    #[test]
    fn tr_unroll_requires_shared_tr_id() {
        let mut aligner = new_ed_aligner();

        let mut a = test_utils::from_parts(0, "a".to_string(), SvType::INSERTION, 100.0, 0.0)
            .expect("variant should build");
        let mut b = test_utils::from_parts(1, "b".to_string(), SvType::INSERTION, 101.0, 0.0)
            .expect("variant should build");

        a.sequence = Some(b"ACGTACGT".to_vec());
        b.sequence = Some(b"CGTACGTA".to_vec());
        a.trid = Some(TrId {
            id: "TR1".to_string(),
            motif_len: 4,
        });
        b.trid = Some(TrId {
            id: "TR2".to_string(),
            motif_len: 4,
        });

        assert!(!a.passes_seq_similarity(&b, &mut aligner, 0.95));
    }

    #[test]
    fn test_distance_preserves_1bp_resolution_at_large_coordinates() {
        let start1 = 3_000_000_000.0f64;
        let start2 = start1 + 1.0f64;
        assert_ne!(
            start1, start2,
            "f64 must represent 1bp deltas at this scale"
        );

        let var1 = test_utils::from_parts(0, "var1".to_string(), SvType::INSERTION, start1, 100.0)
            .unwrap();
        let var2 = test_utils::from_parts(0, "var2".to_string(), SvType::INSERTION, start2, 100.0)
            .unwrap();

        assert!(
            var1.distance(&var2) >= 1.0,
            "Distance should reflect a 1bp delta"
        );
    }

    #[test]
    fn test_passes_overlap() {
        let min_recip_overlap = 0.4f32;

        let var1 =
            test_utils::from_parts(0, "var1".to_string(), SvType::DELETION, 100.0, 200.0).unwrap();

        // 1: Basic overlap
        let var2 =
            test_utils::from_parts(0, "var2".to_string(), SvType::DELETION, 150.0, 250.0).unwrap();

        // 2: No overlap
        let var3 =
            test_utils::from_parts(0, "var3".to_string(), SvType::DELETION, 300.0, 400.0).unwrap();

        // 3: Exact same interval
        let var4 =
            test_utils::from_parts(0, "var4".to_string(), SvType::DELETION, 100.0, 200.0).unwrap();

        // 4: One variant without interval (e.g. INSERTION)
        let var5 =
            test_utils::from_parts(0, "var5".to_string(), SvType::INSERTION, 500.0, 550.0).unwrap();

        assert!(
            var1.passes_overlap(&var2, min_recip_overlap),
            "Overlapping variants should pass"
        );
        assert!(
            !var1.passes_overlap(&var3, min_recip_overlap),
            "Non-overlapping variants should not pass"
        );
        assert!(
            var1.passes_overlap(&var4, min_recip_overlap),
            "Same interval variants should pass"
        );
        assert!(
            var1.passes_overlap(&var5, min_recip_overlap),
            "Variant without interval should always pass"
        );
        assert!(
            var5.passes_overlap(&var1, min_recip_overlap),
            "Variant without interval should always pass (reverse)"
        );
    }

    #[test]
    fn test_min_recip_overlap_threshold() {
        // Test case with 50% overlap:
        // var1: 100-200 (length 100)
        // var2: 150-250 (length 100)
        // overlap region: 150-200 (length 50)
        // overlap ratio = 50/100 = 0.5
        let var1 =
            test_utils::from_parts(0, "var1".to_string(), SvType::DELETION, 100.0, 200.0).unwrap();
        let var2 =
            test_utils::from_parts(0, "var2".to_string(), SvType::DELETION, 150.0, 250.0).unwrap();

        let test_thresholds = vec![
            (0.3, true),  // Should pass: 50% overlap > 30% required
            (0.5, true),  // Should pass: 50% overlap = 50% required
            (0.6, false), // Should fail: 50% overlap < 60% required
            (0.8, false), // Should fail: 50% overlap < 80% required
        ];

        for (threshold, expected_result) in test_thresholds {
            assert_eq!(
                var1.passes_overlap(&var2, threshold),
                expected_result,
                "Failed with min_recip_overlap = {}: expected {}, got {}",
                threshold,
                expected_result,
                !expected_result
            );
        }
    }
}
