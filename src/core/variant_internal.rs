use super::{
    aligner::WFAligner,
    containers::interval_tree::{Interval, IntervalTree},
    svtype::SvType,
};
use crate::{
    cli::MergeArgsInner,
    io::bed_reader::TrId,
    utils::util::{Result, MISSING_INTEGER, VECTOR_END_INTEGER},
    CONFIG, DISTANCE_OFFSET,
};
use anyhow::anyhow;
use rust_htslib::bcf::{
    self,
    record::{Genotype, GenotypeAllele},
};
use std::{collections::HashSet, fmt, mem};

// TODO: Refactor
#[derive(Debug, Clone)]
pub struct VariantInternal {
    pub start: f32,
    pub end: f32, // abs(svlen)
    pub interval: Option<[f32; 2]>,
    pub svlen: f32,
    pub vcf_id: usize,
    pub svtype: SvType,
    pub index: usize,
    pub max_dist: f32,
    pub min_seqid: f32,
    pub trid: Option<TrId>,
    pub sequence: Option<Vec<u8>>,
    pub sample_id: usize, // TODO: This is not used for now, currently assumes we have 1 sample per VCF
    pub id: String,
    pub is_cn: bool,
    // NOTE: All fields specific for writing accessed from htslib
    pub rid: Option<u32>,
    pub pos: i64,
    pub alleles: Vec<Vec<u8>>,
    pub gt: Genotype,
    pub gq: i32,
    pub pl: Vec<i32>,
    pub ad: Vec<i32>,
}

// TODO: Currently we assume a variant corresponds to a single sample
impl VariantInternal {
    pub fn from_vcf_record(
        record: &bcf::Record,
        vcf_id: usize,
        args: &MergeArgsInner,
        tr_it: &Option<&IntervalTree<u32, TrId>>,
    ) -> Result<Self> {
        let svtype = match record.info(b"SVTYPE").string().unwrap_or(None) {
            Some(svtype) => SvType::from_u8(svtype[0])?,
            None => return Err(anyhow!("SVTYPE missing in VCF record")),
        };

        let sequence = if svtype == SvType::INSERTION {
            match record.alleles() {
                alleles if alleles.len() > 1 => Some(alleles[1].to_vec()),
                _ => None,
            }
        } else {
            None
        };

        let svlen = if svtype == SvType::BND {
            0.0
        } else {
            match record.info(b"SVLEN").integer().unwrap_or(None) {
                Some(svlen) => svlen[0] as f32,
                None => return Err(anyhow!("SVLEN missing in VCF record")),
            }
        };

        let start = record.pos() as f32;
        let end = svlen.abs();

        let raw_id_bytes = record.id();
        let vcf_id_str = std::str::from_utf8(&raw_id_bytes).unwrap_or_default();

        let mut is_cn = false;
        let parts: Vec<&str> = vcf_id_str.split(':').collect();
        if parts.len() > 1 && parts[1] == "CNV" {
            is_cn = true;
        }

        let id = format!("{}_", vcf_id) + vcf_id_str;

        let mut max_dist = args.max_dist;
        let min_seqid = args.min_sequence_similarity;

        // there may be a per-variant distance threshold (i.e., a VCF is given that comes from svx)
        // Alternatively there may be a per-sample distance threshold
        // max_dist = PER_SAMPLE_DISTS[sample]

        // Alternatively, there is a length-based threshold
        if args.use_linear_threshold && args.max_dist_linear > 0.0 {
            max_dist = (args.max_dist_linear * svlen.abs()) as i32;
            // TODO: Some code that will still use max_dist but ONLY if it was explicitly set
            if args.min_dist != -1 {
                max_dist = max_dist.max(args.min_dist);
            }
        }

        let interval = match svtype {
            SvType::DELETION | SvType::INVERSION | SvType::DUPLICATION => {
                if args.min_recip_overlap > 0.0 {
                    let end_pos = match record.info(b"END").integer().unwrap_or(None) {
                        Some(end) if !end.is_empty() => end[0] as f32,
                        _ => start + svlen.abs(),
                    };
                    Some([start, end_pos])
                } else {
                    None
                }
            }
            _ => None,
        };

        // htslib related extraction
        // TODO: Refactor
        let alleles: Vec<Vec<u8>> = record.alleles().iter().map(|s| s.to_vec()).collect();

        let gt_buffer = record.genotypes().expect("Error reading genotypes!");
        let gt = gt_buffer.get(0);

        let first_gq: i32 = match record.format(b"GQ").integer() {
            Ok(buffer_backed) => buffer_backed
                .iter()
                .next()
                .and_then(|slice| slice.first())
                .copied()
                .unwrap_or(0),
            Err(_) => 0,
        };

        let all_pls: Vec<i32> = match record.format(b"PL").integer() {
            Ok(pl_values_per_sample) => match pl_values_per_sample.first() {
                Some(pl_slice) if !pl_slice.is_empty() => pl_slice.to_vec(),
                _ => vec![0],
            },
            Err(_) => vec![0],
        };

        let all_ads: Vec<i32> = match record.format(b"AD").integer() {
            Ok(ad_values_per_sample) => {
                match ad_values_per_sample.first() {
                    Some(ad_slice) if !ad_slice.is_empty() => {
                        if ad_slice.len() == 1 && ad_slice[0] == MISSING_INTEGER {
                            // Special case for single MISSING_INTEGER
                            vec![MISSING_INTEGER, VECTOR_END_INTEGER]
                        } else {
                            ad_slice.to_vec()
                        }
                    }
                    _ => vec![MISSING_INTEGER, VECTOR_END_INTEGER],
                }
            }
            Err(_) => vec![MISSING_INTEGER, VECTOR_END_INTEGER],
        };

        let trid = annotate_tr_containment(start, svlen, tr_it, args.tr_span_ratio_threshold);

        Ok(Self {
            start,
            end,
            interval,
            svlen,
            vcf_id,
            sample_id: vcf_id,
            svtype,
            id,
            index: 0,
            min_seqid,
            max_dist: max_dist as f32,
            trid,
            sequence,
            rid: record.rid(),
            pos: record.pos(),
            alleles,
            gt,
            gq: first_gq,
            pl: all_pls,
            ad: all_ads,
            is_cn,
        })
    }

    pub fn from_parts(
        vcf_id: usize,
        id: String,
        svtype: SvType,
        start: f32,
        end: f32,
    ) -> Result<Self> {
        let is_cn = false;
        let alleles: Vec<GenotypeAllele> = Vec::new();
        let g: Genotype = unsafe { mem::transmute(alleles) };
        let interval = match svtype {
            SvType::DELETION | SvType::INVERSION | SvType::DUPLICATION => Some([start, end]),
            _ => None,
        };

        Ok(Self {
            start,
            end,
            interval,
            svlen: end,
            vcf_id,
            sample_id: vcf_id,
            svtype,
            id,
            index: 0,
            min_seqid: 0.0,
            max_dist: 1000.0,
            trid: None,
            sequence: None,
            rid: None,
            pos: 0,
            alleles: Vec::new(),
            gt: g,
            gq: 0,
            pl: Vec::new(),
            ad: Vec::new(),
            is_cn,
        })
    }

    pub fn point(&self) -> [f64; 2] {
        [self.start as f64, self.end as f64]
        // [(self.start as f64).ln(), (self.end as f64).ln()]
    }

    pub fn point_distance(p1: &[f64; 2], p2: &[f64; 2]) -> f32 {
        let d_start = p1[0] - p2[0];
        let d_end = p1[1] - p2[1];
        let kd_tree_norm = CONFIG.kd_tree_norm;
        if kd_tree_norm == 2 {
            ((d_start * d_start + d_end * d_end) as f32).sqrt()
        } else {
            let pow_sum =
                (d_start.abs().powi(kd_tree_norm) + d_end.abs().powi(kd_tree_norm)) as f32;
            pow_sum.powf(1.0 / kd_tree_norm as f32)
        }
    }

    #[inline]
    pub fn distance(&self, other: &VariantInternal) -> f32 {
        let d_start = self.start - other.start;
        let d_end = self.end - other.end;
        let kd_tree_norm = CONFIG.kd_tree_norm;
        if kd_tree_norm == 2 {
            (d_start * d_start + d_end * d_end).sqrt()
        } else {
            let pow_sum = d_start.abs().powi(kd_tree_norm) + d_end.abs().powi(kd_tree_norm);
            pow_sum.powf(1.0 / kd_tree_norm as f32)
        }
    }

    pub fn passes_seq_similarity(
        &self,
        other: &VariantInternal,
        aligner: &mut WFAligner,
        required_similarity: f32,
    ) -> bool {
        if self.sequence.is_none() || other.sequence.is_none() {
            return true;
        }

        let s = self.sequence.as_ref().unwrap();
        let t = other.sequence.as_ref().unwrap();

        let s_len = s.len();
        let t_len = t.len();

        let min_length = s_len.min(t_len);
        let max_length = s_len + t_len - min_length;

        if max_length == 0 {
            return required_similarity <= 0.0 + DISTANCE_OFFSET;
        }

        // if the shorter sequence is too short relative to the longer one to ever meet the required_similarity, even with zero edits.
        // Max possible similarity = min_length / max_length
        if (min_length as f32 / max_length as f32) < required_similarity - DISTANCE_OFFSET {
            return false;
        }

        let _status = aligner.align_end_to_end(s, t);
        let edit_distance = aligner.score() as f32;
        let dist = 1.0 - (edit_distance / max_length as f32);

        dist >= required_similarity - DISTANCE_OFFSET
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

        (min_end - max_start + DISTANCE_OFFSET) >= max_interval_size * min_recip_overlap
    }

    pub fn passes_kmer_jaccard_similarity(
        &self,
        other: &VariantInternal,
        required_jaccard_similarity: f32,
    ) -> bool {
        if self.sequence.is_none() || other.sequence.is_none() {
            return true;
        }

        let s_seq = self.sequence.as_ref().unwrap();
        let t_seq = other.sequence.as_ref().unwrap();

        // let trid_self = self.trid.as_ref().expect("Context: self.trid must be Some");
        // let k = trid_self.motif_len;
        let k = 9;
        if s_seq.len() < k || t_seq.len() < k {
            // If either sequence is shorter than k, Jaccard index is 0 unless both are empty
            // and k is also 0 (which is guarded by motif_len > 0).
            // If required_jaccard_similarity is 0.0, this would pass.
            return 0.0 >= required_jaccard_similarity - DISTANCE_OFFSET;
        }

        let s_kmers = generate_kmers(s_seq, k);
        let t_kmers = generate_kmers(t_seq, k);

        let ji = jaccard_index(&s_kmers, &t_kmers);

        ji >= required_jaccard_similarity - DISTANCE_OFFSET
    }
}

fn generate_kmers(sequence: &[u8], k: usize) -> HashSet<Vec<u8>> {
    let mut kmers = HashSet::new();
    if k == 0 || sequence.len() < k {
        return kmers;
    }
    for i in 0..=(sequence.len() - k) {
        kmers.insert(sequence[i..i + k].to_vec());
    }
    kmers
}

fn jaccard_index(set1: &HashSet<Vec<u8>>, set2: &HashSet<Vec<u8>>) -> f32 {
    if set1.is_empty() && set2.is_empty() {
        return 1.0;
    }
    if set1.is_empty() || set2.is_empty() {
        return 0.0;
    }

    let intersection_size = set1.intersection(set2).count() as f32;
    let union_size = (set1.len() + set2.len()) as f32 - intersection_size;

    if union_size == 0.0 {
        if intersection_size == 0.0 {
            1.0
        } else {
            0.0
        }
    } else {
        intersection_size / union_size
    }
}

impl fmt::Display for VariantInternal {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "(id: {}, sample: {}, start: {}, end: {})",
            self.id, self.sample_id, self.start, self.end
        )
    }
}

fn annotate_tr_containment(
    start: f32,
    svlen: f32,
    tr_it: &Option<&IntervalTree<u32, TrId>>,
    tr_span_ratio_threshold: f32,
) -> Option<TrId> {
    let mut trid = None;

    if let Some(interval_tree) = tr_it {
        let overlaps = interval_tree.find_containing(start as u32, start as u32);
        let mut best_tr_candidate: Option<&Interval<u32, TrId>> = None;
        let mut min_start_distance: u32 = u32::MAX;

        if !overlaps.is_empty() {
            for tr_interval in &overlaps {
                let distance = (start as u32).abs_diff(tr_interval.start);
                if distance < min_start_distance {
                    min_start_distance = distance;
                    best_tr_candidate = Some(tr_interval);
                }
            }
        }

        if let Some(chosen_tr) = best_tr_candidate {
            let tr_span = chosen_tr.stop.saturating_sub(chosen_tr.start);
            let sv_span_abs = svlen.abs();

            // Some other ideas:
            //   motif matching fraction, i.e., match sequence with HMM of TR
            //   anchor uncertainty, i.e., allow some variability in the TR definition start/end this needs to happen in the interval tree
            // Length-ratio
            if tr_span > 0 && sv_span_abs / (tr_span as f32) <= 1.0 + tr_span_ratio_threshold {
                trid = Some(chosen_tr.value.clone());
            } else {
                trid = None;
            }
        } else {
            trid = None;
        }
    }
    trid
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{init_config, utils::util::init_logger, SvxConfig};

    #[test]
    #[ignore]
    fn test_passes_overlap() {
        init_logger();
        init_config(SvxConfig {
            kd_tree_norm: 2,
            dump: false,
        });
        let min_recip_overlap = 0.4f32;

        let var1 =
            VariantInternal::from_parts(0, "var1".to_string(), SvType::DELETION, 100.0, 200.0)
                .unwrap();

        // 1: Basic overlap
        let var2 =
            VariantInternal::from_parts(0, "var2".to_string(), SvType::DELETION, 150.0, 250.0)
                .unwrap();

        // 2: No overlap
        let var3 =
            VariantInternal::from_parts(0, "var3".to_string(), SvType::DELETION, 300.0, 400.0)
                .unwrap();

        // 3: Exact same interval
        let var4 =
            VariantInternal::from_parts(0, "var4".to_string(), SvType::DELETION, 100.0, 200.0)
                .unwrap();

        // 4: One variant without interval (e.g. INSERTION)
        let var5 =
            VariantInternal::from_parts(0, "var5".to_string(), SvType::INSERTION, 500.0, 550.0)
                .unwrap();

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
    #[ignore]
    fn test_min_recip_overlap_threshold() {
        init_logger();

        init_config(SvxConfig {
            kd_tree_norm: 2,
            dump: false,
        });

        // Test case with 50% overlap:
        // var1: 100-200 (length 100)
        // var2: 150-250 (length 100)
        // overlap region: 150-200 (length 50)
        // overlap ratio = 50/100 = 0.5
        let var1 =
            VariantInternal::from_parts(0, "var1".to_string(), SvType::DELETION, 100.0, 200.0)
                .unwrap();
        let var2 =
            VariantInternal::from_parts(0, "var2".to_string(), SvType::DELETION, 150.0, 250.0)
                .unwrap();

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
