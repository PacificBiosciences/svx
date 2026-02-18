use super::{BndEventData, VariantInternal};
use crate::{cli::MergeArgsInner, utils::util::Result};
use rust_htslib::bcf;

impl VariantInternal {
    fn parse_bnd_alt(alt: &[u8]) -> Result<(String, i64)> {
        let first_bracket_idx = alt
            .iter()
            .position(|&b| b == b'[' || b == b']')
            .ok_or_else(|| {
                crate::svx_error!(
                    "BND ALT is missing '[' or ']': {:?}",
                    String::from_utf8_lossy(alt)
                )
            })?;
        let bracket = alt[first_bracket_idx];
        let second_bracket_idx = alt
            .iter()
            .skip(first_bracket_idx + 1)
            .position(|&b| b == bracket)
            .map(|rel| rel + first_bracket_idx + 1)
            .ok_or_else(|| {
                crate::svx_error!(
                    "BND ALT is missing closing bracket '{}': {:?}",
                    bracket as char,
                    String::from_utf8_lossy(alt)
                )
            })?;

        let mate = &alt[first_bracket_idx + 1..second_bracket_idx];
        let mate = std::str::from_utf8(mate)
            .map_err(|e| crate::svx_error!("Invalid UTF-8 in BND ALT: {e}"))?;
        let (contig, pos) = mate.rsplit_once(':').ok_or_else(|| {
            crate::svx_error!("BND ALT mate does not look like contig:pos: {mate:?}")
        })?;
        let pos_1based: i64 = pos
            .parse()
            .map_err(|e| crate::svx_error!("Invalid BND ALT mate position {pos:?}: {e}"))?;
        if pos_1based <= 0 {
            return Err(crate::svx_error!(
                "Invalid BND ALT mate position (must be >= 1): {pos_1based}"
            ));
        }

        Ok((contig.to_string(), pos_1based - 1))
    }

    fn parse_bnd_strands_from_alt(alt: &[u8]) -> Option<[u8; 2]> {
        let &first = alt.first()?;

        if first == b'[' {
            return Some(*b"+-");
        }
        if first == b']' {
            return Some(*b"--");
        }
        if alt.contains(&b'[') {
            return Some(*b"++");
        }
        if alt.contains(&b']') {
            return Some(*b"-+");
        }

        None
    }

    pub(super) fn parse_bnd_mate_from_record(
        record: &bcf::Record,
        alt: &[u8],
    ) -> Result<(String, i64)> {
        let (alt_chr2, alt_pos0) = Self::parse_bnd_alt(alt)?;

        let chr2 = if let Some(chr2) = Self::parse_info_string(record, b"CHR2")? {
            if chr2 != alt_chr2 {
                return Err(crate::svx_error!(
                    "BND CHR2 ({chr2}) does not match ALT mate contig ({alt_chr2})"
                ));
            }
            chr2
        } else {
            alt_chr2
        };

        let end_1based = if let Some(end) = Self::parse_info_i64(record, b"END")? {
            end
        } else {
            alt_pos0 + 1
        };

        if end_1based <= 0 {
            return Err(crate::svx_error!(
                "Invalid BND END position (must be >= 1): {end_1based}"
            ));
        }
        let end_0based = end_1based - 1;

        if end_0based != alt_pos0 {
            return Err(crate::svx_error!(
                "BND END ({end_1based}) does not match ALT mate position ({})",
                alt_pos0 + 1
            ));
        }

        Ok((chr2, end_0based))
    }

    pub(super) fn parse_bnd_strands_from_record(
        record: &bcf::Record,
        alt: &[u8],
    ) -> Result<[u8; 2]> {
        if let Some(strands) = Self::parse_info_string(record, b"STRANDS")? {
            let strands_b = strands.as_bytes();
            if strands_b.len() != 2 {
                return Err(crate::svx_error!(
                    "BND STRANDS must be 2 characters, got {strands:?}"
                ));
            }
            return Ok([strands_b[0], strands_b[1]]);
        }

        Self::parse_bnd_strands_from_alt(alt).ok_or_else(|| {
            crate::svx_error!(
                "Failed to infer BND strands from ALT in BND record {:?}",
                String::from_utf8_lossy(&record.id())
            )
        })
    }

    pub(super) fn bnd_merge_coords(
        record: &bcf::Record,
        contig: &str,
        mate_contig: &str,
        mate_pos0: i64,
    ) -> Result<(f64, f64)> {
        let pos0 = record.pos();
        let pos1 = pos0 + 1;
        let mate_pos1 = mate_pos0 + 1;

        let avg_start = Self::parse_info_f64(record, b"AVG_START")?;
        let avg_end = Self::parse_info_f64(record, b"AVG_END")?;

        let swap = contig > mate_contig || (contig == mate_contig && pos0 > mate_pos0);

        let (first_1based, second_1based) = if swap {
            (
                avg_end.unwrap_or(mate_pos1 as f64),
                avg_start.unwrap_or(pos1 as f64),
            )
        } else {
            (
                avg_start.unwrap_or(pos1 as f64),
                avg_end.unwrap_or(mate_pos1 as f64),
            )
        };

        Ok((first_1based - 1.0, second_1based - 1.0))
    }

    pub fn bnd_graph_id(&self) -> Result<String> {
        let bnd = self.bnd.as_ref().ok_or_else(|| {
            crate::svx_error!("bnd_graph_id called on non-BND variant {}", self.id)
        })?;
        let (a, b) = if bnd.contig <= bnd.mate_contig {
            (&bnd.contig, &bnd.mate_contig)
        } else {
            (&bnd.mate_contig, &bnd.contig)
        };
        let strands = std::str::from_utf8(&bnd.strands)
            .map_err(|e| crate::svx_error!("BND strands are not valid UTF-8: {e}"))?;
        Ok(format!("{a}_{b}_TRA_{strands}"))
    }

    pub fn bnd_event_graph_id(&self) -> Result<String> {
        let bnd = self.bnd_event.as_ref().ok_or_else(|| {
            crate::svx_error!(
                "bnd_event_graph_id called on variant without BND event payload {}",
                self.id
            )
        })?;
        let a_strands = std::str::from_utf8(&bnd.a_strands)
            .map_err(|e| crate::svx_error!("BND a_strands are not valid UTF-8: {e}"))?;
        let b_strands = std::str::from_utf8(&bnd.b_strands)
            .map_err(|e| crate::svx_error!("BND b_strands are not valid UTF-8: {e}"))?;
        Ok(format!(
            "{}_{}_TRA_{a_strands}_{b_strands}",
            bnd.a_contig, bnd.b_contig
        ))
    }

    pub fn from_bnd_pair(
        v1: VariantInternal,
        v2: VariantInternal,
        args: &MergeArgsInner,
    ) -> Result<Self> {
        use crate::core::svtype::SvType;

        if v1.svtype != SvType::BND || v2.svtype != SvType::BND {
            return Err(crate::svx_error!(
                "from_bnd_pair expects two BND variants, got {} and {}",
                v1.svtype,
                v2.svtype
            ));
        }
        if v1.vcf_id != v2.vcf_id {
            return Err(crate::svx_error!(
                "BND mates must come from the same input VCF, got vcf_id {} and {}",
                v1.vcf_id,
                v2.vcf_id
            ));
        }

        let b1 = v1
            .bnd
            .as_ref()
            .ok_or_else(|| crate::svx_error!("BND variant {} is missing breakend data", v1.id))?;
        let b2 = v2
            .bnd
            .as_ref()
            .ok_or_else(|| crate::svx_error!("BND variant {} is missing breakend data", v2.id))?;

        if b1.mate_id != v2.id {
            return Err(crate::svx_error!(
                "BND mate mismatch: {} has MATEID {}, expected {}",
                v1.id,
                b1.mate_id,
                v2.id
            ));
        }
        if b2.mate_id != v1.id {
            return Err(crate::svx_error!(
                "BND mate mismatch: {} has MATEID {}, expected {}",
                v2.id,
                b2.mate_id,
                v1.id
            ));
        }

        // Canonicalize: a is the lexicographically smaller endpoint by (contig, pos0).
        let v1_key = (&b1.contig, b1.pos0);
        let v2_key = (&b2.contig, b2.pos0);
        let (mut a, mut b) = if v1_key <= v2_key { (v1, v2) } else { (v2, v1) };

        let a_bnd = a.bnd.take().unwrap();
        let b_bnd = b.bnd.take().unwrap();

        // Validate contig/pos cross-consistency (both directions). Some callers (e.g. Sawfish)
        // report microhomology length (HOMLEN) and encode mate coordinates that differ from the
        // mate breakend POS by exactly HOMLEN. Accept this case.
        let a_to_b_diff = (a_bnd.mate_pos - b_bnd.pos0).abs();
        let b_to_a_diff = (b_bnd.mate_pos - a_bnd.pos0).abs();

        if a_bnd.mate_contig != b_bnd.contig
            || !(a_to_b_diff == 0
                || a_bnd.homlen == Some(a_to_b_diff)
                || b_bnd.homlen == Some(a_to_b_diff))
        {
            return Err(crate::svx_error!(
                "BND mate coordinates inconsistent: {} points to {}:{}, but mate is {}:{} (|delta|={}; homlen_a={:?}; homlen_b={:?})",
                a.id,
                a_bnd.mate_contig,
                a_bnd.mate_pos + 1,
                b_bnd.contig,
                b_bnd.pos0 + 1,
                a_to_b_diff,
                a_bnd.homlen,
                b_bnd.homlen
            ));
        }
        if b_bnd.mate_contig != a_bnd.contig
            || !(b_to_a_diff == 0
                || b_bnd.homlen == Some(b_to_a_diff)
                || a_bnd.homlen == Some(b_to_a_diff))
        {
            return Err(crate::svx_error!(
                "BND mate coordinates inconsistent: {} points to {}:{}, but mate is {}:{} (|delta|={}; homlen_a={:?}; homlen_b={:?})",
                b.id,
                b_bnd.mate_contig,
                b_bnd.mate_pos + 1,
                a_bnd.contig,
                a_bnd.pos0 + 1,
                b_to_a_diff,
                a_bnd.homlen,
                b_bnd.homlen
            ));
        }

        let a_vcf = a
            .vcf
            .take()
            .ok_or_else(|| crate::svx_error!("BND variant {} is missing VCF write data", a.id))?;
        let b_vcf = b
            .vcf
            .take()
            .ok_or_else(|| crate::svx_error!("BND variant {} is missing VCF write data", b.id))?;

        let a_sample_gts = if let Some(format_data) = &a_vcf.format_data {
            format_data.sample_gts.clone()
        } else {
            vec![a_vcf.gt.clone()]
        };
        let b_sample_gts = if let Some(format_data) = &b_vcf.format_data {
            format_data.sample_gts.clone()
        } else {
            vec![b_vcf.gt.clone()]
        };
        if a_sample_gts != b_sample_gts {
            return Err(crate::svx_error!(
                "BND mates {} and {} have per-sample GT mismatch",
                a.id,
                b.id
            ));
        }

        let a_ref_base = a_vcf
            .alleles
            .first()
            .and_then(|ref_allele| ref_allele.first())
            .copied()
            .unwrap_or(b'N');
        let b_ref_base = b_vcf
            .alleles
            .first()
            .and_then(|ref_allele| ref_allele.first())
            .copied()
            .unwrap_or(b'N');
        let inv_event_id = if a_bnd.has_inv_breakpoint_filter
            && b_bnd.has_inv_breakpoint_filter
            && a_bnd.inv_event_id.is_some()
            && a_bnd.inv_event_id == b_bnd.inv_event_id
        {
            a_bnd.inv_event_id.clone()
        } else {
            None
        };

        if a.support_calls != b.support_calls {
            return Err(crate::svx_error!(
                "BND mates {} and {} have inconsistent SUPP_CALLS values: {} vs {}",
                a.id,
                b.id,
                a.support_calls,
                b.support_calls
            ));
        }

        let max_dist = f64::from(args.max_dist);
        let id = format!("{}/{}", a.id, b.id);
        let mut id_list = a.id_list.clone();
        for source_id in &b.id_list {
            if !id_list.contains(source_id) {
                id_list.push(source_id.clone());
            }
        }
        let mut support_mask = vec![0u64; a.support_mask.len().max(b.support_mask.len())];
        for (word_idx, support_word) in support_mask.iter_mut().enumerate() {
            let a_word = a.support_mask.get(word_idx).copied().unwrap_or(0);
            let b_word = b.support_mask.get(word_idx).copied().unwrap_or(0);
            *support_word = a_word | b_word;
        }

        Ok(Self {
            start: a.start,
            end: a.end,
            interval: None,
            svlen: 0.0,
            vcf_id: a.vcf_id,
            svtype: SvType::BND,
            index: 0,
            max_dist,
            info_hash: a.info_hash ^ b.info_hash,
            sequence: None,
            support_mask,
            support_calls: a.support_calls,
            sample_id: a.sample_id,
            start_mean: a.start_mean,
            start_variance: a.start_variance,
            svlen_mean: b.start_mean,
            svlen_variance: b.start_variance,
            id,
            id_list,
            svclaim: None,
            svclaims: Vec::new(),
            bnd: None,
            bnd_event: Some(BndEventData {
                a_contig: a_bnd.contig,
                a_pos0: a_bnd.pos0,
                a_strands: a_bnd.strands,
                a_ref_base,
                a_vcf,
                b_contig: b_bnd.contig,
                b_pos0: b_bnd.pos0,
                b_strands: b_bnd.strands,
                b_ref_base,
                b_vcf,
                inv_event_id,
            }),
            trid: None,
            vcf: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::svtype::SvType;
    use crate::core::variant::test_utils::make_temp_vcf;
    use rust_htslib::bcf::Read;

    fn parse_variant(
        record: &rust_htslib::bcf::Record,
        vcf_id: usize,
        contig: &str,
        args: &MergeArgsInner,
    ) -> Result<VariantInternal> {
        let format_cache = VariantInternal::build_header_format_cache(record.header())?;
        VariantInternal::from_vcf_record(record, vcf_id, contig, args, &None, &format_cache, None)
    }

    #[test]
    fn sawfish_bnd_alt_parses_mate_pos_to_end() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr5>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10001\tsawfish:0:0:0:0:0\tT\t]chr9:10362]T\t66\tPASS\tSVTYPE=BND;MATEID=sawfish:0:0:0:0:1\tGT:GQ:PL:AD\t1/1:6:100,6,0:0,2
chr1\t91917\tsawfish:0:10:0:0:0\tT\t[chr5:181462520[T\t101\tPASS\tSVTYPE=BND;MATEID=sawfish:0:10:0:0:1\tGT:GQ:PL:AD\t0/1:134:134,0,999:49,6
chr1\t139778\tsawfish:0:14:0:0:0\tT\tT[chr1:736702[\t999\tPASS\tSVTYPE=BND;MATEID=sawfish:0:14:0:0:1\tGT:GQ:PL:AD\t0/1:999:999,0,999:55,44
chr1\t717\tsawfish:0:577:0:0:0\tA\tA]chr1:149549563]\t999\tPASS\tSVTYPE=BND;MATEID=sawfish:0:577:0:0:1\tGT:GQ:PL:AD\t0/1:999:999,0,999:82,33
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let args = MergeArgsInner::default();

        let expected_mate_pos_0based = [10361_i64, 181_462_519, 736_701, 149_549_562];
        let expected_mate_contig = ["chr9", "chr5", "chr1", "chr1"];
        let expected_mate_id = [
            "0_sawfish:0:0:0:0:1",
            "0_sawfish:0:10:0:0:1",
            "0_sawfish:0:14:0:0:1",
            "0_sawfish:0:577:0:0:1",
        ];
        let expected_strands = [*b"--", *b"+-", *b"++", *b"-+"];
        for (idx, record) in reader.records().enumerate() {
            let record = record.unwrap();
            let v = parse_variant(&record, 0, "chr1", &args).unwrap();
            assert_eq!(v.svtype, SvType::BND);
            let bnd = v.bnd.as_ref().unwrap();
            assert_eq!(bnd.pos0, record.pos());
            assert_eq!(bnd.mate_pos, expected_mate_pos_0based[idx]);
            assert_eq!(bnd.mate_contig, expected_mate_contig[idx]);
            assert_eq!(bnd.mate_id, expected_mate_id[idx]);
            assert_eq!(bnd.strands, expected_strands[idx]);
        }
    }

    #[test]
    fn bnd_requires_mateid() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\tb1\tA\t]chr9:200]A\t.\tPASS\tSVTYPE=BND\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let err = parse_variant(&record, 0, "chr1", &args).unwrap_err();
        assert!(err.to_string().contains("MATEID missing in BND record"));
    }

    #[test]
    fn bnd_rejects_malformed_alt_even_with_chr2_and_end() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##INFO=<ID=CHR2,Number=1,Type=String,Description=\"\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">
##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\tb1\tA\tchr9:200\t.\tPASS\tSVTYPE=BND;CHR2=chr9;END=200;STRANDS=+-;MATEID=b2\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let err = parse_variant(&record, 0, "chr1", &args).unwrap_err();
        assert!(
            err.to_string().contains("BND ALT is missing '[' or ']'"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn bnd_pair_rejects_inconsistent_support_calls_between_mates() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr9>
##contig=<ID=chr1>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##INFO=<ID=SUPP_CALLS,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tbnd_a\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b;SUPP_CALLS=2\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
chr1\t200\tbnd_b\tA\t]chr9:100]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_a;SUPP_CALLS=3\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let mut records = reader.records();
        let a = records.next().unwrap().unwrap();
        let b = records.next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let bnd_a = parse_variant(&a, 0, "chr9", &args).unwrap();
        let bnd_b = parse_variant(&b, 0, "chr1", &args).unwrap();

        let err = VariantInternal::from_bnd_pair(bnd_a, bnd_b, &args).unwrap_err();
        assert!(
            err.to_string().contains("inconsistent SUPP_CALLS"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn bnd_pair_allows_homlen_explained_mate_pos_mismatch() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tbnd_a\tA\t]chr1:205]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b;HOMLEN=5\tGT\t0/1
chr1\t200\tbnd_b\tA\t]chr9:105]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_a;HOMLEN=5\tGT\t0/1
";

        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let args = MergeArgsInner::default();
        let mut it = reader.records();
        let r0 = it.next().unwrap().unwrap();
        let r1 = it.next().unwrap().unwrap();

        let a = parse_variant(&r0, 0, "chr9", &args).unwrap();
        let b = parse_variant(&r1, 0, "chr1", &args).unwrap();
        let event = VariantInternal::from_bnd_pair(a, b, &args).unwrap();

        let e = event.bnd_event.as_ref().unwrap();
        assert_eq!(e.a_contig, "chr1");
        assert_eq!(e.b_contig, "chr9");
        assert_eq!(e.a_pos0, 199);
        assert_eq!(e.b_pos0, 99);
    }

    #[test]
    fn bnd_coords_swap_when_contig_order_is_flipped() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tb1\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=b2\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let v = parse_variant(&record, 0, "chr9", &args).unwrap();
        assert_eq!(v.svtype, SvType::BND);
        assert_eq!(v.start as i64, 199);
        assert_eq!(v.end as i64, 99);
        assert_eq!(v.bnd_graph_id().unwrap(), "chr1_chr9_TRA_--");
    }

    #[test]
    fn bnd_coords_use_avg_fields_with_swap_semantics() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##INFO=<ID=AVG_START,Number=1,Type=String,Description=\"\">
##INFO=<ID=AVG_END,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tb1\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=b2;AVG_START=1000.0;AVG_END=2000.0\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let v = parse_variant(&record, 0, "chr9", &args).unwrap();
        assert_eq!(v.svtype, SvType::BND);
        assert_eq!(v.start, 1999.0);
        assert_eq!(v.end, 999.0);
    }

    #[test]
    fn bnd_strands_prefers_info_over_alt_inference() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tb1\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=b2;STRANDS=+-\tGT\t0/1
";
        let path = make_temp_vcf(vcf);
        let mut reader = bcf::Reader::from_path(&path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let args = MergeArgsInner::default();

        let v = parse_variant(&record, 0, "chr9", &args).unwrap();
        let bnd = v.bnd.as_ref().unwrap();
        assert_eq!(bnd.strands, *b"+-");
        assert_eq!(v.bnd_graph_id().unwrap(), "chr1_chr9_TRA_+-");
    }
}
