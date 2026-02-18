use super::VariantInternal;
use crate::core::svtype::SvType;
use crate::utils::util::Result;
use std::{
    fs,
    sync::atomic::{AtomicU64, Ordering},
    time::SystemTime,
};

static TEMP_VCF_COUNTER: AtomicU64 = AtomicU64::new(0);

pub fn make_temp_vcf(contents: &str) -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    let nanos = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let counter = TEMP_VCF_COUNTER.fetch_add(1, Ordering::Relaxed);
    path.push(format!("svx_test_variant_{nanos}_{counter}.vcf"));
    fs::write(&path, contents).unwrap();
    path
}

pub fn from_parts(
    vcf_id: usize,
    id: String,
    svtype: SvType,
    start: f64,
    end: f64,
) -> Result<VariantInternal> {
    let interval = match svtype {
        SvType::DELETION | SvType::INVERSION | SvType::DUPLICATION | SvType::CNV => {
            Some([start, end])
        }
        _ => None,
    };

    let id_list = vec![id.clone()];

    Ok(VariantInternal {
        start,
        end,
        interval,
        svlen: end,
        vcf_id,
        support_mask: {
            let word = vcf_id / (u64::BITS as usize);
            let bit = vcf_id % (u64::BITS as usize);
            let mut mask = vec![0u64; word + 1];
            mask[word] |= 1u64 << bit;
            mask
        },
        support_calls: 1,
        sample_id: vcf_id,
        start_mean: start + 1.0,
        start_variance: 0.0,
        svlen_mean: end,
        svlen_variance: 0.0,
        svtype,
        id,
        id_list,
        svclaim: None,
        svclaims: Vec::new(),
        bnd: None,
        bnd_event: None,
        index: 0,
        max_dist: 1000.0,
        info_hash: 0,
        trid: None,
        sequence: None,
        vcf: None,
    })
}

pub fn insertion_fixture_variants(svtype: SvType) -> Vec<VariantInternal> {
    vec![
        from_parts(0, "var1".to_string(), svtype, 10.0, 5.0).expect("var1 should build"),
        from_parts(0, "var2".to_string(), svtype, 1.0, 5.0).expect("var2 should build"),
        from_parts(0, "var3".to_string(), svtype, 18.0, 5.0).expect("var3 should build"),
        from_parts(1, "var4".to_string(), svtype, 12.0, 7.0).expect("var4 should build"),
        from_parts(1, "var5".to_string(), svtype, 10.0, 5.0).expect("var5 should build"),
        from_parts(1, "var6".to_string(), svtype, 30.0, 30.0).expect("var6 should build"),
        from_parts(1, "var7".to_string(), svtype, 0.0, 0.0).expect("var7 should build"),
        from_parts(2, "var8".to_string(), svtype, 12.0, 12.0).expect("var8 should build"),
        from_parts(2, "var9".to_string(), svtype, 15.0, 15.0).expect("var9 should build"),
        from_parts(2, "var10".to_string(), svtype, 20.0, 20.0).expect("var10 should build"),
        from_parts(2, "var11".to_string(), svtype, 28.0, 28.0).expect("var11 should build"),
        from_parts(3, "var12".to_string(), svtype, 25.0, 25.0).expect("var12 should build"),
        from_parts(4, "var13".to_string(), svtype, 22.0, 22.0).expect("var13 should build"),
    ]
}

pub fn insertion_fixture_variants_disjoint(svtype: SvType) -> Vec<VariantInternal> {
    vec![
        from_parts(0, "var1".to_string(), svtype, 10.0, 5.0).expect("var1 should build"),
        from_parts(1, "var2".to_string(), svtype, 1.0, 5.0).expect("var2 should build"),
        from_parts(2, "var3".to_string(), svtype, 18.0, 5.0).expect("var3 should build"),
        from_parts(3, "var4".to_string(), svtype, 12.0, 7.0).expect("var4 should build"),
        from_parts(4, "var5".to_string(), svtype, 10.0, 5.0).expect("var5 should build"),
        from_parts(5, "var6".to_string(), svtype, 30.0, 30.0).expect("var6 should build"),
        from_parts(6, "var7".to_string(), svtype, 0.0, 0.0).expect("var7 should build"),
        from_parts(7, "var8".to_string(), svtype, 12.0, 12.0).expect("var8 should build"),
        from_parts(8, "var9".to_string(), svtype, 15.0, 15.0).expect("var9 should build"),
        from_parts(9, "var10".to_string(), svtype, 20.0, 20.0).expect("var10 should build"),
        from_parts(10, "var11".to_string(), svtype, 28.0, 28.0).expect("var11 should build"),
        from_parts(11, "var12".to_string(), svtype, 25.0, 25.0).expect("var12 should build"),
        from_parts(12, "var13".to_string(), svtype, 22.0, 22.0).expect("var13 should build"),
    ]
}
