use super::{
    bnd::{BndMateSelectionScope, bucket_bnd_variants, pair_bnd_breakends},
    merge,
    progress::{
        PipelineQueueSnapshot, QueueDepthSnapshot, QueueDepthTracker, build_merge_progress_message,
        compute_pipeline_queue_capacities, compute_progress_enabled,
        resolve_pipeline_queue_capacities,
    },
    shutdown::finalize_merge_threads,
};
use crate::cli::{Cli, Command, MergeArgs, MergeArgsInner};
use crate::core::variant::VariantInternal;
use crate::io::positions_reader::PositionTuple;
use crate::utils::util::init_logger;
use clap::Parser;
use log::LevelFilter;
use rust_htslib::bcf::{self, Header, Read, record::GenotypeAllele};
use std::{collections::HashSet, fs, path::PathBuf, thread};

fn make_temp_path(stem: &str, ext: &str) -> PathBuf {
    let prefix = format!("svx_test_{stem}_");
    let suffix = format!(".{ext}");
    let (_, path) = tempfile::Builder::new()
        .prefix(&prefix)
        .suffix(&suffix)
        .tempfile_in(std::env::temp_dir())
        .expect("temp file should be creatable")
        .keep()
        .expect("temp file should be persistable");
    path
}

fn make_temp_vcf(contents: &str) -> PathBuf {
    let path = make_temp_path("merge_bucket_bnd", "vcf");
    fs::write(&path, contents).expect("test VCF should be writable");
    path
}

fn parse_merge_args(args: &[&str]) -> MergeArgs {
    let parsed = Cli::try_parse_from(args).expect("CLI parse should succeed");
    let Command::Merge(args) = parsed.command;
    args
}

fn parse_single_variant(vcf: &str, contig: &str, args: &MergeArgsInner) -> VariantInternal {
    let path = make_temp_vcf(vcf);
    let mut reader =
        rust_htslib::bcf::Reader::from_path(&path).expect("reader should open test VCF");
    let record = reader
        .records()
        .next()
        .expect("record should exist")
        .expect("record should parse");
    VariantInternal::from_vcf_record(&record, 0, contig, args, &None).expect("record should parse")
}

fn parse_variants_with_vcf_id(
    vcf: &str,
    vcf_id: usize,
    args: &MergeArgsInner,
) -> Vec<VariantInternal> {
    let path = make_temp_vcf(vcf);
    let mut reader =
        rust_htslib::bcf::Reader::from_path(&path).expect("reader should open test VCF");
    let records = reader
        .records()
        .map(|record| record.expect("record should parse"))
        .collect::<Vec<_>>();
    records
        .iter()
        .map(|record| {
            let contig = std::str::from_utf8(
                reader
                    .header()
                    .rid2name(record.rid().expect("record rid should exist"))
                    .expect("rid should map to contig"),
            )
            .expect("contig name should be UTF-8");
            VariantInternal::from_vcf_record(record, vcf_id, contig, args, &None)
                .expect("breakend should parse")
        })
        .collect::<Vec<_>>()
}

fn make_sample_bnd_breakends(
    vcf_id: usize,
    breakend_a_pos_1based: i64,
    breakend_b_pos_1based: i64,
    args: &MergeArgsInner,
) -> Vec<VariantInternal> {
    let raw_a_id = format!("sample{vcf_id}_a");
    let raw_b_id = format!("sample{vcf_id}_b");
    let vcf = format!(
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t{breakend_a_pos_1based}\t{raw_a_id}\tA\t]chr2:{breakend_b_pos_1based}]A\t.\tPASS\tSVTYPE=BND;MATEID={raw_b_id}\tGT\t0/1
chr2\t{breakend_b_pos_1based}\t{raw_b_id}\tA\t]chr1:{breakend_a_pos_1based}]A\t.\tPASS\tSVTYPE=BND;MATEID={raw_a_id}\tGT\t0/1
"
    );
    parse_variants_with_vcf_id(&vcf, vcf_id, args)
}

fn make_temp_indexed_vcf_for_svtype_filter() -> PathBuf {
    let path = make_temp_path("merge_svtype_filter", "vcf");

    let mut header = Header::new();
    header.push_record(br#"##fileformat=VCFv4.2"#);
    header.push_record(br#"##contig=<ID=chr1,length=1000000>"#);
    header.push_record(br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">"#);
    header.push_record(br#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">"#);
    header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    header.push_record(br#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#);
    header.push_record(
        br#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype Likelihoods">"#,
    );
    header.push_record(br#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">"#);
    header.push_sample(b"sample1");

    let mut writer = bcf::Writer::from_path(&path, &header, false, bcf::Format::Vcf)
        .expect("test writer should initialize");

    let mut ins = writer.empty_record();
    ins.set_rid(Some(0));
    ins.set_pos(100);
    ins.set_id(b"ins1").expect("INS id should be set");
    ins.set_alleles(&[b"N", b"<INS>"])
        .expect("INS alleles should be set");
    ins.push_info_string(b"SVTYPE", &[b"INS"])
        .expect("INS SVTYPE should be set");
    ins.push_info_integer(b"SVLEN", &[10])
        .expect("INS SVLEN should be set");
    ins.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
        .expect("INS genotype should be set");
    ins.push_format_integer(b"GQ", &[0])
        .expect("INS GQ should be set");
    ins.push_format_integer(b"PL", &[0, 0, 0])
        .expect("INS PL should be set");
    ins.push_format_integer(b"AD", &[0, 0])
        .expect("INS AD should be set");
    writer.write(&ins).expect("INS record should be written");

    // DEL intentionally lacks SVLEN so it is invalid if parsed.
    let mut del = writer.empty_record();
    del.set_rid(Some(0));
    del.set_pos(200);
    del.set_id(b"del1").expect("DEL id should be set");
    del.set_alleles(&[b"N", b"<DEL>"])
        .expect("DEL alleles should be set");
    del.push_info_string(b"SVTYPE", &[b"DEL"])
        .expect("DEL SVTYPE should be set");
    del.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
        .expect("DEL genotype should be set");
    del.push_format_integer(b"GQ", &[0])
        .expect("DEL GQ should be set");
    del.push_format_integer(b"PL", &[0, 0, 0])
        .expect("DEL PL should be set");
    del.push_format_integer(b"AD", &[0, 0])
        .expect("DEL AD should be set");
    writer.write(&del).expect("DEL record should be written");

    drop(writer);
    bcf::index::build(&path, None, 1, bcf::index::Type::Tbx)
        .expect("tabix index should be created");
    path
}

fn make_temp_indexed_vcf_for_cnv_svtype_filter() -> PathBuf {
    let path = make_temp_path("merge_cnv_svtype_filter", "vcf");

    let mut header = Header::new();
    header.push_record(br#"##fileformat=VCFv4.2"#);
    header.push_record(br#"##contig=<ID=chr1,length=1000000>"#);
    header.push_record(br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">"#);
    header.push_record(br#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">"#);
    header.push_record(br#"##INFO=<ID=SVCLAIM,Number=1,Type=String,Description="">"#);
    header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    header.push_record(br#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#);
    header.push_record(
        br#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype Likelihoods">"#,
    );
    header.push_record(br#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">"#);
    header.push_record(br#"##FORMAT=<ID=CN,Number=1,Type=Float,Description="">"#);
    header.push_record(br#"##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="">"#);
    header.push_sample(b"sample1");

    let mut writer = bcf::Writer::from_path(&path, &header, false, bcf::Format::Vcf)
        .expect("test writer should initialize");

    let mut ins = writer.empty_record();
    ins.set_rid(Some(0));
    ins.set_pos(100);
    ins.set_id(b"ins1").expect("INS id should be set");
    ins.set_alleles(&[b"N", b"<INS>"])
        .expect("INS alleles should be set");
    ins.push_info_string(b"SVTYPE", &[b"INS"])
        .expect("INS SVTYPE should be set");
    ins.push_info_integer(b"SVLEN", &[10])
        .expect("INS SVLEN should be set");
    ins.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
        .expect("INS genotype should be set");
    ins.push_format_integer(b"GQ", &[0])
        .expect("INS GQ should be set");
    ins.push_format_integer(b"PL", &[0, 0, 0])
        .expect("INS PL should be set");
    ins.push_format_integer(b"AD", &[0, 0])
        .expect("INS AD should be set");
    writer.write(&ins).expect("INS record should be written");

    let mut cnv = writer.empty_record();
    cnv.set_rid(Some(0));
    cnv.set_pos(200);
    cnv.set_id(b"sawfish:CNV:0:0:200")
        .expect("CNV id should be set");
    cnv.set_alleles(&[b"N", b"<DEL>"])
        .expect("CNV alleles should be set");
    cnv.push_info_string(b"SVTYPE", &[b"DEL"])
        .expect("CNV SVTYPE should be set");
    cnv.push_info_integer(b"SVLEN", &[1000])
        .expect("CNV SVLEN should be set");
    cnv.push_info_string(b"SVCLAIM", &[b"D"])
        .expect("CNV SVCLAIM should be set");
    cnv.push_genotypes(&[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)])
        .expect("CNV genotype should be set");
    cnv.push_format_float(b"CN", &[0.0])
        .expect("CNV CN should be set");
    cnv.push_format_float(b"CNQ", &[50.0])
        .expect("CNV CNQ should be set");
    writer.write(&cnv).expect("CNV record should be written");

    drop(writer);
    bcf::index::build(&path, None, 1, bcf::index::Type::Tbx)
        .expect("tabix index should be created");
    path
}

#[test]
fn merge_with_svtype_filter_skips_unselected_invalid_records() {
    init_logger();

    let in_vcf = make_temp_indexed_vcf_for_svtype_filter();
    let out_vcf = make_temp_path("merge_svtype_filter_out", "vcf");

    let args = parse_merge_args(&[
        "svx",
        "merge",
        "--vcf",
        in_vcf.to_str().expect("input path should be UTF-8"),
        "--force-single",
        "--svtype",
        "INS",
        "--output",
        out_vcf.to_str().expect("output path should be UTF-8"),
    ]);

    let result = merge(args);
    assert!(
        result.is_ok(),
        "Expected merge to succeed because invalid DEL record is unselected by --svtype INS, got: {:?}",
        result.err()
    );
}

#[test]
fn merge_with_svtype_filter_processes_cnv_records() {
    init_logger();

    let in_vcf = make_temp_indexed_vcf_for_cnv_svtype_filter();
    let out_vcf = make_temp_path("merge_svtype_cnv_filter_out", "vcf");

    let args = parse_merge_args(&[
        "svx",
        "merge",
        "--vcf",
        in_vcf.to_str().expect("input path should be UTF-8"),
        "--force-single",
        "--svtype",
        "CNV",
        "--output",
        out_vcf.to_str().expect("output path should be UTF-8"),
    ]);

    let result = merge(args);
    assert!(
        result.is_ok(),
        "Expected merge to succeed for --svtype CNV, got: {:?}",
        result.err()
    );

    let out = fs::read_to_string(&out_vcf).expect("output VCF should be readable");
    let non_header_lines: Vec<&str> = out.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        non_header_lines.len(),
        1,
        "Expected exactly one CNV record in output for --svtype CNV"
    );
    assert!(
        out.contains("sawfish:CNV:"),
        "Expected CNV record ID to be present in output"
    );
    let format_col = non_header_lines[0]
        .split('\t')
        .nth(8)
        .expect("expected FORMAT column");
    assert_eq!(
        format_col, "GT:CN:CNQ",
        "Expected CNV output to use GT:CN:CNQ format"
    );
    let info_col = non_header_lines[0]
        .split('\t')
        .nth(7)
        .expect("expected INFO column");
    assert!(
        info_col.contains("SVCLAIM=D"),
        "Expected CNV output INFO to retain SVCLAIM for unanimous merged CNV inputs, got INFO={info_col}"
    );
}

#[test]
fn merge_with_svtype_filter_all_and_cnv_processes_both_classes() {
    init_logger();

    let in_vcf = make_temp_indexed_vcf_for_cnv_svtype_filter();
    let out_vcf = make_temp_path("merge_svtype_all_cnv_filter_out", "vcf");

    let args = parse_merge_args(&[
        "svx",
        "merge",
        "--vcf",
        in_vcf.to_str().expect("input path should be UTF-8"),
        "--force-single",
        "--svtype",
        "ALL,CNV",
        "--output",
        out_vcf.to_str().expect("output path should be UTF-8"),
    ]);

    let result = merge(args);
    assert!(
        result.is_ok(),
        "Expected merge to succeed for --svtype ALL,CNV, got: {:?}",
        result.err()
    );

    let out = fs::read_to_string(&out_vcf).expect("output VCF should be readable");
    let non_header_lines: Vec<&str> = out.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        non_header_lines.len(),
        2,
        "Expected INS and CNV records in output for --svtype ALL,CNV"
    );
    assert!(out.contains("ins1"));
    assert!(out.contains("sawfish:CNV:"));
    let mut ins_format: Option<&str> = None;
    let mut cnv_format: Option<&str> = None;
    for line in non_header_lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols[2].contains("ins1") {
            ins_format = cols.get(8).copied();
        }
        if cols[2].contains("sawfish:CNV:") {
            cnv_format = cols.get(8).copied();
        }
    }
    assert_eq!(ins_format, Some("GT:GQ:PL:AD"));
    assert_eq!(cnv_format, Some("GT:CN:CNQ"));
}

#[test]
fn merge_completes_with_tiny_pipeline_queues() {
    init_logger();

    let in_vcf = make_temp_indexed_vcf_for_svtype_filter();
    let out_vcf = make_temp_path("merge_tiny_pipeline_queue", "vcf");

    let args = parse_merge_args(&[
        "svx",
        "merge",
        "--vcf",
        in_vcf.to_str().expect("input path should be UTF-8"),
        "--force-single",
        "--svtype",
        "INS",
        "--blob-queue-capacity",
        "1",
        "--result-queue-capacity",
        "1",
        "--output",
        out_vcf.to_str().expect("output path should be UTF-8"),
    ]);

    let result = merge(args);
    assert!(
        result.is_ok(),
        "Expected merge to complete with tiny pipeline queues, got: {:?}",
        result.err()
    );

    let out = fs::read_to_string(&out_vcf).expect("output VCF should be readable");
    let non_header_lines: Vec<&str> = out.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        non_header_lines.len(),
        1,
        "Expected one INS record when merging with tiny pipeline queues"
    );
}

#[test]
fn bucket_bnd_variants_groups_across_contigs_by_tra_graph_id() {
    init_logger();

    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tbnd_a\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
chr1\t200\tbnd_b\tA\t]chr9:100]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
";
    let path = make_temp_vcf(vcf);
    let mut reader =
        rust_htslib::bcf::Reader::from_path(&path).expect("reader should open test VCF");
    let args = MergeArgsInner::default();

    let mut it = reader.records();
    let r0 = it
        .next()
        .expect("first record should exist")
        .expect("first record should parse");
    let r1 = it
        .next()
        .expect("second record should exist")
        .expect("second record should parse");

    let v0 =
        VariantInternal::from_vcf_record(&r0, 0, "chr9", &args, &None).expect("v0 should parse");
    let v1 =
        VariantInternal::from_vcf_record(&r1, 0, "chr1", &args, &None).expect("v1 should parse");

    let events = pair_bnd_breakends(vec![v0, v1], &args).expect("BND pairing should succeed");
    let buckets = bucket_bnd_variants(events).expect("BND bucketing should succeed");
    assert_eq!(buckets.len(), 1);
    let (k, v) = buckets.into_iter().next().expect("one BND bucket expected");
    assert_eq!(k, "chr1_chr9_TRA_--_--");
    assert_eq!(v.len(), 1);
}

#[test]
fn bucket_bnd_variants_splits_by_orientation_key() {
    init_logger();

    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\tx1_a\tA\tA[chr9:200[\t.\tPASS\tSVTYPE=BND;MATEID=x1_b\tGT\t0/1
chr9\t200\tx1_b\tA\t]chr1:100]A\t.\tPASS\tSVTYPE=BND;MATEID=x1_a\tGT\t0/1
chr1\t100\tx2_a\tA\t]chr9:200]A\t.\tPASS\tSVTYPE=BND;MATEID=x2_b\tGT\t0/1
chr9\t200\tx2_b\tA\tA[chr1:100[\t.\tPASS\tSVTYPE=BND;MATEID=x2_a\tGT\t0/1
";

    let path = make_temp_vcf(vcf);
    let mut reader =
        rust_htslib::bcf::Reader::from_path(&path).expect("reader should open test VCF");
    let args = MergeArgsInner::default();

    let records = reader
        .records()
        .map(|record| record.expect("record should parse"))
        .collect::<Vec<_>>();
    let breakends = records
        .iter()
        .map(|record| {
            let contig = std::str::from_utf8(
                reader
                    .header()
                    .rid2name(record.rid().expect("record rid should exist"))
                    .expect("rid should map to contig"),
            )
            .expect("contig name should be UTF-8");
            VariantInternal::from_vcf_record(record, 0, contig, &args, &None)
                .expect("breakend should parse")
        })
        .collect::<Vec<_>>();

    let events = pair_bnd_breakends(breakends, &args).expect("BND pairing should succeed");
    assert_eq!(events.len(), 2);

    let buckets = bucket_bnd_variants(events).expect("BND bucketing should succeed");
    assert_eq!(
        buckets.len(),
        2,
        "Expected different orientation translocations to be bucketed separately"
    );
}

#[test]
fn bnd_pairing_orders_events_stably_by_vcf_id() {
    init_logger();

    let args = MergeArgsInner::default();
    let mut breakends = Vec::new();
    for vcf_id in (0..10).rev() {
        let breakend_a_pos_1based = 100 + vcf_id as i64 * 10;
        let breakend_b_pos_1based = 200 + vcf_id as i64 * 10;
        breakends.extend(make_sample_bnd_breakends(
            vcf_id,
            breakend_a_pos_1based,
            breakend_b_pos_1based,
            &args,
        ));
    }

    let mut observed_orders: HashSet<Vec<usize>> = HashSet::new();
    for _ in 0..64 {
        let events =
            pair_bnd_breakends(breakends.clone(), &args).expect("BND pairing should succeed");
        observed_orders.insert(events.iter().map(|event| event.vcf_id).collect::<Vec<_>>());
    }

    assert_eq!(
        observed_orders.len(),
        1,
        "BND event ordering should be stable across repeated runs, observed orders: {observed_orders:?}"
    );
    assert_eq!(
        observed_orders
            .into_iter()
            .next()
            .expect("expected one observed ordering"),
        (0..10).collect::<Vec<_>>()
    );
}

#[test]
fn bnd_pairing_warns_and_skips_if_mate_missing_even_if_alt_gt() {
    init_logger();

    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr9>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr9\t100\tbnd_a\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT:GQ:PL:AD\t0/1:0:0,0,0:0,0
";
    let args = MergeArgsInner::default();
    let v = parse_single_variant(vcf, "chr9", &args);

    let events = pair_bnd_breakends(vec![v], &args).expect("BND pairing should succeed");
    assert_eq!(events.len(), 0);
}

#[test]
fn bnd_pairing_skips_self_pointing_orphan_if_nocall() {
    init_logger();

    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr4>
##FILTER=<ID=MinQUAL,Description=\"\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr4\t92649032\tbnd_a\tA\t[chr4:92649032[A\t0\tMinQUAL\tSVTYPE=BND;MATEID=bnd_b\tGT:GQ:PL:AD\t./.:.:0,0,0:0,0
";
    let args = MergeArgsInner::default();
    let v = parse_single_variant(vcf, "chr4", &args);

    let events = pair_bnd_breakends(vec![v], &args).expect("BND pairing should succeed");
    assert_eq!(events.len(), 0);
}

#[test]
fn bnd_pairing_skips_orphan_if_no_alt_gt() {
    init_logger();

    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chrEBV>
##FILTER=<ID=MinQUAL,Description=\"\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chrEBV\t39214\tbnd_a\tG\t[chrEBV:39224[G\t0\tMinQUAL\tSVTYPE=BND;MATEID=bnd_b\tGT:GQ:PL:AD\t0/0:999:0,999,999:693,0
";
    let args = MergeArgsInner::default();
    let v = parse_single_variant(vcf, "chrEBV", &args);

    let events = pair_bnd_breakends(vec![v], &args).expect("BND pairing should succeed");
    assert_eq!(events.len(), 0);
}

#[test]
fn progress_auto_mode_requires_tty_and_info_logging() {
    assert!(compute_progress_enabled(
        false,
        false,
        true,
        LevelFilter::Info
    ));
    assert!(!compute_progress_enabled(
        false,
        false,
        false,
        LevelFilter::Info
    ));
    assert!(!compute_progress_enabled(
        false,
        false,
        true,
        LevelFilter::Debug
    ));
    assert!(!compute_progress_enabled(
        false,
        false,
        true,
        LevelFilter::Trace
    ));
}

#[test]
fn progress_flag_precedence_and_constraints() {
    assert!(!compute_progress_enabled(
        true,
        true,
        true,
        LevelFilter::Info
    ));
    assert!(compute_progress_enabled(
        true,
        false,
        true,
        LevelFilter::Info
    ));
    assert!(!compute_progress_enabled(
        true,
        false,
        false,
        LevelFilter::Info
    ));
    assert!(!compute_progress_enabled(
        true,
        false,
        true,
        LevelFilter::Debug
    ));
}

#[test]
fn pipeline_queue_capacities_follow_expected_bounds() {
    assert_eq!(compute_pipeline_queue_capacities(0), (4, 8));
    assert_eq!(compute_pipeline_queue_capacities(1), (4, 8));
    assert_eq!(compute_pipeline_queue_capacities(4), (8, 16));
    assert_eq!(compute_pipeline_queue_capacities(32), (64, 128));
}

#[test]
fn pipeline_queue_capacity_overrides_take_precedence() {
    assert_eq!(
        resolve_pipeline_queue_capacities(4, Some(1), Some(2)),
        (1, 2)
    );
    assert_eq!(resolve_pipeline_queue_capacities(4, Some(3), None), (3, 16));
    assert_eq!(resolve_pipeline_queue_capacities(4, None, Some(5)), (8, 5));
}

#[test]
fn queue_depth_tracker_tracks_current_and_peak() {
    let tracker = QueueDepthTracker::default();

    assert_eq!(tracker.snapshot().current, 0);
    assert_eq!(tracker.snapshot().peak, 0);

    tracker.increment();
    tracker.increment();
    assert_eq!(tracker.snapshot().current, 2);
    assert_eq!(tracker.snapshot().peak, 2);

    tracker.decrement();
    assert_eq!(tracker.snapshot().current, 1);
    assert_eq!(tracker.snapshot().peak, 2);

    tracker.decrement();
    tracker.decrement();
    assert_eq!(tracker.snapshot().current, 0);
    assert_eq!(tracker.snapshot().peak, 2);
}

#[test]
fn merge_progress_message_includes_queue_depths() {
    let no_queue = build_merge_progress_message(10, 8, 6, 4, None);
    assert_eq!(no_queue, "read=10 merged=8 written=6 records=4");

    let queue_snapshot = PipelineQueueSnapshot {
        blob: QueueDepthSnapshot {
            current: 3,
            peak: 9,
        },
        result: QueueDepthSnapshot {
            current: 2,
            peak: 7,
        },
    };
    let with_queue = build_merge_progress_message(10, 8, 6, 4, Some(queue_snapshot));
    assert_eq!(
        with_queue,
        "read=10 merged=8 written=6 records=4 blob_q=3/9 result_q=2/7"
    );
}

#[test]
fn finalize_merge_threads_aggregates_failures_from_all_threads() {
    let reader_thread = thread::spawn(|| -> crate::utils::util::Result<()> {
        Err(crate::svx_error!("reader failed"))
    });
    let writer_thread = thread::spawn(|| -> crate::utils::util::Result<()> {
        Err(crate::svx_error!("writer failed"))
    });
    let progress_thread = Some(thread::spawn(|| {
        panic!("progress failed");
    }));

    let result = finalize_merge_threads(reader_thread, writer_thread, None, progress_thread);
    let error = result.expect_err("expected merged shutdown errors");
    let message = error.to_string();
    assert!(message.contains("Multiple thread shutdown errors"));
    assert!(message.contains("Reader thread failed: reader failed"));
    assert!(message.contains("Writer thread failed: writer failed"));
    assert!(message.contains("Progress thread panicked: progress failed"));
}

#[test]
fn bnd_mate_selection_scope_excludes_unselected_contig() {
    let selected_contigs = vec!["chr9".to_string()];
    let scope = BndMateSelectionScope::new(&selected_contigs, None);
    assert!(!scope.mate_is_in_scope("chr1", 200));
}

#[test]
fn bnd_mate_selection_scope_excludes_mate_outside_selected_range() {
    let selected_contigs = vec!["chr1".to_string()];
    let target_positions = vec![PositionTuple {
        contig: "chr1".to_string(),
        start: 99,
        end: Some(199),
    }];
    let scope = BndMateSelectionScope::new(&selected_contigs, Some(&target_positions));
    assert!(!scope.mate_is_in_scope("chr1", 200));
}

#[test]
fn bnd_mate_selection_scope_includes_mate_inside_selected_range() {
    let selected_contigs = vec!["chr1".to_string()];
    let target_positions = vec![PositionTuple {
        contig: "chr1".to_string(),
        start: 99,
        end: Some(199),
    }];
    let scope = BndMateSelectionScope::new(&selected_contigs, Some(&target_positions));
    assert!(scope.mate_is_in_scope("chr1", 199));
}

#[test]
fn variant_blob_is_available_via_types_submodule_path() {
    let _ = std::any::TypeId::of::<super::types::VariantBlob>();
}
