use super::{
    OutputSorter, SortConfig,
    keys::PackedSortKey,
    spill_payload::{
        decode_record_payload, decode_record_payload_into, duplicate_header_view,
        encode_record_payload, parse_payload_prefix,
    },
};
use crate::SpillSortError;
use rust_htslib::bcf::{self, Read};
use std::{
    fs,
    path::{Path, PathBuf},
};
use tempfile::NamedTempFile;

fn test_header() -> bcf::Header {
    let mut header = bcf::Header::new();
    header.push_record(br#"##fileformat=VCFv4.2"#);
    header.push_record(br#"##contig=<ID=chr1,length=1000>"#);
    header.push_record(br#"##contig=<ID=chr2,length=1000>"#);
    header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="">"#);
    header.push_sample(b"S1");
    header
}

fn record_factory() -> bcf::Writer {
    bcf::Writer::from_stdout(&test_header(), true, bcf::Format::Vcf)
        .expect("record factory should initialize")
}

fn make_record(
    factory: &bcf::Writer,
    rid: u32,
    pos: i64,
    id: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> bcf::Record {
    let mut record = factory.empty_record();
    record.set_rid(Some(rid));
    record.set_pos(pos);
    record
        .set_id(id.as_bytes())
        .expect("record id should be set for sorter test");
    record
        .set_alleles(&[ref_allele.as_bytes(), alt_allele.as_bytes()])
        .expect("alleles should be set for sorter test");
    record
        .push_genotypes(&[
            bcf::record::GenotypeAllele::Unphased(0),
            bcf::record::GenotypeAllele::Unphased(1),
        ])
        .expect("genotype should be set for sorter test");
    record
}

fn temp_output_path() -> (NamedTempFile, String) {
    let file = NamedTempFile::new().expect("output file should be creatable");
    let path = file.path().to_string_lossy().to_string();
    (file, path)
}

fn collect_sorted_records(path: &str) -> Vec<(PackedSortKey, String)> {
    let mut reader = bcf::Reader::from_path(path).expect("sorted output should be readable");
    reader
        .records()
        .map(|record| {
            let record = record.expect("record should parse");
            let key = PackedSortKey::from_record(&record).expect("sort key should be derived");
            let id = String::from_utf8(record.id().to_vec()).expect("record ID should be UTF-8");
            (key, id)
        })
        .collect()
}

fn finish_sorter_to_path(sorter: &mut OutputSorter, output_header: &bcf::Header, out_path: &str) {
    let mut writer = bcf::Writer::from_path(out_path, output_header, true, bcf::Format::Vcf)
        .expect("output writer should initialize");
    sorter
        .finish_with(|record| writer.write(record).map_err(SpillSortError::from))
        .expect("sort should finish");
    drop(writer);
}

fn run_tie_heavy_spill_sort_ids(
    blob_count: u64,
    max_open_files: usize,
    merge_fan_in: usize,
) -> Vec<String> {
    let factory = record_factory();
    let header_view = factory.header();
    let sorter_header = bcf::Header::from_template(header_view);
    let output_header = bcf::Header::from_template(header_view);
    let mut sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1, None, max_open_files, merge_fan_in)
            .expect("sort config should initialize for tie-heavy spill stress test"),
    )
    .expect("sorter should initialize");

    for blob_ordinal in 0..blob_count {
        sorter
            .push_blob_sorted_run(
                blob_ordinal,
                vec![make_record(
                    &factory,
                    0,
                    25,
                    &format!("tie_{blob_ordinal:03}"),
                    "A",
                    "C",
                )],
            )
            .expect("tie-heavy blob should be accepted");
    }

    let (_guard, out_path) = temp_output_path();
    finish_sorter_to_path(&mut sorter, &output_header, &out_path);

    collect_sorted_records(&out_path)
        .into_iter()
        .map(|(_, id)| id)
        .collect()
}

fn expected_tie_order_ids(blob_count: u64) -> Vec<String> {
    (0..blob_count)
        .map(|blob_ordinal| format!("tie_{blob_ordinal:03}"))
        .collect()
}

fn spill_file_paths(temp_dir: &Path) -> Vec<PathBuf> {
    let mut paths = fs::read_dir(temp_dir)
        .expect("sort temp directory should be readable")
        .map(|entry| {
            entry
                .expect("sort temp directory entry should be readable")
                .path()
        })
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "spill")
        })
        .collect::<Vec<_>>();
    paths.sort();
    paths
}

#[test]
fn parse_payload_prefix_rejects_invalid_magic() {
    let error = parse_payload_prefix(&[
        b'N', b'O', b'T', b'P', 0, 1, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 9,
    ])
    .expect_err("prefix parse should fail for non-HTSP payload");
    assert!(error.to_string().contains("magic mismatch"));
}

#[test]
fn encode_record_payload_uses_binary_magic_prefix() {
    let factory = record_factory();
    let record = make_record(&factory, 0, 42, "r1", "A", "C");

    let payload =
        encode_record_payload(&record).expect("record payload should encode successfully");
    assert!(
        payload.starts_with(b"HTSP"),
        "expected binary payload magic prefix, got first bytes: {:?}",
        &payload[..payload.len().min(8)]
    );
}

#[test]
fn decode_record_payload_into_roundtrips_binary_payload() {
    let factory = record_factory();
    let record = make_record(&factory, 0, 12, "r1", "A", "G");
    let payload =
        encode_record_payload(&record).expect("record payload should encode successfully");
    let decode_header = duplicate_header_view(&bcf::Header::from_template(factory.header()))
        .expect("duplicate header should initialize");
    let mut decoded = decode_header.empty_record();

    decode_record_payload_into(
        payload.as_slice(),
        record.rid().expect("record rid should be set"),
        record.pos(),
        &mut decoded,
    )
    .expect("binary payload should decode successfully");
    assert_eq!(
        decoded
            .to_vcf_string()
            .expect("decoded record should serialize as VCF"),
        record
            .to_vcf_string()
            .expect("source record should serialize as VCF")
    );
}

#[test]
fn decode_record_payload_into_rejects_truncated_binary_payload() {
    let factory = record_factory();
    let record = make_record(&factory, 0, 14, "r1", "A", "T");
    let mut payload =
        encode_record_payload(&record).expect("record payload should encode successfully");
    payload.pop();

    let decode_header = duplicate_header_view(&bcf::Header::from_template(factory.header()))
        .expect("duplicate header should initialize");
    let mut decoded = decode_header.empty_record();
    let error = decode_record_payload_into(
        payload.as_slice(),
        record.rid().expect("record rid should be set"),
        record.pos(),
        &mut decoded,
    )
    .expect_err("truncated binary payload should fail to decode");
    assert!(error.to_string().contains("truncated"));
}

#[test]
fn packed_record_payload_roundtrip_uses_binary_payload_bytes() {
    let factory = record_factory();
    let record = make_record(&factory, 0, 10, "packed_roundtrip", "A", "C");
    let payload = encode_record_payload(&record).expect("record payload encoding should succeed");
    assert!(
        payload.iter().any(|byte| *byte == 0),
        "packed payload should include binary bytes (at least one NUL byte)"
    );

    let decode_header_source = bcf::Header::from_template(factory.header());
    let decode_header = duplicate_header_view(&decode_header_source)
        .expect("decode header should duplicate for payload roundtrip");
    let decoded = decode_record_payload(payload.as_slice(), &decode_header)
        .expect("payload should decode back into a valid record");
    assert_eq!(
        decoded
            .to_vcf_string()
            .expect("decoded record should render to VCF"),
        record
            .to_vcf_string()
            .expect("original record should render to VCF")
    );
}

#[test]
fn packed_record_payload_decode_rejects_truncated_payload_with_specific_error() {
    let factory = record_factory();
    let record = make_record(&factory, 0, 11, "packed_truncated", "A", "C");
    let payload = encode_record_payload(&record).expect("record payload encoding should succeed");
    let truncation_len = payload.len() / 2;
    let truncated = payload
        .get(..truncation_len)
        .expect("truncated payload slice should exist");
    let decode_header_source = bcf::Header::from_template(factory.header());
    let decode_header = duplicate_header_view(&decode_header_source)
        .expect("decode header should duplicate for truncation test");
    let error = decode_record_payload(truncated, &decode_header)
        .expect_err("decode should fail for truncated packed payload");
    let message = error.to_string().to_ascii_lowercase();
    assert!(
        message.contains("truncat") || message.contains("length"),
        "expected truncation/length error, got: {message}"
    );
}

#[test]
fn packed_record_payload_decode_rejects_corrupt_payload_with_specific_error() {
    let factory = record_factory();
    let record = make_record(&factory, 0, 12, "packed_corrupt", "A", "C");
    let mut payload =
        encode_record_payload(&record).expect("record payload encoding should succeed");
    payload.fill(0);

    let decode_header_source = bcf::Header::from_template(factory.header());
    let decode_header = duplicate_header_view(&decode_header_source)
        .expect("decode header should duplicate for corruption test");
    let error = decode_record_payload(payload.as_slice(), &decode_header)
        .expect_err("decode should fail for corrupted packed payload");
    let message = error.to_string().to_ascii_lowercase();
    assert!(
        message.contains("corrupt") || message.contains("invalid") || message.contains("malformed"),
        "expected corruption/invalid error, got: {message}"
    );
}

#[test]
fn push_blob_interleaves_presorted_batches_before_spill() {
    let factory = record_factory();
    let header_view = factory.header();
    let sorter_header = bcf::Header::from_template(header_view);
    let output_header = bcf::Header::from_template(header_view);
    let mut sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1_000_000, None, 8, 4)
            .expect("sort config should initialize for sorter test"),
    )
    .expect("sorter should initialize");

    sorter
        .push_blob_sorted_run(
            1,
            vec![
                make_record(&factory, 0, 20, "b1", "A", "C"),
                make_record(&factory, 0, 40, "b2", "A", "C"),
            ],
        )
        .expect("blob should be accepted");
    sorter
        .push_blob_sorted_run(
            0,
            vec![
                make_record(&factory, 0, 10, "a1", "A", "C"),
                make_record(&factory, 0, 30, "a2", "A", "C"),
            ],
        )
        .expect("blob should be accepted");

    let (_guard, out_path) = temp_output_path();
    finish_sorter_to_path(&mut sorter, &output_header, &out_path);

    let records = collect_sorted_records(&out_path);
    let observed_ids = records.into_iter().map(|(_, id)| id).collect::<Vec<_>>();
    assert_eq!(observed_ids, vec!["a1", "b1", "a2", "b2"]);
}

#[test]
fn spills_and_kway_merges_without_loading_all_runs() {
    let factory = record_factory();
    let header_view = factory.header();
    let sorter_header = bcf::Header::from_template(header_view);
    let output_header = bcf::Header::from_template(header_view);
    let mut sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1, None, 4, 2).expect("sort config should initialize for spill test"),
    )
    .expect("sorter should initialize");

    for blob_ordinal in 0..6_u64 {
        sorter
            .push_blob_sorted_run(
                blob_ordinal,
                vec![make_record(
                    &factory,
                    0,
                    100 - blob_ordinal as i64,
                    &format!("id{blob_ordinal}"),
                    "A",
                    "C",
                )],
            )
            .expect("blob should be accepted");
    }

    let (_guard, out_path) = temp_output_path();
    finish_sorter_to_path(&mut sorter, &output_header, &out_path);

    let records = collect_sorted_records(&out_path);
    assert_eq!(records.len(), 6);
    assert!(
        records.windows(2).all(|pair| pair[0].0 <= pair[1].0),
        "records are not globally sorted"
    );
}

#[test]
fn spill_pending_runs_resets_pending_state_after_spill() {
    let factory = record_factory();
    let header_view = factory.header();
    let sorter_header = bcf::Header::from_template(header_view);
    let mut sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1_000_000, None, 8, 4)
            .expect("sort config should initialize for spill state reset test"),
    )
    .expect("sorter should initialize");

    sorter
        .push_blob_sorted_run(0, vec![make_record(&factory, 0, 10, "a0", "A", "C")])
        .expect("blob should be accepted");
    sorter
        .push_blob_sorted_run(1, vec![make_record(&factory, 0, 11, "a1", "A", "C")])
        .expect("blob should be accepted");

    assert!(
        !sorter.pending_runs.is_empty(),
        "pending runs should be buffered"
    );
    assert!(
        sorter.pending_bytes > 0,
        "pending bytes should track buffered runs"
    );
    assert!(
        sorter.spill_runs.is_empty(),
        "no spill runs should exist before explicit spill"
    );

    sorter
        .spill_pending_runs()
        .expect("spilling pending runs should succeed");

    assert!(
        sorter.pending_runs.is_empty(),
        "pending runs should be cleared after spill"
    );
    assert_eq!(
        sorter.pending_bytes, 0,
        "pending bytes should reset after spill"
    );
    assert_eq!(
        sorter.spill_runs.len(),
        1,
        "one spill run should be created"
    );
    assert!(
        sorter.spill_runs[0].path.exists(),
        "spill run should exist on disk after spill"
    );
}

#[test]
fn finish_with_cleans_up_spills_when_emit_record_fails() {
    let factory = record_factory();
    let sorter_header = bcf::Header::from_template(factory.header());
    let mut sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1, None, 4, 2)
            .expect("sort config should initialize for emit failure cleanup test"),
    )
    .expect("sorter should initialize");

    for blob_ordinal in 0..6_u64 {
        sorter
            .push_blob_sorted_run(
                blob_ordinal,
                vec![make_record(
                    &factory,
                    0,
                    10 + blob_ordinal as i64,
                    &format!("emit_fail_{blob_ordinal}"),
                    "A",
                    "C",
                )],
            )
            .expect("blob should be accepted");
    }

    assert!(
        !sorter.spill_runs.is_empty(),
        "setup should force spill runs before finish"
    );
    let spill_files_before = spill_file_paths(sorter.temp_dir.path());
    assert!(
        !spill_files_before.is_empty(),
        "spill files should exist before finish failure setup"
    );

    let error = sorter
        .finish_with(|_| Err(SpillSortError::message("injected emit failure".to_string())))
        .expect_err("finish should return emit failure");
    assert!(
        error.to_string().contains("injected emit failure"),
        "finish should preserve emit failure context"
    );

    let spill_files_after = spill_file_paths(sorter.temp_dir.path());
    assert!(
        spill_files_after.is_empty(),
        "spill files should be removed when finish emit fails; remaining files: {:?}",
        spill_files_after
    );
}

#[test]
fn spill_compaction_low_fan_in_matches_high_fan_in_tie_order() {
    let blob_count = 17;
    let expected_ids = expected_tie_order_ids(blob_count);
    let high_fan_in_ids = run_tie_heavy_spill_sort_ids(blob_count, 64, 64);
    let low_fan_in_ids = run_tie_heavy_spill_sort_ids(blob_count, 2, 2);

    assert_eq!(
        high_fan_in_ids, expected_ids,
        "high-fan-in merge changed deterministic tie ordering"
    );
    assert_eq!(
        low_fan_in_ids, expected_ids,
        "low-fan-in spill compaction changed deterministic tie ordering"
    );
}

#[test]
fn rejects_non_monotonic_blob_order() {
    let factory = record_factory();
    let sorter_header = bcf::Header::from_template(factory.header());
    let mut sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1_000_000, None, 8, 4)
            .expect("sort config should initialize for monotonicity test"),
    )
    .expect("sorter should initialize");

    let error = sorter
        .push_blob_sorted_run(
            0,
            vec![
                make_record(&factory, 0, 30, "x2", "A", "C"),
                make_record(&factory, 0, 10, "x1", "A", "C"),
            ],
        )
        .expect_err("non-monotonic blob should fail");
    assert!(error.to_string().contains("not sorted"));
}

#[test]
fn returns_error_when_temp_directory_cannot_be_created() {
    let factory = record_factory();
    let sorter_header = bcf::Header::from_template(factory.header());
    let temp_file = NamedTempFile::new().expect("temp file should be creatable");
    let error = OutputSorter::new(
        sorter_header,
        SortConfig::new(1, Some(temp_file.path().to_path_buf()), 4, 2)
            .expect("sort config should initialize for tmp-dir test"),
    )
    .expect_err("invalid tmp directory root should fail");
    assert!(error.to_string().contains("sort temp directory"));
}

#[test]
fn sorter_fails_when_sort_tmp_dir_is_unusable() {
    let factory = record_factory();
    let sorter_header = bcf::Header::from_template(factory.header());
    let temp_root = tempfile::TempDir::new().expect("temp root should be creatable");
    let file_path = temp_root.path().join("not_a_dir");
    std::fs::write(&file_path, b"x").expect("file should be writable");

    let error = OutputSorter::new(
        sorter_header,
        SortConfig::new(1, Some(file_path), 4, 2)
            .expect("sort config should initialize for tmp-dir usability test"),
    )
    .expect_err("non-directory tmp root should fail");
    assert!(error.to_string().contains("sort temp directory"));
}

#[test]
fn sorter_temp_dir_uses_spill_sort_prefix() {
    let factory = record_factory();
    let sorter_header = bcf::Header::from_template(factory.header());
    let sorter = OutputSorter::new(
        sorter_header,
        SortConfig::new(1, None, 4, 2).expect("sort config should initialize"),
    )
    .expect("sorter should initialize");
    let dir_name = sorter
        .temp_dir
        .path()
        .file_name()
        .expect("temp dir should have file name")
        .to_string_lossy()
        .to_string();
    assert!(
        dir_name.starts_with("spill-sort-"),
        "expected spill-sort temp dir prefix, got {dir_name}"
    );
}
