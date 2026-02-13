use crate::{io::tpool::HtsThreadPool, utils::util::Result};
use rust_htslib::bcf::{self, Header, HeaderRecord, Read};
use std::{
    collections::{HashMap, HashSet},
    ffi::OsString,
    fs::File,
    io::Read as ReadIo,
    path::{Path, PathBuf},
    sync::Arc,
};

#[derive(Debug, Clone)]
pub struct SampleMapping {
    // Maps (vcf_index, original_sample_pos) -> merged_vcf_pos
    pub index_map: HashMap<(usize, usize), usize>,
    // reverse mapping (debug)
    pub reverse_map: HashMap<usize, (usize, usize)>,
}

impl SampleMapping {
    pub fn new() -> Self {
        SampleMapping {
            index_map: HashMap::new(),
            reverse_map: HashMap::new(),
        }
    }
}

impl Default for SampleMapping {
    fn default() -> Self {
        Self::new()
    }
}

pub struct VcfReader {
    pub reader: bcf::IndexedReader,
    pub current_record: bcf::Record,
    pub sample_n: u32,
    pub index: usize,
    _io_tpool: Option<Arc<HtsThreadPool>>,
}

fn add_extension(path: &Path, ext: &str) -> PathBuf {
    let mut out = path.to_path_buf();
    let new_ext: OsString = match path.extension() {
        Some(old) => {
            let mut s = old.to_os_string();
            s.push(".");
            s.push(ext);
            s
        }
        None => OsString::from(ext),
    };
    out.set_extension(new_ext);
    out
}

fn is_indexed(file: &Path) -> bool {
    add_extension(file, "csi").exists() || add_extension(file, "tbi").exists()
}

fn has_gzip_magic(path: &Path) -> Result<bool> {
    let mut f = File::open(path)
        .map_err(|e| crate::svx_error!("Failed to open {}: {e}", path.display()))?;
    let mut m = [0u8; 2];
    let n = std::io::Read::read(&mut f, &mut m)
        .map_err(|e| crate::svx_error!("Failed to read {}: {e}", path.display()))?;
    Ok(n == 2 && m == [0x1f, 0x8b])
}

fn looks_like_bcf(p: &[u8]) -> bool {
    p.len() >= 5 && &p[..3] == b"BCF" && p[3] == 0x02
}

fn looks_like_vcf(p: &[u8]) -> bool {
    let mut p = p;
    // Strip the UTF-8 BOM if present
    if p.starts_with(&[0xEF, 0xBB, 0xBF]) {
        p = &p[3..];
    }
    if p.starts_with(b"##fileformat=VCF") {
        return true;
    }
    false
}

pub fn validate_bgzip_vcf(file: &Path) -> Result<()> {
    let is_bgzf = rust_htslib::bgzf::is_bgzip(file).map_err(|e| {
        crate::svx_error!(
            "Failed to determine whether {} is BGZF-compressed: {e}",
            file.display()
        )
    })?;

    if !is_bgzf {
        if has_gzip_magic(file)? {
            return Err(crate::svx_error!(
                "File {} is gzip-compressed but not BGZF (bgzip). Recompress with bgzip (or `bcftools view -Oz`) and index with `tabix -p vcf` / `bcftools index` as needed.",
                file.display()
            ));
        }
        return Err(crate::svx_error!(
            "File {} is not BGZF (bgzip) compressed",
            file.display()
        ));
    }

    let mut r = rust_htslib::bgzf::Reader::from_path(file)
        .map_err(|e| crate::svx_error!("Failed to open BGZF reader for {}: {e}", file.display()))?;
    let mut buf = vec![0u8; 512];

    let n = r.read(&mut buf).map_err(|e| {
        crate::svx_error!("Failed to read BGZF stream from {}: {e}", file.display())
    })?;
    buf.truncate(n);

    if buf.is_empty() {
        return Err(crate::svx_error!("File {} is empty", file.display()));
    }

    if looks_like_vcf(&buf) || looks_like_bcf(&buf) {
        Ok(())
    } else {
        Err(crate::svx_error!(
            "File {} is BGZF-compressed, but the decompressed header does not look like VCF or BCF",
            file.display()
        ))
    }
}

fn validate_indexed_vcf(file: &Path) -> Result<()> {
    if !is_indexed(file) {
        return Err(crate::svx_error!(
            "VCF file {} is not indexed (.tbi or .csi not found)",
            file.display()
        ));
    }
    validate_bgzip_vcf(file)
}

impl VcfReader {
    pub fn new(file: PathBuf, index: usize, io_tpool: Option<Arc<HtsThreadPool>>) -> Result<Self> {
        log::trace!("Start loading VCF {:?}", &file);
        validate_indexed_vcf(&file)
            .map_err(|e| crate::svx_error!("Error validating VCF: {}", e))?;

        let mut reader = bcf::IndexedReader::from_path(&file)
            .map_err(|e| crate::svx_error!("Failed to open VCF file {}: {}", file.display(), e))?;
        if let Some(ref pool) = io_tpool {
            unsafe {
                pool.attach_to_indexed_reader(&mut reader).map_err(|e| {
                    crate::svx_error!(
                        "Failed to attach shared HTS thread pool to VCF file {}: {}",
                        file.display(),
                        e
                    )
                })?;
                log::trace!("Attached shared HTS thread pool to VCF file {:?}", &file);
            }
        }

        let sample_n = reader.header().sample_count();

        log::trace!("{:?} samples n = {}", file.file_name().unwrap(), sample_n);

        if sample_n != 1 {
            return Err(crate::svx_error!(
                "Unsupported: VCF file {} must contain exactly 1 sample (n = {})",
                file.display(),
                sample_n
            ));
        }

        // TODO: Create a normalized struct for variant records
        let current_record = reader.empty_record();
        log::trace!("Finished loading VCF {:?}", &file);
        Ok(VcfReader {
            reader,
            current_record,
            sample_n,
            index,
            _io_tpool: io_tpool,
        })
    }

    pub fn advance(&mut self) -> Result<bool> {
        match self.reader.read(&mut self.current_record) {
            Some(Ok(())) => Ok(true),
            Some(Err(e)) => Err(crate::svx_error!(
                "Error reading record from VCF[{}]: {e}",
                self.index
            )),
            None => Ok(false),
        }
    }
}

pub struct VcfReaders {
    pub readers: Vec<VcfReader>,
    pub n: usize,
    pub io_tpool: Option<Arc<HtsThreadPool>>,
}

fn parse_contig_header_fields(
    id: Option<&str>,
    length: Option<&str>,
) -> Result<(String, Option<u64>)> {
    let id = id
        .ok_or_else(|| crate::svx_error!("Contig header is missing required 'ID' field"))?
        .to_string();
    let length = match length {
        Some(length) => Some(length.parse::<u64>().map_err(|e| {
            crate::svx_error!("Contig '{id}' has invalid 'length' value '{length}': {e}")
        })?),
        None => None,
    };

    Ok((id, length))
}

impl VcfReaders {
    pub fn new(vcf_files: Vec<PathBuf>, io_threads: usize) -> Result<Self> {
        let io_tpool = if io_threads >= 2 {
            let n_threads = i32::try_from(io_threads).map_err(|_| {
                crate::svx_error!(
                    "Invalid --io-threads value {}: exceeds supported i32 range",
                    io_threads
                )
            })?;
            Some(Arc::new(HtsThreadPool::new(n_threads)?))
        } else {
            None
        };

        let readers = vcf_files
            .into_iter()
            .enumerate()
            .map(|(index, file)| VcfReader::new(file, index, io_tpool.clone()))
            .collect::<Result<Vec<_>>>()?;

        let n = readers.len();

        Ok(VcfReaders {
            readers,
            n,
            io_tpool,
        })
    }

    pub fn get_contig_order(&self) -> Result<Vec<String>> {
        let mut contig_map: HashMap<String, HashSet<u64>> = HashMap::new();
        let mut contig_order = Vec::new();
        for reader in &self.readers {
            for record in reader.reader.header().header_records() {
                if let HeaderRecord::Contig { values, .. } = record {
                    let (id, length) = parse_contig_header_fields(
                        values.get("ID").map(String::as_str),
                        values.get("length").map(String::as_str),
                    )
                    .map_err(|e| {
                        crate::svx_error!(
                            "Malformed contig definition in VCF[{}]: {e}",
                            reader.index
                        )
                    })?;
                    let entry = contig_map.entry(id.clone()).or_insert_with(|| {
                        contig_order.push(id.clone());
                        HashSet::new()
                    });
                    if let Some(length) = length {
                        entry.insert(length);
                        if entry.len() > 1 {
                            return Err(crate::svx_error!(
                                "Inconsistent contig definitions found in VCF headers: Contig '{}' is defined with multiple lengths: {:?}",
                                id,
                                entry
                            ));
                        }
                    }
                }
            }
        }
        Ok(contig_order)
    }

    pub fn merge_headers(&self, dst_header: &mut Header) -> Result<SampleMapping> {
        let mut observed_sample_ids = HashSet::new();
        let mut sample_mapping = SampleMapping::new();
        let mut merged_pos = 0;

        for (vcf_index, reader) in self.readers.iter().enumerate() {
            let src_header = reader.reader.header();
            unsafe {
                dst_header.inner =
                    rust_htslib::htslib::bcf_hdr_merge(dst_header.inner, src_header.inner);
            }

            for (sample_index, sample_id) in src_header.samples().iter().enumerate() {
                if observed_sample_ids.contains(*sample_id) {
                    return Err(crate::svx_error!(
                        "Duplicate sample ID found: {}",
                        String::from_utf8_lossy(sample_id)
                    ));
                }
                observed_sample_ids.insert(sample_id.to_vec());
                dst_header.push_sample(sample_id);

                sample_mapping
                    .index_map
                    .insert((vcf_index, sample_index), merged_pos);
                sample_mapping
                    .reverse_map
                    .insert(merged_pos, (vcf_index, sample_index));
                merged_pos += 1;
            }
        }

        unsafe {
            rust_htslib::htslib::bcf_hdr_sync(dst_header.inner);
        }

        Ok(sample_mapping)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bcf::record::GenotypeAllele;
    use rust_htslib::bcf::Header;
    use std::sync::Arc;
    use tempdir::TempDir;

    fn create_test_reader_with_contig_record(
        samples: &[&str],
        contig_record: &[u8],
        index: usize,
    ) -> VcfReader {
        let temp_dir = TempDir::new("vcf_test").unwrap();
        let vcf_path = temp_dir.path().join("test.vcf");

        let mut header = Header::new();
        for sample in samples {
            header.push_sample(sample.as_bytes());
        }

        header.push_record(contig_record);
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);

        {
            let writer =
                rust_htslib::bcf::Writer::from_path(&vcf_path, &header, false, bcf::Format::Vcf)
                    .unwrap();

            let mut record = writer.empty_record();
            record.set_rid(Some(0));
            record.set_pos(100);
            record.set_alleles(&[b"A", b"T"]).unwrap();

            let sample_count = samples.len();
            let mut genotypes = Vec::new();
            for _ in 0..sample_count {
                genotypes.push(GenotypeAllele::Unphased(0));
                genotypes.push(GenotypeAllele::Unphased(1));
            }
            record.push_genotypes(&genotypes).unwrap();
        }

        {
            bcf::index::build(&vcf_path, None, 1, bcf::index::Type::Tbx).unwrap();
        }

        let reader = rust_htslib::bcf::IndexedReader::from_path(&vcf_path).unwrap();
        let record = reader.empty_record();

        VcfReader {
            reader,
            current_record: record,
            sample_n: samples.len() as u32,
            index,
            _io_tpool: None,
        }
    }

    fn create_test_reader(samples: &[&str], index: usize) -> VcfReader {
        create_test_reader_with_contig_record(
            samples,
            br#"##contig=<ID=1,length=249250621>"#,
            index,
        )
    }

    fn create_test_vcf_path(samples: &[&str], name: &str) -> (TempDir, PathBuf) {
        let temp_dir = TempDir::new("vcf_path_test").unwrap();
        let vcf_path = temp_dir.path().join(format!("{name}.vcf"));

        let mut header = Header::new();
        for sample in samples {
            header.push_sample(sample.as_bytes());
        }

        header.push_record(br#"##contig=<ID=1,length=249250621>"#);
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);

        {
            let writer =
                rust_htslib::bcf::Writer::from_path(&vcf_path, &header, false, bcf::Format::Vcf)
                    .unwrap();

            let mut record = writer.empty_record();
            record.set_rid(Some(0));
            record.set_pos(100);
            record.set_alleles(&[b"A", b"T"]).unwrap();
            if !samples.is_empty() {
                record
                    .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                    .unwrap();
            }
        }

        bcf::index::build(&vcf_path, None, 1, bcf::index::Type::Tbx).unwrap();
        (temp_dir, vcf_path)
    }

    #[test]
    fn test_merge_headers_duplicate_samples() {
        let readers = vec![
            create_test_reader(&["sample1", "sample2"], 0),
            create_test_reader(&["sample2", "sample3"], 1), // sample2 is duplicate
        ];

        let vcf_readers = VcfReaders {
            readers,
            n: 2,
            io_tpool: None,
        };
        let mut merged_header = Header::new();

        let error = vcf_readers
            .merge_headers(&mut merged_header)
            .expect_err("duplicate sample names should be rejected");
        // Duplicate sample names must always cause an error.
        assert!(error.to_string().contains("sample2"));
    }

    #[test]
    fn test_merge_headers_basic() {
        let readers = vec![
            create_test_reader(&["sample1", "sample2"], 0),
            create_test_reader(&["sample3"], 1),
            create_test_reader(&["sample4", "sample5", "sample6"], 2),
        ];

        let vcf_readers = VcfReaders {
            readers,
            n: 3,
            io_tpool: None,
        };
        let mut merged_header = Header::new();

        let sample_mapping = vcf_readers.merge_headers(&mut merged_header).unwrap();

        assert_eq!(sample_mapping.index_map.len(), 6);

        assert_eq!(sample_mapping.index_map.get(&(0, 0)), Some(&0)); // sample1 -> 0
        assert_eq!(sample_mapping.index_map.get(&(0, 1)), Some(&1)); // sample2 -> 1
        assert_eq!(sample_mapping.index_map.get(&(1, 0)), Some(&2)); // sample3 -> 2
        assert_eq!(sample_mapping.index_map.get(&(2, 0)), Some(&3)); // sample4 -> 3
        assert_eq!(sample_mapping.index_map.get(&(2, 1)), Some(&4)); // sample5 -> 4
        assert_eq!(sample_mapping.index_map.get(&(2, 2)), Some(&5)); // sample6 -> 5

        assert_eq!(sample_mapping.reverse_map.get(&0), Some(&(0, 0))); // 0 -> (vcf0, sample0)
        assert_eq!(sample_mapping.reverse_map.get(&1), Some(&(0, 1)));
        assert_eq!(sample_mapping.reverse_map.get(&2), Some(&(1, 0))); // 2 -> (vcf1, sample0)
        assert_eq!(sample_mapping.reverse_map.get(&3), Some(&(2, 0)));
        assert_eq!(sample_mapping.reverse_map.get(&4), Some(&(2, 1)));
        assert_eq!(sample_mapping.reverse_map.get(&5), Some(&(2, 2))); // 5 -> (vcf2, sample2)
    }

    #[test]
    fn test_vcf_readers_share_one_io_pool_when_io_threads_at_least_two() {
        let fixtures = [
            create_test_vcf_path(&["sample1"], "input_a"),
            create_test_vcf_path(&["sample2"], "input_b"),
        ];
        let vcf_files = fixtures
            .iter()
            .map(|(_, path)| path.clone())
            .collect::<Vec<_>>();

        let vcf_readers = VcfReaders::new(vcf_files, 2).unwrap();
        assert_eq!(vcf_readers.readers.len(), 2);
        assert!(vcf_readers.readers[0]._io_tpool.is_some());
        assert!(vcf_readers.readers[1]._io_tpool.is_some());

        let pool0 = vcf_readers.readers[0]._io_tpool.as_ref().unwrap();
        let pool1 = vcf_readers.readers[1]._io_tpool.as_ref().unwrap();
        assert!(Arc::ptr_eq(pool0, pool1));
    }

    #[test]
    fn test_vcf_readers_do_not_attach_io_pool_when_io_threads_is_one() {
        let fixtures = [
            create_test_vcf_path(&["sample1"], "input_c"),
            create_test_vcf_path(&["sample2"], "input_d"),
        ];
        let vcf_files = fixtures
            .iter()
            .map(|(_, path)| path.clone())
            .collect::<Vec<_>>();

        let vcf_readers = VcfReaders::new(vcf_files, 1).unwrap();
        assert_eq!(vcf_readers.readers.len(), 2);
        assert!(
            vcf_readers
                .readers
                .iter()
                .all(|reader| reader._io_tpool.is_none())
        );
    }

    #[test]
    fn test_parse_contig_header_fields_errors_when_id_is_missing() {
        let error = parse_contig_header_fields(None, Some("1000"))
            .expect_err("contig headers without ID should be rejected");
        assert!(error.to_string().contains("ID"));
    }

    #[test]
    fn test_parse_contig_header_fields_errors_when_contig_length_is_non_numeric() {
        let error = parse_contig_header_fields(Some("chr1"), Some("abc"))
            .expect_err("contig headers with non-numeric length should be rejected");
        assert!(error.to_string().contains("length"));
    }

    #[test]
    fn test_parse_contig_header_fields_treats_missing_length_as_unknown() {
        let parsed = parse_contig_header_fields(Some("chr1"), None)
            .expect("missing contig length should be treated as unknown");
        assert_eq!(parsed, ("chr1".to_string(), None));
    }

    #[test]
    fn test_get_contig_order_allows_mixed_known_and_unknown_contig_lengths() {
        let readers = vec![
            create_test_reader_with_contig_record(
                &["sample1"],
                br#"##contig=<ID=1,length=1000>"#,
                0,
            ),
            create_test_reader_with_contig_record(&["sample2"], br#"##contig=<ID=1>"#, 1),
        ];
        let vcf_readers = VcfReaders {
            readers,
            n: 2,
            io_tpool: None,
        };

        let contig_order = vcf_readers
            .get_contig_order()
            .expect("mixed known and unknown contig lengths should be accepted");

        assert_eq!(contig_order, vec!["1".to_string()]);
    }

    #[test]
    fn test_get_contig_order_rejects_multiple_known_contig_lengths_even_with_unknowns() {
        let readers = vec![
            create_test_reader_with_contig_record(
                &["sample1"],
                br#"##contig=<ID=1,length=1000>"#,
                0,
            ),
            create_test_reader_with_contig_record(&["sample2"], br#"##contig=<ID=1>"#, 1),
            create_test_reader_with_contig_record(
                &["sample3"],
                br#"##contig=<ID=1,length=2000>"#,
                2,
            ),
        ];
        let vcf_readers = VcfReaders {
            readers,
            n: 3,
            io_tpool: None,
        };

        let error = vcf_readers
            .get_contig_order()
            .expect_err("conflicting known contig lengths must still be rejected");

        assert!(error.to_string().contains("multiple lengths"));
    }
}
