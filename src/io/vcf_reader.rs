use crate::utils::util::Result;
use anyhow::anyhow;
use rust_htslib::bcf::{self, Header, HeaderRecord, Read};
use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader, Seek},
    path::{Path, PathBuf},
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
    pub header: bcf::header::HeaderView,
    pub current_record: bcf::Record,
    pub sample_n: u32,
    pub index: usize,
}

const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

fn add_extension(path: &Path, ext: &str) -> PathBuf {
    let mut file_name = path.file_name().unwrap().to_os_string();
    file_name.push(".");
    file_name.push(ext);
    path.with_file_name(file_name)
}

fn is_indexed(file: &Path) -> bool {
    let csi_index = add_extension(file, "csi");
    let tbi_index = add_extension(file, "tbi");
    csi_index.exists() || tbi_index.exists()
}

// TODO: Extend logic to also deal with bcf files
fn validate_vcf_file(file: &Path) -> Result<()> {
    let mut fp = std::fs::File::open(file).map_err(|e| anyhow!("Failed to open file: {}", e))?;
    let mut buffer = [0; 2];
    std::io::Read::read_exact(&mut fp, &mut buffer)
        .map_err(|e| anyhow!("Failed to read file: {}", e))?;

    if buffer != GZIP_MAGIC_NUMBER {
        return Err(anyhow!("File {} is not gzip compressed", file.display()));
    }

    fp.rewind()
        .map_err(|e| anyhow!("Failed to rewind file: {}", e))?;

    let gz = flate2::read::GzDecoder::new(fp);
    let mut reader = BufReader::new(gz);
    let mut first_line = String::new();
    reader
        .read_line(&mut first_line)
        .map_err(|e| anyhow!("Failed to read file: {}", e))?;

    if !first_line.trim().starts_with("##fileformat=VCFv") {
        return Err(anyhow!("File {} is not a valid VCF file", file.display()));
    }

    Ok(())
}

impl VcfReader {
    pub fn new(file: PathBuf, index: usize) -> Result<Self> {
        log::trace!("Start loading VCF {:?}", &file);
        validate_vcf_file(&file).map_err(|e| anyhow!("Error validating VCF: {}", e))?;

        if !is_indexed(&file) {
            return Err(anyhow!("VCF file {} is not indexed", file.display()));
        }

        let reader = bcf::IndexedReader::from_path(&file)
            .map_err(|e| anyhow!("Failed to open VCF file {}: {}", file.display(), e))?;
        let header = reader.header().clone();
        let sample_n = header.sample_count();

        log::trace!("{:?} samples n = {}", file.file_name().unwrap(), sample_n);

        if sample_n > 1 {
            return Err(anyhow!(
                "Unsupported: VCF file {} has more than 1 sample (n = {})!",
                file.display(),
                sample_n
            ));
        }

        // TODO: Create a normalized struct for variant records
        let current_record = reader.empty_record();
        log::trace!("Finished loading VCF {:?}", &file);
        Ok(VcfReader {
            reader,
            header,
            current_record,
            sample_n,
            index,
        })
    }

    pub fn advance(&mut self) -> bool {
        match self.reader.read(&mut self.current_record) {
            Some(Ok(())) => true,
            Some(Err(_)) | None => false,
        }
    }
}

pub struct VcfReaders {
    pub readers: Vec<VcfReader>,
    pub n: usize,
}

impl VcfReaders {
    pub fn new(vcf_files: Vec<PathBuf>) -> Result<Self> {
        let readers = vcf_files
            .into_iter()
            .enumerate()
            .map(|(index, file)| VcfReader::new(file, index))
            .collect::<Result<Vec<_>>>()?;

        let n = readers.len();

        Ok(VcfReaders { readers, n })
    }

    pub fn get_contig_order(&self) -> Result<Vec<String>> {
        let mut contig_map: HashMap<String, HashSet<u64>> = HashMap::new();
        let mut contig_order = Vec::new();
        for reader in &self.readers {
            for record in reader.header.header_records() {
                if let HeaderRecord::Contig { values, .. } = record {
                    let id = values.get("ID").unwrap().to_string();
                    let length = values
                        .get("length")
                        .and_then(|l| l.parse::<u64>().ok())
                        .unwrap_or(0);
                    let entry = contig_map.entry(id.clone()).or_insert_with(|| {
                        contig_order.push(id.clone());
                        HashSet::new()
                    });
                    entry.insert(length);
                }
            }
        }
        for id in &contig_order {
            let lengths = contig_map.get(id).unwrap();
            if lengths.len() > 1 {
                return Err(anyhow!(
                    "Inconsistent contig definitions found in VCF headers: Contig '{}' is defined with multiple lengths: {:?}",
                    id, lengths
                ));
            }
        }
        Ok(contig_order)
    }

    pub fn merge_headers(
        &self,
        dst_header: &mut Header,
        force_samples: bool,
    ) -> Result<SampleMapping> {
        let mut observed_sample_ids = HashSet::new();
        let mut sample_mapping = SampleMapping::new();
        let mut merged_pos = 0;

        for (vcf_index, reader) in self.readers.iter().enumerate() {
            let src_header = &reader.header;
            unsafe {
                // Merge the current vcf_header into the output header
                dst_header.inner =
                    rust_htslib::htslib::bcf_hdr_merge(dst_header.inner, src_header.inner);
            }

            for (sample_index, sample_id) in src_header.samples().iter().enumerate() {
                if observed_sample_ids.contains(*sample_id) {
                    if force_samples {
                        continue; // If forcing samples, skip duplicates
                    } else {
                        return Err(anyhow!(
                            "Duplicate sample ID found: {}",
                            String::from_utf8_lossy(sample_id)
                        ));
                    }
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
    use tempdir::TempDir;

    fn create_test_reader(samples: &[&str], index: usize) -> VcfReader {
        let temp_dir = TempDir::new("vcf_test").unwrap();
        let vcf_path = temp_dir.path().join("test.vcf");

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
        let header_view = reader.header().clone();
        let record = reader.empty_record();

        VcfReader {
            reader,
            header: header_view,
            current_record: record,
            sample_n: samples.len() as u32,
            index,
        }
    }

    #[test]
    fn test_merge_headers_duplicate_samples() {
        let readers = vec![
            create_test_reader(&["sample1", "sample2"], 0),
            create_test_reader(&["sample2", "sample3"], 1), // sample2 is duplicate
        ];

        let vcf_readers = VcfReaders { readers, n: 2 };
        let mut merged_header = Header::new();

        // Test that duplicate samples cause an error when force_samples is false
        assert!(vcf_readers
            .merge_headers(&mut merged_header, false)
            .is_err());
    }

    #[test]
    fn test_merge_headers_basic() {
        let readers = vec![
            create_test_reader(&["sample1", "sample2"], 0),
            create_test_reader(&["sample3"], 1),
            create_test_reader(&["sample4", "sample5", "sample6"], 2),
        ];

        let vcf_readers = VcfReaders { readers, n: 3 };
        let mut merged_header = Header::new();

        let sample_mapping = vcf_readers
            .merge_headers(&mut merged_header, false)
            .unwrap();

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
}
