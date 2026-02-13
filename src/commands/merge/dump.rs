use crate::{
    core::{svtype::SvType, variant::VariantInternal, variant_block::VariantBlockResult},
    utils::util::Result,
};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::{Arc, Mutex, MutexGuard},
};

const DUMP_HEADER: &str = "record_type\tcontig\tvariant_type\tgroup_size\tgroup_count\tgroup_index\tvariant_index\tvariant_id\tvariant_vcf_id\tsample_id\tstart\tend\tsvlen\ttrid\tkd_x\tkd_y";

pub(crate) type SharedDumpWriter = Arc<DumpWriter>;

pub(crate) struct DumpWriter {
    writer: Mutex<BufWriter<File>>,
}

impl DumpWriter {
    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        let file = File::create(path).map_err(|error| {
            crate::svx_error!(
                "Failed to create dump file at {}: {}",
                path.display(),
                error
            )
        })?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{DUMP_HEADER}")
            .map_err(|error| crate::svx_error!("Failed to write dump header: {error}"))?;
        Ok(Self {
            writer: Mutex::new(writer),
        })
    }

    pub(crate) fn dump_kd_coordinates(
        &self,
        contig: &str,
        variant_type: SvType,
        variants: &[VariantInternal],
    ) -> Result<()> {
        let contig = sanitize_tsv_field(contig);
        let variant_type = variant_type.to_string();
        let mut writer = self.lock_writer()?;
        for (variant_index, variant) in variants.iter().enumerate() {
            let point = variant.point();
            write_dump_row(
                &mut writer,
                [
                    "kd_coord".to_string(),
                    contig.clone(),
                    variant_type.clone(),
                    String::new(),
                    String::new(),
                    String::new(),
                    variant_index.to_string(),
                    sanitize_tsv_field(&variant.id),
                    variant.vcf_id.to_string(),
                    variant.sample_id.to_string(),
                    variant.start.to_string(),
                    variant.end.to_string(),
                    variant.svlen.to_string(),
                    variant
                        .trid
                        .as_ref()
                        .map_or_else(String::new, |trid| sanitize_tsv_field(&trid.id)),
                    point[0].to_string(),
                    point[1].to_string(),
                ],
            )?;
        }
        writer
            .flush()
            .map_err(|error| crate::svx_error!("Failed to flush dump rows: {error}"))?;
        Ok(())
    }

    pub(crate) fn dump_group_stats(&self, variant_block_result: &VariantBlockResult) -> Result<()> {
        let contig = sanitize_tsv_field(&variant_block_result.contig);
        let variant_type = variant_block_result.variant_type.to_string();
        let mut group_size_counts: HashMap<usize, usize> = HashMap::new();
        let mut writer = self.lock_writer()?;

        for (group_index, group) in variant_block_result.groups.iter().enumerate() {
            if group.is_empty() {
                continue;
            }
            let group_size = group.len();
            *group_size_counts.entry(group_size).or_insert(0) += 1;
            for variant in group {
                write_dump_row(
                    &mut writer,
                    [
                        "group_variant".to_string(),
                        contig.clone(),
                        variant_type.clone(),
                        group_size.to_string(),
                        String::new(),
                        group_index.to_string(),
                        variant.index.to_string(),
                        sanitize_tsv_field(&variant.id),
                        variant.vcf_id.to_string(),
                        variant.sample_id.to_string(),
                        variant.start.to_string(),
                        variant.end.to_string(),
                        variant.svlen.to_string(),
                        variant
                            .trid
                            .as_ref()
                            .map_or_else(String::new, |trid| sanitize_tsv_field(&trid.id)),
                        String::new(),
                        String::new(),
                    ],
                )?;
            }
        }

        let mut sorted_group_size_counts: Vec<_> = group_size_counts.into_iter().collect();
        sorted_group_size_counts.sort_unstable_by_key(|(group_size, _)| *group_size);

        for (group_size, group_count) in sorted_group_size_counts {
            write_dump_row(
                &mut writer,
                [
                    "group_size_count".to_string(),
                    contig.clone(),
                    variant_type.clone(),
                    group_size.to_string(),
                    group_count.to_string(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                ],
            )?;
        }

        writer
            .flush()
            .map_err(|error| crate::svx_error!("Failed to flush dump rows: {error}"))?;
        Ok(())
    }

    fn lock_writer(&self) -> Result<MutexGuard<'_, BufWriter<File>>> {
        self.writer
            .lock()
            .map_err(|_| crate::svx_error!("Dump file writer lock poisoned"))
    }
}

fn sanitize_tsv_field(input: &str) -> String {
    input
        .replace('\t', "\\t")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

fn write_dump_row(writer: &mut BufWriter<File>, fields: [String; 16]) -> Result<()> {
    writeln!(writer, "{}", fields.join("\t"))
        .map_err(|error| crate::svx_error!("Failed to write dump row: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{svtype::SvType, variant::test_utils, variant_block::VariantBlockResult};
    use tempfile::NamedTempFile;

    #[test]
    fn dump_writer_writes_kd_coordinates_and_group_stats() {
        let dump_file = NamedTempFile::new().expect("dump file should be created");
        let writer =
            DumpWriter::from_path(dump_file.path()).expect("dump writer should initialize");

        let mut variant_a =
            test_utils::from_parts(0, "variant_a".to_string(), SvType::INSERTION, 10.0, 11.0)
                .expect("variant a should be created");
        variant_a.index = 0;
        let mut variant_b =
            test_utils::from_parts(1, "variant_b".to_string(), SvType::INSERTION, 20.0, 22.0)
                .expect("variant b should be created");
        variant_b.index = 1;

        writer
            .dump_kd_coordinates(
                "chr1",
                SvType::INSERTION,
                &[variant_a.clone(), variant_b.clone()],
            )
            .expect("KD coordinates should be dumped");
        writer
            .dump_group_stats(&VariantBlockResult {
                groups: vec![vec![variant_a, variant_b], vec![]],
                contig: "chr1".to_string(),
                variant_type: SvType::INSERTION,
                n: 2,
            })
            .expect("group stats should be dumped");
        drop(writer);

        let contents =
            std::fs::read_to_string(dump_file.path()).expect("dump file should be readable");
        assert!(contents.starts_with(DUMP_HEADER));
        assert!(contents.contains("kd_coord\tchr1\tINS"));
        assert!(contents.contains("group_variant\tchr1\tINS\t2"));
        assert!(contents.contains("group_size_count\tchr1\tINS\t2\t1"));
    }
}
