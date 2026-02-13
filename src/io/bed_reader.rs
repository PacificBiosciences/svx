use super::readers::open_catalog_reader;
use crate::{
    core::containers::interval_tree::{Interval, IntervalTree},
    utils::util::Result,
};
use std::{collections::HashMap, io::BufRead, path::Path};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TrId {
    pub id: String,
    pub motif_len: usize,
}

pub fn line_to_interval(line: &str) -> Result<(String, Interval<u32, TrId>)> {
    const EXPECTED_FIELD_COUNT: usize = 4;
    let split_line: Vec<&str> = line.split_whitespace().collect();
    if split_line.len() != EXPECTED_FIELD_COUNT {
        return Err(crate::svx_error!(
            "Expected {} fields in the format 'chrom start end info', found {}: {}",
            EXPECTED_FIELD_COUNT,
            split_line.len(),
            line
        ));
    }

    let (chrom, start, end, info_fields) = match &split_line[..] {
        [chrom, start, end, info_fields] => (*chrom, *start, *end, *info_fields),
        _ => unreachable!(),
    };

    let start: u32 = start
        .parse()
        .map_err(|e| crate::svx_error!("Invalid start position: {}", e))?;
    let end: u32 = end
        .parse()
        .map_err(|e| crate::svx_error!("Invalid end position: {}", e))?;

    let fields = decode_fields(info_fields)?;
    let id = get_field(&fields, "ID")?;
    let motifs: Vec<String> = get_field(&fields, "MOTIFS")?
        .split(',')
        .map(|s| s.to_string())
        .collect();
    let motif_len = motifs.iter().map(|s| s.len()).min().unwrap();

    let info = TrId { id, motif_len };
    Ok((chrom.to_owned(), Interval::new(start, end, info)))
}

pub fn get_field(fields: &HashMap<&str, String>, key: &str) -> Result<String> {
    fields
        .get(key)
        .ok_or_else(|| crate::svx_error!("{} field missing", key))
        .map(|s| s.to_string())
}

pub fn decode_fields(info_fields: &str) -> Result<HashMap<&str, String>> {
    let mut fields = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding)?;
        if fields.insert(name, value.to_string()).is_some() {
            return Err(crate::svx_error!("Duplicate field name: '{}'", name));
        }
    }
    Ok(fields)
}

fn decode_info_field(encoding: &str) -> Result<(&str, &str)> {
    let parts: Vec<&str> = encoding.splitn(2, '=').collect();
    if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
        Err(crate::svx_error!(
            "Field must be in 'name=value' format: '{}'",
            encoding
        ))
    } else {
        Ok((parts[0], parts[1]))
    }
}

pub type BedIntervalTree = IntervalTree<u32, TrId>;

pub struct BedMap {
    pub interval_map: HashMap<String, BedIntervalTree>,
}

impl BedMap {
    pub fn new<P: AsRef<Path>>(bed_path: P) -> Result<Self> {
        let reader = open_catalog_reader(bed_path.as_ref())?;
        let mut interval_map: HashMap<String, Vec<Interval<u32, TrId>>> = HashMap::new();

        for (line_number, result_line) in reader.lines().enumerate() {
            let line = result_line
                .map_err(|e| crate::svx_error!("Error at BED line {}: {}", line_number + 1, e))?;
            let (chrom, interval) = line_to_interval(&line)
                .map_err(|e| crate::svx_error!("Error at BED line {}: {}", line_number + 1, e))?;

            interval_map.entry(chrom).or_default().push(interval);
        }

        let interval_trees: HashMap<String, BedIntervalTree> = interval_map
            .into_iter()
            .map(|(chrom, intervals)| (chrom, IntervalTree::new(intervals)))
            .collect();

        Ok(BedMap {
            interval_map: interval_trees,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_line_to_interval_valid() {
        let line_with_id = "chr1\t100\t200\tID=RFC1;MOTIFS=abc";
        let (_, interval_with_id) = line_to_interval(line_with_id).unwrap();
        assert_eq!(interval_with_id.value.id, "RFC1");
    }

    #[test]
    fn test_line_to_interval_invalid() {
        let line_too_short = "chr1\t100";
        assert!(line_to_interval(line_too_short).is_err());
        let line_invalid_start = "chr1\tabc\t200";
        assert!(line_to_interval(line_invalid_start).is_err());
        let line_invalid_end = "chr1\t100\txyz";
        assert!(line_to_interval(line_invalid_end).is_err());
    }

    #[test]
    fn test_bed_map_new() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "chr1\t10\t20\tID=data1;MOTIFS=ATC")?;
        writeln!(temp_file, "chr2\t30\t40\tID=data2;MOTIFS=ATC")?;
        writeln!(temp_file, "chr1\t50\t60\tID=data3;MOTIFS=ATC")?;
        temp_file.flush()?;
        let bed_map = BedMap::new(temp_file.path())?;
        assert_eq!(bed_map.interval_map.len(), 2);

        let chr1_tree = bed_map.interval_map.get("chr1").unwrap();
        let overlapping_chr1 = chr1_tree.find_overlapping(15, 55);
        assert_eq!(overlapping_chr1.len(), 2);
        assert!(
            overlapping_chr1
                .iter()
                .any(|iv| iv.start == 10 && iv.stop == 20 && iv.value.id == "data1")
        );
        assert!(
            overlapping_chr1
                .iter()
                .any(|iv| iv.start == 50 && iv.stop == 60 && iv.value.id == "data3")
        );

        let contained_chr1 = chr1_tree.find_contained(5, 65);
        assert_eq!(contained_chr1.len(), 2);

        let chr2_tree = bed_map.interval_map.get("chr2").unwrap();
        let overlapping_chr2 = chr2_tree.find_overlapping(35, 36);
        assert_eq!(overlapping_chr2.len(), 1);
        assert_eq!(overlapping_chr2[0].start, 30);
        assert_eq!(overlapping_chr2[0].stop, 40);
        assert_eq!(overlapping_chr2[0].value.id, "data2");

        Ok(())
    }

    #[test]
    fn test_bed_map_new_empty_file() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let bed_map = BedMap::new(temp_file.path())?;
        assert!(bed_map.interval_map.is_empty());
        Ok(())
    }
}
