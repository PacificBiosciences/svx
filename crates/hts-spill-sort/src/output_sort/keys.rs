use super::spill_format::SpillFormatEntry;
use crate::{Result, SpillSortError};
use rust_htslib::bcf;
use std::{cmp::Ordering, sync::Arc};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PackedSortKey {
    rid: u32,
    pos: i64,
    packed_alleles: Arc<[u8]>,
}

impl PackedSortKey {
    pub fn from_record(record: &bcf::Record) -> Result<Self> {
        let rid = record
            .rid()
            .ok_or_else(|| SpillSortError::message("Sort key requires record RID".to_string()))?;
        let pos = record.pos();
        let packed_alleles = pack_alleles_case_sensitive(record.alleles().as_slice())?;
        Ok(Self {
            rid,
            pos,
            packed_alleles,
        })
    }

    pub fn packed_allele_bytes_len(&self) -> usize {
        self.packed_alleles.len()
    }

    pub fn rid(&self) -> u32 {
        self.rid
    }

    pub fn pos(&self) -> i64 {
        self.pos
    }

    pub fn packed_alleles(&self) -> Arc<[u8]> {
        self.packed_alleles.clone()
    }
}

impl Ord for PackedSortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid
            .cmp(&other.rid)
            .then_with(|| self.pos.cmp(&other.pos))
            .then_with(|| compare_packed_alleles(&self.packed_alleles, &other.packed_alleles))
    }
}

impl PartialOrd for PackedSortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn pack_alleles_case_sensitive(alleles: &[&[u8]]) -> Result<Arc<[u8]>> {
    let capacity = alleles.iter().fold(0usize, |acc, allele| {
        acc.saturating_add(4).saturating_add(allele.len())
    });
    let mut packed = Vec::with_capacity(capacity);
    for allele in alleles {
        pack_single_allele_case_sensitive(&mut packed, allele)?;
    }
    Ok(Arc::from(packed))
}

fn pack_single_allele_case_sensitive(packed: &mut Vec<u8>, allele: &[u8]) -> Result<()> {
    let len = u32::try_from(allele.len()).map_err(|_| {
        SpillSortError::message(format!(
            "Sort key allele length {} exceeds u32::MAX",
            allele.len()
        ))
    })?;
    packed.extend_from_slice(&len.to_be_bytes());
    packed.extend_from_slice(allele);
    Ok(())
}

fn compare_packed_alleles(mut left: &[u8], mut right: &[u8]) -> Ordering {
    while !left.is_empty() {
        if right.is_empty() {
            return Ordering::Greater;
        }
        let (left_allele, left_remaining) = unpack_packed_allele(left);
        let (right_allele, right_remaining) = unpack_packed_allele(right);
        let cmp = left_allele.cmp(right_allele);
        if cmp != Ordering::Equal {
            return cmp;
        }
        left = left_remaining;
        right = right_remaining;
    }
    if right.is_empty() {
        Ordering::Equal
    } else {
        Ordering::Less
    }
}

fn unpack_packed_allele(packed: &[u8]) -> (&[u8], &[u8]) {
    assert!(
        packed.len() >= 4,
        "packed sort key must include allele length bytes"
    );
    let len = {
        let bytes = [packed[0], packed[1], packed[2], packed[3]];
        u32::from_be_bytes(bytes) as usize
    };
    let payload = &packed[4..];
    assert!(
        payload.len() >= len,
        "packed sort key allele payload length is truncated"
    );
    payload.split_at(len)
}

pub fn packed_sort_key_from_spill_entry(entry: &SpillFormatEntry) -> PackedSortKey {
    PackedSortKey {
        rid: entry.rid,
        pos: entry.pos,
        packed_alleles: entry.packed_alleles.clone(),
    }
}

pub fn estimate_record_bytes(key: &PackedSortKey, record: &bcf::Record) -> usize {
    128usize
        .saturating_add(key.packed_allele_bytes_len())
        .saturating_add(record.id().len())
        .saturating_add(8)
}

pub fn ensure_records_sorted_for_output(records: &mut Vec<bcf::Record>) -> Result<()> {
    if records.len() < 2 {
        return Ok(());
    }

    let keys = records
        .iter()
        .map(PackedSortKey::from_record)
        .collect::<Result<Vec<_>>>()?;
    if keys.windows(2).all(|pair| pair[0] <= pair[1]) {
        return Ok(());
    }

    let mut keyed_records = records
        .drain(..)
        .zip(keys)
        .enumerate()
        .map(|(record_ordinal, (record, key))| (key, record_ordinal, record))
        .collect::<Vec<_>>();
    keyed_records.sort_by(|(left_key, left_idx, _), (right_key, right_idx, _)| {
        left_key
            .cmp(right_key)
            .then_with(|| left_idx.cmp(right_idx))
    });
    records.extend(keyed_records.into_iter().map(|(_, _, record)| record));
    Ok(())
}

#[cfg(test)]
#[derive(Clone, Debug, Eq, PartialEq)]
struct SortKey {
    rid: u32,
    pos: i64,
    alleles: Vec<Vec<u8>>,
}

#[cfg(test)]
impl SortKey {
    fn from_record(record: &bcf::Record) -> Result<Self> {
        let rid = record
            .rid()
            .ok_or_else(|| SpillSortError::message("Sort key requires record RID".to_string()))?;
        let pos = record.pos();
        let alleles = record
            .alleles()
            .into_iter()
            .map(|allele| allele.to_vec())
            .collect::<Vec<_>>();
        Ok(Self { rid, pos, alleles })
    }
}

#[cfg(test)]
impl Ord for SortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid
            .cmp(&other.rid)
            .then_with(|| self.pos.cmp(&other.pos))
            .then_with(|| compare_allele_vectors(&self.alleles, &other.alleles))
    }
}

#[cfg(test)]
impl PartialOrd for SortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
fn compare_allele_vectors(left: &[Vec<u8>], right: &[Vec<u8>]) -> Ordering {
    for (idx, left_allele) in left.iter().enumerate() {
        if idx >= right.len() {
            return Ordering::Greater;
        }
        let cmp = compare_alleles_case_sensitive(left_allele.as_slice(), right[idx].as_slice());
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    left.len().cmp(&right.len())
}

#[cfg(test)]
fn compare_alleles_case_sensitive(left: &[u8], right: &[u8]) -> Ordering {
    left.cmp(right)
}

#[cfg(test)]
mod tests {
    use super::{PackedSortKey, SortKey, ensure_records_sorted_for_output};
    use rust_htslib::bcf;

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

    #[test]
    fn compare_sort_key_matches_bcftools_packed_sort_semantics() {
        let factory = record_factory();
        let first = make_record(&factory, 0, 10, "r1", "a", "G");
        let second = make_record(&factory, 0, 10, "r2", "A", "g");
        let third = make_record(&factory, 0, 10, "r3", "A", "GA");

        let first_key = SortKey::from_record(&first).expect("sort key should be derived");
        let second_key = SortKey::from_record(&second).expect("sort key should be derived");
        let third_key = SortKey::from_record(&third).expect("sort key should be derived");

        assert_eq!(first_key.cmp(&second_key), std::cmp::Ordering::Greater);
        assert_eq!(second_key.cmp(&third_key), std::cmp::Ordering::Greater);
    }

    #[test]
    fn packed_sort_key_ordering_matches_sort_key_ordering() {
        let factory = record_factory();
        let first = make_record(&factory, 0, 10, "r1", "A", "C");
        let second = make_record(&factory, 0, 10, "r2", "A", "GA");
        let third = make_record(&factory, 1, 5, "r3", "T", "G");

        let first_key = SortKey::from_record(&first).expect("sort key should be derived");
        let second_key = SortKey::from_record(&second).expect("sort key should be derived");
        let third_key = SortKey::from_record(&third).expect("sort key should be derived");

        let first_packed =
            PackedSortKey::from_record(&first).expect("packed sort key should be derived");
        let second_packed =
            PackedSortKey::from_record(&second).expect("packed sort key should be derived");
        let third_packed =
            PackedSortKey::from_record(&third).expect("packed sort key should be derived");

        assert_eq!(first_key.cmp(&second_key), first_packed.cmp(&second_packed));
        assert_eq!(first_key.cmp(&third_key), first_packed.cmp(&third_packed));
        assert_eq!(second_key.cmp(&third_key), second_packed.cmp(&third_packed));
    }

    #[test]
    fn ensure_records_sorted_for_output_reorders_non_monotonic_batch() {
        let factory = record_factory();
        let mut records = vec![
            make_record(&factory, 1, 20, "high_rid", "A", "C"),
            make_record(&factory, 0, 10, "low_rid", "A", "C"),
        ];

        ensure_records_sorted_for_output(&mut records)
            .expect("record sorting should succeed for non-monotonic batch");

        let ids = records
            .iter()
            .map(|record| {
                String::from_utf8(record.id().to_vec()).expect("record ID should be UTF-8")
            })
            .collect::<Vec<_>>();
        assert_eq!(ids, vec!["low_rid", "high_rid"]);
    }
}
