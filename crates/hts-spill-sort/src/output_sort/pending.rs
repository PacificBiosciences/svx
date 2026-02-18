use super::keys::PackedSortKey;
use crate::{Result, SpillSortError};
use rust_htslib::bcf;
use std::{cmp::Ordering, collections::BinaryHeap};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct SourceOrder {
    pub blob_ordinal: u64,
    pub record_ordinal: usize,
}

#[derive(Debug)]
pub struct SortEntry {
    pub key: PackedSortKey,
    pub source: SourceOrder,
    pub record: bcf::Record,
}

#[derive(Debug)]
pub struct PendingBlobRun {
    pub entries: Vec<SortEntry>,
    pub cursor: usize,
}

impl PendingBlobRun {
    fn current(&self) -> Option<&SortEntry> {
        self.entries.get(self.cursor)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct PackedPendingHeapItem {
    key: PackedSortKey,
    source: SourceOrder,
    run_index: usize,
}

impl Ord for PackedPendingHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .cmp(&self.key)
            .then_with(|| other.source.cmp(&self.source))
    }
}

impl PartialOrd for PackedPendingHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn drain_pending_runs_into_sink<F>(
    pending_runs: &mut Vec<PendingBlobRun>,
    pending_bytes: &mut usize,
    mut emit: F,
) -> Result<()>
where
    F: FnMut(&SortEntry) -> Result<()>,
{
    if pending_runs.is_empty() {
        return Ok(());
    }

    let mut heap = BinaryHeap::new();
    for (run_index, run) in pending_runs.iter().enumerate() {
        if let Some(entry) = run.current() {
            heap.push(PackedPendingHeapItem {
                key: entry.key.clone(),
                source: entry.source,
                run_index,
            });
        }
    }

    while let Some(item) = heap.pop() {
        let run = &mut pending_runs[item.run_index];
        let entry = run.entries.get(run.cursor).ok_or_else(|| {
            SpillSortError::message(
                "Pending run cursor out of bounds during pre-spill interleaving".to_string(),
            )
        })?;
        emit(entry)?;
        run.cursor = run.cursor.saturating_add(1);

        if let Some(next) = run.current() {
            heap.push(PackedPendingHeapItem {
                key: next.key.clone(),
                source: next.source,
                run_index: item.run_index,
            });
        }
    }

    pending_runs.clear();
    *pending_bytes = 0;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{PackedPendingHeapItem, SourceOrder};
    use crate::output_sort::keys::PackedSortKey;
    use rust_htslib::bcf;
    use std::collections::BinaryHeap;

    fn test_header() -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(br#"##fileformat=VCFv4.2"#);
        header.push_record(br#"##contig=<ID=chr1,length=1000>"#);
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="">"#);
        header.push_sample(b"S1");
        header
    }

    fn record_factory() -> bcf::Writer {
        bcf::Writer::from_stdout(&test_header(), true, bcf::Format::Vcf)
            .expect("record factory should initialize")
    }

    fn make_record(factory: &bcf::Writer, id: &str) -> bcf::Record {
        let mut record = factory.empty_record();
        record.set_rid(Some(0));
        record.set_pos(10);
        record
            .set_id(id.as_bytes())
            .expect("record id should be set for sorter test");
        record
            .set_alleles(&[b"A", b"C"])
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
    fn packed_pending_heap_tie_breaks_by_source_order() {
        let factory = record_factory();
        let first = make_record(&factory, "r1");
        let second = make_record(&factory, "r2");
        let first_packed =
            PackedSortKey::from_record(&first).expect("packed sort key should be derived");
        let second_packed =
            PackedSortKey::from_record(&second).expect("packed sort key should be derived");

        let lower_source = PackedPendingHeapItem {
            key: first_packed,
            source: SourceOrder {
                blob_ordinal: 0,
                record_ordinal: 0,
            },
            run_index: 0,
        };
        let higher_source = PackedPendingHeapItem {
            key: second_packed,
            source: SourceOrder {
                blob_ordinal: 1,
                record_ordinal: 0,
            },
            run_index: 1,
        };

        let mut heap = BinaryHeap::new();
        heap.push(higher_source);
        heap.push(lower_source.clone());

        let first_item = heap.pop().expect("heap should contain two items");
        assert_eq!(first_item.source, lower_source.source);
    }
}
