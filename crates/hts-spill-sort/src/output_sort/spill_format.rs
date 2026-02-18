use crate::Result;
use std::sync::Arc;

pub const SPILL_MAGIC: [u8; 4] = *b"HTSS";
pub const SPILL_VERSION: u16 = 2;
pub const SPILL_FRAME_HEADER_LEN: usize = 14;

const MAGIC_OFFSET: usize = 0;
const VERSION_OFFSET: usize = MAGIC_OFFSET + SPILL_MAGIC.len();
const PAYLOAD_LEN_OFFSET: usize = VERSION_OFFSET + std::mem::size_of::<u16>();
const FIXED_FIELDS_LEN: usize = std::mem::size_of::<u32>() // rid
    + std::mem::size_of::<i64>() // pos
    + std::mem::size_of::<u64>() // blob_ordinal
    + std::mem::size_of::<u64>() // record_ordinal
    + std::mem::size_of::<u64>() // run_id
    + std::mem::size_of::<u32>() // packed allele bytes len
    + std::mem::size_of::<u32>(); // record payload len

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SpillSourceOrder {
    pub blob_ordinal: u64,
    pub record_ordinal: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SpillFormatEntry {
    pub rid: u32,
    pub pos: i64,
    pub packed_alleles: Arc<[u8]>,
    pub source: SpillSourceOrder,
    pub run_id: u64,
    pub record_payload: Vec<u8>,
}

#[cfg(test)]
pub fn encode_entry(entry: &SpillFormatEntry) -> Result<Vec<u8>> {
    let mut encoded = Vec::new();
    encode_entry_into(entry, &mut encoded)?;
    Ok(encoded)
}

#[cfg(test)]
pub fn encode_entry_into(entry: &SpillFormatEntry, encoded: &mut Vec<u8>) -> Result<()> {
    encode_entry_with_run_id_into(entry, entry.run_id, encoded)
}

pub fn encode_entry_with_run_id_into(
    entry: &SpillFormatEntry,
    run_id: u64,
    encoded: &mut Vec<u8>,
) -> Result<()> {
    let packed_alleles_len = u32::try_from(entry.packed_alleles.len()).map_err(|_| {
        crate::SpillSortError::message(format!(
            "spill entry packed-alleles length {} exceeds u32::MAX",
            entry.packed_alleles.len()
        ))
    })?;
    let record_payload_len = u32::try_from(entry.record_payload.len()).map_err(|_| {
        crate::SpillSortError::message(format!(
            "spill entry record payload length {} exceeds u32::MAX",
            entry.record_payload.len()
        ))
    })?;
    let record_ordinal = u64::try_from(entry.source.record_ordinal).map_err(|_| {
        crate::SpillSortError::message(format!(
            "spill entry source record ordinal {} exceeds u64::MAX",
            entry.source.record_ordinal
        ))
    })?;
    let packed_alleles_len_usize = usize::try_from(packed_alleles_len).map_err(|_| {
        crate::SpillSortError::message(
            "spill entry packed-alleles length exceeds usize::MAX".to_string(),
        )
    })?;
    let record_payload_len_usize = usize::try_from(record_payload_len).map_err(|_| {
        crate::SpillSortError::message(
            "spill entry record payload length exceeds usize::MAX".to_string(),
        )
    })?;
    let payload_len = FIXED_FIELDS_LEN
        .saturating_add(packed_alleles_len_usize)
        .saturating_add(record_payload_len_usize);
    let payload_len_u64 = u64::try_from(payload_len).map_err(|_| {
        crate::SpillSortError::message(
            "spill entry payload length exceeds u64::MAX when encoding".to_string(),
        )
    })?;

    encoded.clear();
    encoded.reserve(SPILL_FRAME_HEADER_LEN.saturating_add(payload_len));
    encoded.extend_from_slice(&SPILL_MAGIC);
    encoded.extend_from_slice(&SPILL_VERSION.to_be_bytes());
    encoded.extend_from_slice(&payload_len_u64.to_be_bytes());

    encoded.extend_from_slice(&entry.rid.to_be_bytes());
    encoded.extend_from_slice(&entry.pos.to_be_bytes());
    encoded.extend_from_slice(&entry.source.blob_ordinal.to_be_bytes());
    encoded.extend_from_slice(&record_ordinal.to_be_bytes());
    encoded.extend_from_slice(&run_id.to_be_bytes());
    encoded.extend_from_slice(&packed_alleles_len.to_be_bytes());
    encoded.extend_from_slice(&record_payload_len.to_be_bytes());
    encoded.extend_from_slice(entry.packed_alleles.as_ref());
    encoded.extend_from_slice(entry.record_payload.as_slice());
    Ok(())
}

pub fn decode_frame_payload_len(frame_header: &[u8]) -> Result<usize> {
    if frame_header.len() != SPILL_FRAME_HEADER_LEN {
        return Err(crate::SpillSortError::message(format!(
            "spill frame header length {} does not match expected {}",
            frame_header.len(),
            SPILL_FRAME_HEADER_LEN
        )));
    }

    let observed_magic = &frame_header[MAGIC_OFFSET..VERSION_OFFSET];
    if observed_magic != SPILL_MAGIC.as_slice() {
        return Err(crate::SpillSortError::message(format!(
            "spill frame magic mismatch: expected {:?}, observed {:?}",
            SPILL_MAGIC, observed_magic
        )));
    }

    let version = {
        let mut bytes = [0_u8; std::mem::size_of::<u16>()];
        bytes.copy_from_slice(&frame_header[VERSION_OFFSET..PAYLOAD_LEN_OFFSET]);
        u16::from_be_bytes(bytes)
    };
    if version != SPILL_VERSION {
        return Err(crate::SpillSortError::message(format!(
            "spill frame version mismatch: expected {}, observed {}",
            SPILL_VERSION, version
        )));
    }

    let payload_len = {
        let mut bytes = [0_u8; std::mem::size_of::<u64>()];
        bytes.copy_from_slice(&frame_header[PAYLOAD_LEN_OFFSET..SPILL_FRAME_HEADER_LEN]);
        u64::from_be_bytes(bytes)
    };
    usize::try_from(payload_len).map_err(|_| {
        crate::SpillSortError::message(format!(
            "spill frame payload length {} exceeds platform usize::MAX",
            payload_len
        ))
    })
}

#[cfg(test)]
pub fn decode_entry(encoded: &[u8]) -> Result<SpillFormatEntry> {
    if encoded.len() < SPILL_FRAME_HEADER_LEN {
        return Err(crate::SpillSortError::message(format!(
            "spill frame truncated before header: got {} bytes, need at least {}",
            encoded.len(),
            SPILL_FRAME_HEADER_LEN
        )));
    }

    let (frame_header, payload) = encoded.split_at(SPILL_FRAME_HEADER_LEN);
    decode_entry_from_parts(frame_header, payload)
}

pub fn decode_entry_from_parts(frame_header: &[u8], payload: &[u8]) -> Result<SpillFormatEntry> {
    let payload_len = decode_frame_payload_len(frame_header)?;
    if payload.len() != payload_len {
        return Err(crate::SpillSortError::message(format!(
            "spill frame length mismatch: payload_len={} requires payload {} bytes, observed {}",
            payload_len,
            payload_len,
            payload.len()
        )));
    }
    if payload_len < FIXED_FIELDS_LEN {
        return Err(crate::SpillSortError::message(format!(
            "spill frame payload too short: {} bytes, need at least {}",
            payload_len, FIXED_FIELDS_LEN
        )));
    }

    let mut cursor = SpillDecodeCursor::new(payload);

    let rid = cursor.read_u32()?;
    let pos = cursor.read_i64()?;
    let blob_ordinal = cursor.read_u64()?;
    let record_ordinal_u64 = cursor.read_u64()?;
    let run_id = cursor.read_u64()?;
    let packed_alleles_len = cursor.read_u32()?;
    let record_payload_len = cursor.read_u32()?;

    let packed_alleles_len = usize::try_from(packed_alleles_len).map_err(|_| {
        crate::SpillSortError::message(
            "spill packed-alleles length exceeds platform usize::MAX".to_string(),
        )
    })?;
    let record_payload_len = usize::try_from(record_payload_len).map_err(|_| {
        crate::SpillSortError::message(
            "spill record payload length exceeds platform usize::MAX".to_string(),
        )
    })?;
    let total_variable_len = packed_alleles_len.saturating_add(record_payload_len);
    if cursor.remaining_len() != total_variable_len {
        return Err(crate::SpillSortError::message(format!(
            "spill frame variable payload mismatch: expected {} bytes, observed {} bytes",
            total_variable_len,
            cursor.remaining_len()
        )));
    }

    let packed_alleles = Arc::<[u8]>::from(cursor.read_bytes(packed_alleles_len)?);
    let record_payload = cursor.read_bytes(record_payload_len)?.to_vec();
    let record_ordinal = usize::try_from(record_ordinal_u64).map_err(|_| {
        crate::SpillSortError::message(format!(
            "spill source record ordinal {} exceeds platform usize::MAX",
            record_ordinal_u64
        ))
    })?;

    Ok(SpillFormatEntry {
        rid,
        pos,
        packed_alleles,
        source: SpillSourceOrder {
            blob_ordinal,
            record_ordinal,
        },
        run_id,
        record_payload,
    })
}

struct SpillDecodeCursor<'a> {
    payload: &'a [u8],
    offset: usize,
}

impl<'a> SpillDecodeCursor<'a> {
    fn new(payload: &'a [u8]) -> Self {
        Self { payload, offset: 0 }
    }

    fn read_u32(&mut self) -> Result<u32> {
        let bytes = self.read_bytes(std::mem::size_of::<u32>())?;
        let mut array = [0_u8; std::mem::size_of::<u32>()];
        array.copy_from_slice(bytes);
        Ok(u32::from_be_bytes(array))
    }

    fn read_u64(&mut self) -> Result<u64> {
        let bytes = self.read_bytes(std::mem::size_of::<u64>())?;
        let mut array = [0_u8; std::mem::size_of::<u64>()];
        array.copy_from_slice(bytes);
        Ok(u64::from_be_bytes(array))
    }

    fn read_i64(&mut self) -> Result<i64> {
        let bytes = self.read_bytes(std::mem::size_of::<i64>())?;
        let mut array = [0_u8; std::mem::size_of::<i64>()];
        array.copy_from_slice(bytes);
        Ok(i64::from_be_bytes(array))
    }

    fn read_bytes(&mut self, len: usize) -> Result<&'a [u8]> {
        let end = self.offset.saturating_add(len);
        if end > self.payload.len() {
            return Err(crate::SpillSortError::message(format!(
                "spill frame truncated while decoding: need {} more bytes, have {}",
                len,
                self.payload.len().saturating_sub(self.offset)
            )));
        }
        let bytes = &self.payload[self.offset..end];
        self.offset = end;
        Ok(bytes)
    }

    fn remaining_len(&self) -> usize {
        self.payload.len().saturating_sub(self.offset)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_entry() -> SpillFormatEntry {
        SpillFormatEntry {
            rid: 3,
            pos: 1234,
            packed_alleles: Arc::from(b"\x00\x00\x00\x01a\x00\x00\x00\x01t".as_slice()),
            source: SpillSourceOrder {
                blob_ordinal: 11,
                record_ordinal: 7,
            },
            run_id: 19,
            record_payload: b"packed-record".to_vec(),
        }
    }

    #[test]
    fn encode_decode_roundtrip_preserves_payload_and_tie_break_fields() {
        let input = sample_entry();

        let encoded =
            encode_entry(&input).expect("spill format should encode roundtrip sample entry");
        let decoded =
            decode_entry(encoded.as_slice()).expect("spill format should decode encoded entry");

        assert_eq!(decoded, input);
    }

    #[test]
    fn encode_entry_with_run_id_override_rewrites_only_run_id() {
        let input = sample_entry();
        let mut encoded = Vec::new();
        encode_entry_with_run_id_into(&input, 27, &mut encoded)
            .expect("spill format should encode with run_id override");
        let decoded =
            decode_entry(encoded.as_slice()).expect("spill format should decode encoded entry");

        let mut expected = input;
        expected.run_id = 27;
        assert_eq!(decoded, expected);
    }

    #[test]
    fn decode_rejects_corrupt_magic_prefix() {
        let mut corrupted = Vec::new();
        corrupted.extend_from_slice(b"BADS");
        corrupted.extend_from_slice(&SPILL_VERSION.to_be_bytes());

        let error = decode_entry(corrupted.as_slice())
            .expect_err("decoding should fail for non-HTSS spill magic");
        assert!(!error.to_string().is_empty());
    }

    #[test]
    fn decode_rejects_truncated_frame() {
        let mut truncated = Vec::new();
        truncated.extend_from_slice(&SPILL_MAGIC);
        truncated.extend_from_slice(&SPILL_VERSION.to_be_bytes());
        truncated.push(0);

        let error = decode_entry(truncated.as_slice())
            .expect_err("decoding should fail for truncated spill frame");
        assert!(!error.to_string().is_empty());
    }
}
