use crate::Result;
use rust_htslib::{bcf, htslib};

const PACKED_RECORD_PAYLOAD_MAGIC: [u8; 4] = *b"HTSP";
const PACKED_RECORD_PAYLOAD_VERSION: u16 = 1;
const PACKED_RECORD_PAYLOAD_HEADER_LEN: usize =
    PACKED_RECORD_PAYLOAD_MAGIC.len() + std::mem::size_of::<u16>();
#[cfg(test)]
const PACKED_RECORD_PAYLOAD_LOCUS_FIELDS_LEN: usize =
    std::mem::size_of::<u32>() + std::mem::size_of::<i64>();
const PACKED_RECORD_PAYLOAD_FIXED_FIELDS_LEN: usize = std::mem::size_of::<u32>() // rid
    + std::mem::size_of::<i64>() // pos
    + std::mem::size_of::<i64>() // rlen
    + std::mem::size_of::<u32>() // qual bits
    + std::mem::size_of::<u32>() // n_info
    + std::mem::size_of::<u32>() // n_allele
    + std::mem::size_of::<u32>() // n_fmt
    + std::mem::size_of::<u32>() // n_sample
    + std::mem::size_of::<u32>() // shared bytes len
    + std::mem::size_of::<u32>(); // indiv bytes len

pub fn encode_record_payload(record: &bcf::Record) -> Result<Vec<u8>> {
    let synced = unsafe { htslib::bcf_dup(record.inner) };
    if synced.is_null() {
        return Err(crate::SpillSortError::message(
            "Failed to encode spill record payload: bcf_dup returned null".to_string(),
        ));
    }

    let result = {
        let inner = unsafe { &*synced };
        let rid = u32::try_from(inner.rid).map_err(|_| {
            crate::SpillSortError::message(format!(
                "Failed to encode spill record payload: record RID {} is invalid",
                inner.rid
            ))
        })?;
        let shared = kstring_payload_slice(&inner.shared, "shared")?;
        let indiv = kstring_payload_slice(&inner.indiv, "indiv")?;
        let shared_len = u32::try_from(shared.len()).map_err(|_| {
            crate::SpillSortError::message(format!(
                "Failed to encode spill record payload: shared length {} exceeds u32::MAX",
                shared.len()
            ))
        })?;
        let indiv_len = u32::try_from(indiv.len()).map_err(|_| {
            crate::SpillSortError::message(format!(
                "Failed to encode spill record payload: indiv length {} exceeds u32::MAX",
                indiv.len()
            ))
        })?;

        let payload_len = PACKED_RECORD_PAYLOAD_HEADER_LEN
            .saturating_add(PACKED_RECORD_PAYLOAD_FIXED_FIELDS_LEN)
            .saturating_add(shared.len())
            .saturating_add(indiv.len());
        let mut payload = Vec::with_capacity(payload_len);
        payload.extend_from_slice(&PACKED_RECORD_PAYLOAD_MAGIC);
        payload.extend_from_slice(&PACKED_RECORD_PAYLOAD_VERSION.to_be_bytes());
        payload.extend_from_slice(&rid.to_be_bytes());
        payload.extend_from_slice(&inner.pos.to_be_bytes());
        payload.extend_from_slice(&inner.rlen.to_be_bytes());
        payload.extend_from_slice(&inner.qual.to_bits().to_be_bytes());
        payload.extend_from_slice(&inner.n_info().to_be_bytes());
        payload.extend_from_slice(&inner.n_allele().to_be_bytes());
        payload.extend_from_slice(&inner.n_fmt().to_be_bytes());
        payload.extend_from_slice(&inner.n_sample().to_be_bytes());
        payload.extend_from_slice(&shared_len.to_be_bytes());
        payload.extend_from_slice(&indiv_len.to_be_bytes());
        payload.extend_from_slice(shared);
        payload.extend_from_slice(indiv);
        Ok(payload)
    };

    unsafe {
        htslib::bcf_destroy(synced);
    }
    result
}

pub fn duplicate_header_view(header: &bcf::Header) -> Result<bcf::header::HeaderView> {
    let duplicated_header = unsafe { htslib::bcf_hdr_dup(header.inner) };
    if duplicated_header.is_null() {
        return Err(crate::SpillSortError::message(
            "Failed to duplicate BCF header for spill payload decoding".to_string(),
        ));
    }
    Ok(bcf::header::HeaderView::new(duplicated_header))
}

#[cfg(test)]
pub fn parse_payload_prefix(record_payload: &[u8]) -> Result<(u32, i64)> {
    let (_, rid, pos) = parse_payload_prefix_with_min_len(
        record_payload,
        PACKED_RECORD_PAYLOAD_HEADER_LEN.saturating_add(PACKED_RECORD_PAYLOAD_LOCUS_FIELDS_LEN),
    )?;
    Ok((rid, pos))
}

#[cfg(test)]
pub fn decode_record_payload(
    record_payload: &[u8],
    decode_header: &bcf::header::HeaderView,
) -> Result<bcf::Record> {
    let (rid, pos) = parse_payload_prefix(record_payload)?;
    let mut record = decode_header.empty_record();
    decode_record_payload_into(record_payload, rid, pos, &mut record)?;
    Ok(record)
}

pub fn decode_record_payload_into(
    record_payload: &[u8],
    expected_rid: u32,
    expected_pos: i64,
    record: &mut bcf::Record,
) -> Result<()> {
    let (mut cursor, rid, pos) = parse_payload_prefix_with_min_len(
        record_payload,
        PACKED_RECORD_PAYLOAD_HEADER_LEN.saturating_add(PACKED_RECORD_PAYLOAD_FIXED_FIELDS_LEN),
    )?;
    if rid != expected_rid || pos != expected_pos {
        return Err(crate::SpillSortError::message(format!(
            "corrupt packed spill record payload: key mismatch payload_rid={} payload_pos={} expected_rid={} expected_pos={}",
            rid, pos, expected_rid, expected_pos
        )));
    }

    let rlen = cursor.read_i64()?;
    let qual_bits = cursor.read_u32()?;
    let n_info = cursor.read_u32()?;
    let n_allele = cursor.read_u32()?;
    let n_fmt = cursor.read_u32()?;
    let n_sample = cursor.read_u32()?;
    let shared_len = cursor.read_len_u32("shared")?;
    let indiv_len = cursor.read_len_u32("indiv")?;
    let shared = cursor.read_bytes(shared_len)?;
    let indiv = cursor.read_bytes(indiv_len)?;
    if cursor.remaining_len() != 0 {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: trailing {} bytes after decode",
            cursor.remaining_len()
        )));
    }

    validate_bcf_counter(n_info, "n_info", u32::from(u16::MAX))?;
    validate_bcf_counter(n_allele, "n_allele", u32::from(u16::MAX))?;
    validate_bcf_counter(n_fmt, "n_fmt", u32::from(u8::MAX))?;
    validate_bcf_counter(n_sample, "n_sample", 0x00ff_ffff)?;
    if rid > i32::MAX as u32 {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: rid {} exceeds i32::MAX",
            rid
        )));
    }

    record.clear();
    record.set_rid(Some(rid));
    record.set_pos(pos);
    let inner = record.inner_mut();
    inner.rlen = rlen;
    inner.qual = f32::from_bits(qual_bits);
    inner.set_n_info(n_info);
    inner.set_n_allele(n_allele);
    inner.set_n_fmt(n_fmt);
    inner.set_n_sample(n_sample);
    replace_kstring_bytes(&mut inner.shared, shared, "shared")?;
    replace_kstring_bytes(&mut inner.indiv, indiv, "indiv")?;
    inner.unpacked = 0;
    inner.errcode = 0;

    let unpack_result = unsafe { htslib::bcf_unpack(record.inner, htslib::BCF_UN_ALL as i32) };
    if unpack_result != 0 {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: bcf_unpack returned {}",
            unpack_result
        )));
    }
    if record.inner().errcode != 0 {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: decoded record errcode={}",
            record.inner().errcode
        )));
    }
    Ok(())
}

fn parse_payload_prefix_with_min_len(
    record_payload: &[u8],
    min_payload_len: usize,
) -> Result<(PackedRecordDecodeCursor<'_>, u32, i64)> {
    if record_payload.len() < min_payload_len {
        return Err(crate::SpillSortError::message(format!(
            "truncated packed spill record payload: payload too short ({} bytes)",
            record_payload.len()
        )));
    }

    let mut cursor = PackedRecordDecodeCursor::new(record_payload);
    let observed_magic = cursor.read_bytes(PACKED_RECORD_PAYLOAD_MAGIC.len())?;
    if observed_magic != PACKED_RECORD_PAYLOAD_MAGIC.as_slice() {
        return Err(crate::SpillSortError::message(format!(
            "corrupt packed spill record payload: magic mismatch expected {:?}, observed {:?}",
            PACKED_RECORD_PAYLOAD_MAGIC, observed_magic
        )));
    }

    let payload_version = cursor.read_u16()?;
    if payload_version != PACKED_RECORD_PAYLOAD_VERSION {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: version mismatch expected {}, observed {}",
            PACKED_RECORD_PAYLOAD_VERSION, payload_version
        )));
    }

    let rid = cursor.read_u32()?;
    let pos = cursor.read_i64()?;
    Ok((cursor, rid, pos))
}

fn kstring_payload_slice<'a>(payload: &'a htslib::kstring_t, field_name: &str) -> Result<&'a [u8]> {
    if payload.l == 0 {
        return Ok(&[]);
    }
    if payload.s.is_null() {
        return Err(crate::SpillSortError::message(format!(
            "Failed to encode spill record payload: {} data pointer is null with non-zero length {}",
            field_name, payload.l
        )));
    }

    Ok(unsafe { std::slice::from_raw_parts(payload.s.cast::<u8>(), payload.l) })
}

fn replace_kstring_bytes(
    target: &mut htslib::kstring_t,
    source: &[u8],
    field_name: &str,
) -> Result<()> {
    if source.len() > target.m {
        let requested = u64::try_from(source.len()).map_err(|_| {
            crate::SpillSortError::message(format!(
                "invalid packed spill record payload: {} length {} exceeds u64::MAX",
                field_name,
                source.len()
            ))
        })?;
        let resized = unsafe { htslib::realloc(target.s.cast(), requested) };
        if resized.is_null() {
            return Err(crate::SpillSortError::message(format!(
                "invalid packed spill record payload: realloc failed for {} data ({} bytes)",
                field_name,
                source.len()
            )));
        }
        target.s = resized.cast();
        target.m = source.len();
    }
    target.l = source.len();
    if source.is_empty() {
        return Ok(());
    }
    if target.s.is_null() {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: {} destination pointer is null for non-empty data",
            field_name
        )));
    }

    unsafe {
        std::ptr::copy_nonoverlapping(source.as_ptr(), target.s.cast::<u8>(), source.len());
    }
    Ok(())
}

fn validate_bcf_counter(value: u32, field_name: &str, max: u32) -> Result<()> {
    if value > max {
        return Err(crate::SpillSortError::message(format!(
            "invalid packed spill record payload: {}={} exceeds {}",
            field_name, value, max
        )));
    }
    Ok(())
}

struct PackedRecordDecodeCursor<'a> {
    payload: &'a [u8],
    offset: usize,
}

impl<'a> PackedRecordDecodeCursor<'a> {
    fn new(payload: &'a [u8]) -> Self {
        Self { payload, offset: 0 }
    }

    fn read_u16(&mut self) -> Result<u16> {
        let bytes = self.read_bytes(std::mem::size_of::<u16>())?;
        let mut array = [0_u8; std::mem::size_of::<u16>()];
        array.copy_from_slice(bytes);
        Ok(u16::from_be_bytes(array))
    }

    fn read_u32(&mut self) -> Result<u32> {
        let bytes = self.read_bytes(std::mem::size_of::<u32>())?;
        let mut array = [0_u8; std::mem::size_of::<u32>()];
        array.copy_from_slice(bytes);
        Ok(u32::from_be_bytes(array))
    }

    fn read_i64(&mut self) -> Result<i64> {
        let bytes = self.read_bytes(std::mem::size_of::<i64>())?;
        let mut array = [0_u8; std::mem::size_of::<i64>()];
        array.copy_from_slice(bytes);
        Ok(i64::from_be_bytes(array))
    }

    fn read_len_u32(&mut self, field_name: &str) -> Result<usize> {
        let len = self.read_u32()?;
        usize::try_from(len).map_err(|_| {
            crate::SpillSortError::message(format!(
                "invalid packed spill record payload: {} length {} exceeds usize::MAX",
                field_name, len
            ))
        })
    }

    fn read_bytes(&mut self, len: usize) -> Result<&'a [u8]> {
        let end = self.offset.saturating_add(len);
        if end > self.payload.len() {
            return Err(crate::SpillSortError::message(format!(
                "truncated packed spill record payload: need {} bytes, have {}",
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
