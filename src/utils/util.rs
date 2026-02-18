use crate::error::SvxResult;
use log;
use std::{
    fmt::{Binary, Display},
    path::Path,
    sync::Once,
};

pub type Result<T> = SvxResult<T>;

pub const MISSING_INTEGER: i32 = i32::MIN;
pub const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
pub const MISSING_FLOAT: f32 = f32::from_bits(0x7F80_0001);
pub const VECTOR_END_FLOAT: f32 = f32::from_bits(0x7F80_0002);

#[allow(unused)]
static INIT_LOG: Once = Once::new();

#[allow(unused)]
pub fn init_logger() {
    INIT_LOG.call_once(|| {
        env_logger::builder()
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .init();
    });
}

pub fn handle_error_and_exit(err: impl Display) -> ! {
    log::error!("{err}");
    std::process::exit(1);
}

pub fn try_exists(path: &Path) -> Result<()> {
    if !path.exists() {
        return Err(crate::svx_error!(
            "Path/File does not exist: {}",
            path.display()
        ));
    }
    Ok(())
}

pub fn log_warning<T>(err: impl Display, default: T) -> T {
    log::warn!("{err}");
    default
}

pub fn stable_hash64(bytes: &[u8]) -> u64 {
    const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

    let mut hash = FNV_OFFSET_BASIS;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

pub fn stable_hash(bytes: &[u8]) -> i32 {
    const MOD: i64 = 1_000_000_007;
    let mut res: i64 = 0;
    for &b in bytes {
        res = (res * 17 + i64::from(b)) % MOD;
    }
    res as i32
}

pub fn round_to_i64(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

pub fn to_info_i32(value: i64, label: &str) -> Result<i32> {
    i32::try_from(value).map_err(|_| {
        crate::svx_error!(
            "Cannot write {label}={value} to INFO as i32: value is outside supported range"
        )
    })
}

pub fn format_number_with_commas<T>(n: T) -> String
where
    T: Display + Binary,
{
    let s = n.to_string();
    let (sign, digits) = s.strip_prefix('-').map_or(("", s.as_str()), |d| ("-", d));

    if let 0..=3 = digits.len() {
        return s;
    }

    let mut result = String::with_capacity(digits.len() + (digits.len() - 1) / 3 + sign.len());
    for (digit_count, c) in digits.chars().rev().enumerate() {
        if digit_count > 0 && digit_count % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }

    result = result.chars().rev().collect();
    if !sign.is_empty() {
        result.insert_str(0, sign);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_number_with_commas_i8() {
        assert_eq!(format_number_with_commas(0i8), "0");
        assert_eq!(format_number_with_commas(100i8), "100");
        assert_eq!(format_number_with_commas(-100i8), "-100");
        assert_eq!(format_number_with_commas(i8::MAX), "127");
        assert_eq!(format_number_with_commas(i8::MIN), "-128");
    }

    #[test]
    fn test_format_number_with_commas_u8() {
        assert_eq!(format_number_with_commas(0u8), "0");
        assert_eq!(format_number_with_commas(100u8), "100");
        assert_eq!(format_number_with_commas(u8::MAX), "255");
    }

    #[test]
    fn test_format_number_with_commas_i16() {
        assert_eq!(format_number_with_commas(0i16), "0");
        assert_eq!(format_number_with_commas(1_000i16), "1,000");
        assert_eq!(format_number_with_commas(-1_000i16), "-1,000");
        assert_eq!(format_number_with_commas(i16::MAX), "32,767");
        assert_eq!(format_number_with_commas(i16::MIN), "-32,768");
    }

    #[test]
    fn test_format_number_with_commas_u16() {
        assert_eq!(format_number_with_commas(0u16), "0");
        assert_eq!(format_number_with_commas(1_000u16), "1,000");
        assert_eq!(format_number_with_commas(u16::MAX), "65,535");
    }

    #[test]
    fn test_format_number_with_commas_i32() {
        assert_eq!(format_number_with_commas(0i32), "0");
        assert_eq!(format_number_with_commas(10_000i32), "10,000");
        assert_eq!(format_number_with_commas(-10_000i32), "-10,000");
        assert_eq!(format_number_with_commas(i32::MAX), "2,147,483,647");
        assert_eq!(format_number_with_commas(i32::MIN), "-2,147,483,648");
    }

    #[test]
    fn test_format_number_with_commas_u32() {
        assert_eq!(format_number_with_commas(0u32), "0");
        assert_eq!(format_number_with_commas(10_000u32), "10,000");
        assert_eq!(format_number_with_commas(u32::MAX), "4,294,967,295");
    }

    #[test]
    fn test_format_number_with_commas_i64() {
        assert_eq!(format_number_with_commas(0i64), "0");
        assert_eq!(format_number_with_commas(1_000_000i64), "1,000,000");
        assert_eq!(format_number_with_commas(-1_000_000i64), "-1,000,000");
        assert_eq!(
            format_number_with_commas(i64::MAX),
            "9,223,372,036,854,775,807"
        );
        assert_eq!(
            format_number_with_commas(i64::MIN),
            "-9,223,372,036,854,775,808"
        );
    }

    #[test]
    fn test_format_number_with_commas_u64() {
        assert_eq!(format_number_with_commas(0u64), "0");
        assert_eq!(format_number_with_commas(1_000_000u64), "1,000,000");
        assert_eq!(
            format_number_with_commas(u64::MAX),
            "18,446,744,073,709,551,615"
        );
    }

    #[test]
    fn test_format_number_with_commas_usize() {
        assert_eq!(format_number_with_commas(0usize), "0");
        assert_eq!(
            format_number_with_commas(1_234_567_890usize),
            "1,234,567,890"
        );
        assert_eq!(
            format_number_with_commas(usize::MAX),
            "18,446,744,073,709,551,615"
        );
    }

    #[test]
    fn test_format_number_with_commas_isize() {
        assert_eq!(format_number_with_commas(0isize), "0");
        assert_eq!(
            format_number_with_commas(-1_234_567_890isize),
            "-1,234,567,890"
        );
        assert_eq!(
            format_number_with_commas(isize::MAX),
            "9,223,372,036,854,775,807"
        );
        assert_eq!(
            format_number_with_commas(isize::MIN),
            "-9,223,372,036,854,775,808"
        );
    }
}
