use crate::utils::util::Result;
use anyhow::anyhow;

#[derive(Debug, Clone, PartialEq, Copy)]
pub enum SvType {
    INSERTION,
    DELETION,
    INVERSION,
    DUPLICATION,
    BND,
}

impl SvType {
    pub fn from_u8(bytes: &[u8]) -> Result<Self> {
        match bytes {
            b"INS" => Ok(SvType::INSERTION),
            b"DEL" => Ok(SvType::DELETION),
            b"INV" => Ok(SvType::INVERSION),
            b"DUP" => Ok(SvType::DUPLICATION),
            b"BND" => Ok(SvType::BND),
            _ => Err(anyhow!(
                "Invalid SVTYPE: {:?}",
                String::from_utf8_lossy(bytes)
            )),
        }
    }
}

impl std::str::FromStr for SvType {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        Self::from_u8(s.as_bytes())
    }
}

impl std::fmt::Display for SvType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SvType::INSERTION => write!(f, "INS"),
            SvType::DELETION => write!(f, "DEL"),
            SvType::INVERSION => write!(f, "INV"),
            SvType::DUPLICATION => write!(f, "DUP"),
            SvType::BND => write!(f, "BND"),
        }
    }
}
