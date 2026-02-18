use std::fmt;

#[derive(Clone, Copy)]
pub struct VariantLabel<'a> {
    pub idx: usize,
    pub id: &'a str,
}

impl<'a> VariantLabel<'a> {
    pub fn new(idx: usize, id: &'a str) -> Self {
        Self { idx, id }
    }
}

impl fmt::Display for VariantLabel<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.idx, self.id)
    }
}
