use std::fmt;

#[derive(Clone, Copy)]
pub(super) struct VariantLabel<'a> {
    pub(super) idx: usize,
    pub(super) id: &'a str,
}

impl<'a> VariantLabel<'a> {
    pub(super) fn new(idx: usize, id: &'a str) -> Self {
        Self { idx, id }
    }
}

impl fmt::Display for VariantLabel<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.idx, self.id)
    }
}
