use crate::{
    cli::MergeArgsInner,
    core::{svtype::SvType, variant::VariantInternal},
};

use super::dump::SharedDumpWriter;

pub struct VariantBlob {
    pub variants: Vec<VariantInternal>,
    pub contig: String,
    pub variant_type: SvType,
    pub args: MergeArgsInner,
    pub(crate) dump_writer: Option<SharedDumpWriter>,
}
