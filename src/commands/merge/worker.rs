use crate::core::variant_block::{VariantBlock, VariantBlockResult};
use crossbeam_channel::Sender;
use std::thread;

use super::{
    progress::{PipelineQueueMetrics, ProgressEvent, send_progress_event},
    types::VariantBlob,
};

pub fn process_blob(
    variant_blob: VariantBlob,
    sender: &Sender<VariantBlockResult>,
    progress_sender: Option<&Sender<ProgressEvent>>,
    queue_metrics: &PipelineQueueMetrics,
) -> crate::utils::util::Result<()> {
    queue_metrics.blob.decrement();
    let current_thread = thread::current();
    let worker_name = current_thread.name().unwrap_or("unnamed");
    log::debug!(
        "Worker [{worker_name}]: Processing blob: contig={}, type={}, variants={}",
        variant_blob.contig,
        variant_blob.variant_type,
        variant_blob.variants.len()
    );

    if variant_blob.variants.is_empty() {
        log::debug!(
            "Worker [{worker_name}]: Received empty variant blob for {} - {}, skipping.",
            variant_blob.contig,
            variant_blob.variant_type
        );
        return Ok(());
    }

    let dump_writer = variant_blob.dump_writer.clone();
    if let Some(dump_writer) = dump_writer.as_ref() {
        dump_writer.dump_kd_coordinates(
            &variant_blob.contig,
            variant_blob.variant_type,
            &variant_blob.variants,
        )?;
    }

    let mut variant_block = VariantBlock::new(variant_blob);
    log::debug!(
        "Worker [{worker_name}]: Starting merge for {}/{}",
        variant_block.contig,
        variant_block.variant_type
    );
    variant_block.merge_block();
    log::debug!(
        "Worker [{worker_name}]: Merge complete for {}/{}, getting groups",
        variant_block.contig,
        variant_block.variant_type
    );
    let variant_block_result = variant_block.get_groups();

    if let Some(dump_writer) = dump_writer.as_ref() {
        dump_writer.dump_group_stats(&variant_block_result)?;
    }

    log::debug!(
        "Worker [{worker_name}]: Sending result to writer: contig={}, type={}, groups={}",
        variant_block_result.contig,
        variant_block_result.variant_type,
        variant_block_result.groups.len()
    );
    send_progress_event(
        progress_sender,
        ProgressEvent::VariantsMerged {
            variants: variant_block_result.n as u64,
        },
    );
    queue_metrics.result.increment();
    if let Err(error) = sender.send(variant_block_result) {
        queue_metrics.result.decrement();
        log::error!("Worker [{worker_name}]: Failed to send result to writer thread: {error}");
        return Err(crate::svx_error!(
            "Failed to send result to writer thread: {error}"
        ));
    }
    Ok(())
}
