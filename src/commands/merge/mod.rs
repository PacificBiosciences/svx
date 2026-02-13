use crate::{
    cli::MergeArgs,
    core::{svtype::SvType, variant::VariantInternal, variant_block::VariantBlockResult},
    io::{
        bed_reader::BedMap,
        merge_reader::load_contig,
        merge_writer::{create_output_header, write_variants},
        positions_reader::positions_to_interval_trees,
        vcf_reader::VcfReaders,
        vcf_writer::VcfWriter,
    },
    utils::util::Result,
};
use crossbeam_channel::{Receiver, Sender, bounded, unbounded};
use rayon::{ThreadPoolBuilder, prelude::*};
use std::{
    collections::HashSet,
    io::{self, IsTerminal},
    path::Path,
    sync::Arc,
    thread,
};

mod bnd;
mod dump;
mod progress;
mod shutdown;
mod types;
mod worker;

pub use types::VariantBlob;

use bnd::{
    BndMateSelectionScope, bucket_bnd_variants, pair_bnd_breakends, pair_bnd_breakends_with_scope,
};
use dump::DumpWriter;
use progress::{
    PipelineQueueMetrics, ProgressEvent, compute_progress_enabled,
    resolve_pipeline_queue_capacities, run_progress_ui, send_progress_event,
};
use shutdown::finalize_merge_threads;
use worker::{count_records_to_write, count_variants_in_groups, process_blob};

#[cfg(test)]
mod tests;

fn get_contig_order(vcf_readers: &VcfReaders, args: &MergeArgs) -> Result<Vec<String>> {
    let mut contig_order = vcf_readers.get_contig_order()?;
    if let Some(ref user_contigs) = args.contigs {
        let user_contig_set: HashSet<&String> = user_contigs.iter().collect();
        let original_contig_set: HashSet<&String> = contig_order.iter().collect();
        if !user_contig_set.is_subset(&original_contig_set) {
            let missing_contigs: Vec<&&String> =
                user_contig_set.difference(&original_contig_set).collect();
            return Err(crate::svx_error!(
                "The following user-specified contigs do not exist in the VCF files: {:?}",
                missing_contigs
            ));
        }
        contig_order.retain(|contig| user_contig_set.contains(contig));
    }
    // Filter based on exact positions
    if let Some(ref target_positions) = args.target_positions {
        let target_contigs: HashSet<&String> = target_positions
            .iter()
            .map(|pos_tuple| &pos_tuple.contig)
            .collect();
        contig_order.retain(|contig| target_contigs.contains(contig));
    }
    Ok(contig_order)
}

pub fn merge(args: MergeArgs) -> Result<()> {
    let vcf_paths = args
        .process_vcf_paths()
        .map_err(|error| crate::svx_error!("{error}"))?;
    let selected_svtypes = args.selected_svtypes();
    let vcf_readers = VcfReaders::new(vcf_paths, args.num_io_threads)?;
    if vcf_readers.readers.len() == 1 && !args.force_single {
        return Err(crate::svx_error!(
            "Expected two or more files to merge, got only one. Use --force-single to proceed anyway"
        ));
    }

    let contig_order = get_contig_order(&vcf_readers, &args)?;
    let progress_enabled = compute_progress_enabled(
        args.progress,
        args.no_progress,
        io::stderr().is_terminal(),
        log::max_level(),
    );

    // TODO: Move into reader thread
    let bed_map = if let Some(ref bed_path) = args.tr_bed_path {
        Some(BedMap::new(bed_path)?)
    } else {
        None
    };

    let (out_header, sample_mapping) = create_output_header(&vcf_readers, &args)?;
    let writer_io_tpool = vcf_readers.io_tpool.clone();
    let writer = VcfWriter::new(
        &out_header,
        &args.output_type,
        args.output.as_ref(),
        writer_io_tpool,
    )?;

    let dump_writer = if let Some(dump_path) = args.merge_args.dump_path.as_ref() {
        if matches!(
            args.output.as_ref(),
            Some(output_path) if Path::new(output_path) == dump_path
        ) {
            return Err(crate::svx_error!(
                "The dump path and output path must be different: {}",
                dump_path.display()
            ));
        }

        Some(Arc::new(DumpWriter::from_path(dump_path)?))
    } else {
        None
    };

    if args.print_header {
        return Ok(());
    }

    // TODO: This is the n of VCFs not the samples, ok for now since we enforce single sample VCFs
    let n_samples = vcf_readers.n;
    let mut readers = vcf_readers.readers;

    let (blob_queue_capacity, result_queue_capacity) = resolve_pipeline_queue_capacities(
        args.num_threads,
        args.blob_queue_capacity,
        args.result_queue_capacity,
    );
    log::debug!(
        "Pipeline queue capacities: blob={} result={}",
        blob_queue_capacity,
        result_queue_capacity
    );
    let queue_metrics = Arc::new(PipelineQueueMetrics::default());

    let (blob_sender, blob_receiver): (Sender<VariantBlob>, Receiver<VariantBlob>) =
        bounded(blob_queue_capacity);
    let (result_sender, result_receiver): (
        Sender<VariantBlockResult>,
        Receiver<VariantBlockResult>,
    ) = bounded(result_queue_capacity);

    let total_contigs = contig_order.len() as u64;
    let mut progress_thread = None;
    let progress_sender = if progress_enabled {
        let (sender, receiver): (Sender<ProgressEvent>, Receiver<ProgressEvent>) = unbounded();
        let queue_metrics_progress = Arc::clone(&queue_metrics);
        progress_thread = Some(thread::spawn(move || {
            run_progress_ui(receiver, total_contigs, Some(queue_metrics_progress));
        }));
        Some(sender)
    } else {
        None
    };

    let progress_sender_reader = progress_sender.clone();
    let queue_metrics_reader = Arc::clone(&queue_metrics);
    let dump_writer_reader = dump_writer.clone();
    let reader_thread = thread::spawn(move || -> Result<()> {
        let tx = blob_sender;

        let pos_map = if let Some(ref target_positions) = args.target_positions {
            Some(positions_to_interval_trees(target_positions)?)
        } else {
            None
        };
        let bnd_mate_selection_scope = if args.contigs.is_some() || args.target_positions.is_some()
        {
            Some(BndMateSelectionScope::new(
                &contig_order,
                args.target_positions.as_deref(),
            ))
        } else {
            None
        };

        let mut all_bnd_variants: Vec<VariantInternal> = Vec::new();

        log::debug!("Reader thread started.");
        for contig_name in contig_order {
            match load_contig(
                &contig_name,
                &args.merge_args,
                &selected_svtypes,
                &bed_map,
                &pos_map,
                &mut readers,
            ) {
                Ok(mut variants_by_type_for_contig) => {
                    if variants_by_type_for_contig.len() != SvType::BUCKET_COUNT {
                        return Err(crate::svx_error!(
                            "Unexpected SV type bucket count {} while processing contig {}",
                            variants_by_type_for_contig.len(),
                            contig_name
                        ));
                    }

                    let bnd_variants = {
                        let bnd_bucket = variants_by_type_for_contig
                            .get_mut(SvType::BND.bucket_index())
                            .ok_or_else(|| {
                                crate::svx_error!(
                                    "Missing BND bucket while processing contig {}",
                                    contig_name
                                )
                            })?;
                        std::mem::take(bnd_bucket)
                    };

                    for variant_type in SvType::NON_BND_BUCKET_TYPES {
                        let variants = {
                            let bucket = variants_by_type_for_contig
                                .get_mut(variant_type.bucket_index())
                                .ok_or_else(|| {
                                    crate::svx_error!(
                                        "Missing {} bucket while processing contig {}",
                                        variant_type,
                                        contig_name
                                    )
                                })?;
                            std::mem::take(bucket)
                        };
                        if variants.is_empty() {
                            continue;
                        }
                        let blob = VariantBlob {
                            variants,
                            contig: contig_name.clone(),
                            variant_type,
                            args: args.merge_args.clone(),
                            dump_writer: dump_writer_reader.clone(),
                        };
                        log::debug!(
                            "Reader: Sending blob for Contig={}, Type={}, Variants={}",
                            blob.contig,
                            blob.variant_type,
                            blob.variants.len()
                        );
                        let variant_count = blob.variants.len() as u64;
                        queue_metrics_reader.blob.increment();
                        if tx.send(blob).is_err() {
                            queue_metrics_reader.blob.decrement();
                            log::error!(
                                "Failed to send variant blob for contig {} / type {}. Receiver closed.",
                                contig_name,
                                variant_type
                            );
                            return Err(crate::svx_error!(
                                "Channel receiver closed unexpectedly in reader thread"
                            ));
                        }
                        send_progress_event(
                            progress_sender_reader.as_ref(),
                            ProgressEvent::VariantsQueued {
                                variants: variant_count,
                            },
                        );
                    }

                    if !bnd_variants.is_empty() {
                        log::debug!(
                            "Reader: Collected {} BND variants from contig {}",
                            bnd_variants.len(),
                            contig_name
                        );
                        all_bnd_variants.extend(bnd_variants);
                    }
                    send_progress_event(
                        progress_sender_reader.as_ref(),
                        ProgressEvent::ContigCompleted {
                            contig: contig_name.clone(),
                        },
                    );
                }
                Err(e) => {
                    log::error!(
                        "Error processing contig {} in reader thread: {}",
                        contig_name,
                        e
                    );
                    return Err(e);
                }
            }
        }
        log::debug!(
            "Reader thread finished processing all contigs. Total BNDs collected: {}",
            all_bnd_variants.len()
        );

        if !all_bnd_variants.is_empty() {
            let events = if let Some(scope) = bnd_mate_selection_scope.as_ref() {
                pair_bnd_breakends_with_scope(all_bnd_variants, &args.merge_args, Some(scope))?
            } else {
                pair_bnd_breakends(all_bnd_variants, &args.merge_args)?
            };
            let by_graph = bucket_bnd_variants(events)?;

            for (graph_id, variants) in by_graph {
                if variants.is_empty() {
                    continue;
                }
                let blob = VariantBlob {
                    variants,
                    contig: graph_id.clone(),
                    variant_type: SvType::BND,
                    args: args.merge_args.clone(),
                    dump_writer: dump_writer_reader.clone(),
                };
                log::debug!(
                    "Reader: Sending BND blob for Graph={}, Variants={}",
                    blob.contig,
                    blob.variants.len()
                );
                let variant_count = blob.variants.len() as u64;
                queue_metrics_reader.blob.increment();
                if tx.send(blob).is_err() {
                    queue_metrics_reader.blob.decrement();
                    log::error!(
                        "Failed to send BND variant blob for graph {}. Receiver closed.",
                        graph_id
                    );
                    return Err(crate::svx_error!(
                        "Channel receiver closed unexpectedly in reader thread"
                    ));
                }
                send_progress_event(
                    progress_sender_reader.as_ref(),
                    ProgressEvent::VariantsQueued {
                        variants: variant_count,
                    },
                );
            }
        }

        send_progress_event(progress_sender_reader.as_ref(), ProgressEvent::ReaderDone);
        Ok(())
    });

    let progress_sender_writer = progress_sender.clone();
    let queue_metrics_writer = Arc::clone(&queue_metrics);
    let writer_thread = thread::spawn(move || -> Result<()> {
        let mut writer_instance = writer;
        log::debug!("Writer thread started.");

        for payload in result_receiver {
            queue_metrics_writer.result.decrement();
            log::debug!(
                "Writer: Processing result: contig={}, type={}, groups={}",
                payload.contig,
                payload.variant_type,
                payload.groups.len()
            );

            log::debug!(
                "Writer: Writing variants for {}/{}",
                payload.contig,
                payload.variant_type
            );
            let written_variants = count_variants_in_groups(&payload) as u64;
            let written_records = count_records_to_write(&payload) as u64;
            write_variants(payload, &sample_mapping, n_samples, &mut writer_instance)?;
            send_progress_event(
                progress_sender_writer.as_ref(),
                ProgressEvent::WriterAdvanced {
                    variants: written_variants,
                    records: written_records,
                },
            );
        }

        log::debug!("Writer thread finished.");
        Ok(())
    });

    log::debug!(
        "Initializing merge thread pool with {} threads...",
        args.num_threads
    );
    let pool = ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .thread_name(|i| format!("svx-merge-{i}"))
        .build()
        .map_err(|e| crate::svx_error!("Failed to initialize merge thread pool: {e}"))?;

    let progress_sender_workers = progress_sender.clone();
    let queue_metrics_workers = Arc::clone(&queue_metrics);
    let worker_result: Result<()> = pool.install(|| {
        blob_receiver.into_iter().par_bridge().try_for_each_with(
            (
                result_sender,
                progress_sender_workers,
                queue_metrics_workers,
            ),
            |state, variant_blob| {
                process_blob(variant_blob, &state.0, state.1.as_ref(), state.2.as_ref())
            },
        )
    });

    let shutdown_result = finalize_merge_threads(
        reader_thread,
        writer_thread,
        progress_sender,
        progress_thread,
    );

    let queue_snapshot = queue_metrics.snapshot();
    log::info!(
        "Pipeline queue peak depth: blob={} result={}",
        queue_snapshot.blob.peak,
        queue_snapshot.result.peak
    );

    worker_result?;
    shutdown_result
}
