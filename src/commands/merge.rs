use crate::{
    cli::{MergeArgs, MergeArgsInner},
    core::{
        svtype::SvType,
        variant_block::{VariantBlock, VariantBlockResult},
        variant_internal::VariantInternal,
    },
    init_config,
    io::{
        bed_reader::BedMap,
        merge_reader::load_contig,
        merge_writer::{create_output_header, write_variants},
        positions_reader::positions_to_interval_trees,
        vcf_reader::VcfReaders,
        vcf_writer::VcfWriter,
    },
    utils::util::Result,
    SvxConfig,
};
use anyhow::anyhow;
use crossbeam_channel::{unbounded, Receiver, Sender};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::{
    collections::{HashMap, HashSet},
    thread,
};

pub struct VariantBlob {
    pub variants: Vec<VariantInternal>,
    pub contig: String,
    pub variant_type: SvType,
    pub args: MergeArgsInner,
}

pub fn print_group_stats(variant_block_result: &VariantBlockResult) {
    let mut category_counts: HashMap<usize, usize> = HashMap::new();
    for list in &variant_block_result.groups {
        let group_size = list.len();
        if group_size == 0 {
            continue;
        }
        *category_counts.entry(group_size).or_insert(0) += 1;
        println!("Variant:");
        for v in list {
            println!(
                "id: {}, sample: {}, start: {}, end: {}, svlen: {}",
                v.id, v.sample_id, v.start, v.end, v.svlen
            );
        }
    }
    println!("\nGroups:");
    let mut sorted_category_counts: Vec<_> = category_counts.iter().collect();
    sorted_category_counts.sort_by_key(|&(group_size, _)| group_size);
    for (group_size, count) in sorted_category_counts {
        println!("  Groups of size {:>3}: {}", group_size, count);
    }
}

fn get_contig_order(vcf_readers: &VcfReaders, args: &MergeArgs) -> Result<Vec<String>> {
    let mut contig_order = vcf_readers.get_contig_order()?;
    if let Some(ref user_contigs) = args.contigs {
        let user_contig_set: HashSet<&String> = user_contigs.iter().collect();
        let original_contig_set: HashSet<&String> = contig_order.iter().collect();
        if !user_contig_set.is_subset(&original_contig_set) {
            let missing_contigs: Vec<&&String> =
                user_contig_set.difference(&original_contig_set).collect();
            return Err(anyhow!(
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
    // TODO: Get rid of globals
    init_config(SvxConfig {
        kd_tree_norm: args.merge_args.kd_tree_norm,
        dump: args.merge_args.dump,
    });

    let vcf_paths = args.process_vcf_paths()?;
    // TODO: This happens twice currently, which is unneeded
    let vcf_readers = VcfReaders::new(vcf_paths.clone())?;
    if vcf_readers.readers.len() == 1 && !args.force_single {
        return Err(anyhow!("Expected two or more files to merge, got only one. Use --force-single to proceed anyway"));
    }

    let contig_order = get_contig_order(&vcf_readers, &args)?;

    // TODO: Move into reader thread
    let bed_map = if let Some(ref bed_path) = args.tr_bed_path {
        Some(BedMap::new(bed_path)?)
    } else {
        None
    };

    let (out_header, sample_mapping) = create_output_header(&vcf_readers, &args)?;
    let writer = VcfWriter::new(&out_header, &args.output_type, args.output.as_ref())?;

    if args.print_header {
        return Ok(());
    }

    // TODO: This is the n of VCFs not the samples, ok for now since we enforce single sample VCFs
    let n_samples = vcf_readers.n;

    let (blob_sender, blob_receiver): (Sender<VariantBlob>, Receiver<VariantBlob>) = unbounded();
    let (result_sender, result_receiver): (
        Sender<VariantBlockResult>,
        Receiver<VariantBlockResult>,
    ) = unbounded();

    let reader_thread = thread::spawn(move || -> Result<()> {
        let mut readers = VcfReaders::new(vcf_paths.clone())?.readers;
        let tx = blob_sender;

        let pos_map = if let Some(ref target_positions) = args.target_positions {
            Some(positions_to_interval_trees(target_positions)?)
        } else {
            None
        };

        let mut all_bnd_variants: Vec<VariantInternal> = Vec::new();
        let sv_types_for_blob = [
            SvType::INSERTION,
            SvType::DELETION,
            SvType::INVERSION,
            SvType::DUPLICATION,
        ];

        log::debug!("Reader thread started.");
        for contig_name in contig_order {
            match load_contig(
                &contig_name,
                &args.merge_args,
                &bed_map,
                &pos_map,
                &mut readers,
            ) {
                Ok(mut variants_by_type_for_contig) => {
                    for (type_index, variants) in variants_by_type_for_contig.drain(..).enumerate()
                    {
                        if variants.is_empty() {
                            continue;
                        }

                        if type_index < sv_types_for_blob.len() {
                            // Non-BND types
                            let variant_type = sv_types_for_blob[type_index];
                            let blob = VariantBlob {
                                variants,
                                contig: contig_name.clone(),
                                variant_type,
                                args: args.merge_args.clone(),
                            };
                            log::debug!(
                                "Reader: Sending blob for Contig={}, Type={}, Variants={}",
                                blob.contig,
                                blob.variant_type,
                                blob.variants.len()
                            );
                            if tx.send(blob).is_err() {
                                log::error!(
                                    "Failed to send variant blob for contig {} / type {}. Receiver closed.",
                                    contig_name,
                                    variant_type
                                );
                                return Err(anyhow!(
                                    "Channel receiver closed unexpectedly in reader thread"
                                ));
                            }
                        } else {
                            // BND type (index 4)
                            log::debug!(
                                "Reader: Collected {} BND variants from contig {}",
                                variants.len(),
                                contig_name
                            );
                            all_bnd_variants.extend(variants);
                        }
                    }
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
        // TODO: Process all_bnd_variants
        // for v in all_bnd_variants {
        //     println!("{}", v);
        // }

        Ok(())
    });

    let writer_thread = thread::spawn(move || -> Result<()> {
        let mut writer_instance = writer;
        log::debug!("Writer thread started.");

        for payload in result_receiver {
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
            write_variants(payload, &sample_mapping, n_samples, &mut writer_instance)?;
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
        .thread_name(|i| format!("svx-merge-{}", i))
        .build()
        .map_err(|e| anyhow!("Failed to initialize merge thread pool: {}", e))?;

    pool.install(|| {
        blob_receiver
            .into_iter()
            .par_bridge()
            .for_each_with(result_sender, |s, variant_blob| {
                process_blob(variant_blob, s);
            });
    });

    if let Err(e) = reader_thread.join().unwrap() {
        return Err(anyhow!("Reader thread failed: {}", e));
    }
    log::debug!("Reader thread joined successfully.");

    if let Err(e) = writer_thread.join().unwrap() {
        return Err(anyhow!("Writer thread failed: {}", e));
    }
    log::debug!("Writer thread joined successfully.");

    Ok(())
}

fn process_blob(variant_blob: VariantBlob, sender: &Sender<VariantBlockResult>) {
    log::debug!(
        "Worker [{}]: Processing blob: Contig={}, Type={}, Variants={}",
        std::thread::current().name().unwrap_or("unnamed"),
        variant_blob.contig,
        variant_blob.variant_type,
        variant_blob.variants.len()
    );

    if variant_blob.variants.is_empty() {
        log::debug!(
            "Worker [{}]: Received empty variant blob for {} - {}, skipping.",
            std::thread::current().name().unwrap_or("unnamed"),
            variant_blob.contig,
            variant_blob.variant_type
        );
        return;
    }

    let dump = variant_blob.args.dump;
    let mut variant_block = VariantBlock::new(variant_blob);
    log::debug!(
        "Worker [{}]: Starting merge for {}/{}",
        std::thread::current().name().unwrap_or("unnamed"),
        variant_block.contig,
        variant_block.variant_type
    );
    variant_block.merge_block();
    log::debug!(
        "Worker [{}]: Merge complete for {}/{}, getting groups",
        std::thread::current().name().unwrap_or("unnamed"),
        variant_block.contig,
        variant_block.variant_type
    );
    let variant_block_result = variant_block.get_groups();

    if dump {
        print_group_stats(&variant_block_result);
    }

    log::debug!(
        "Worker [{}]: Sending result to writer: contig={}, type={}, groups={}",
        std::thread::current().name().unwrap_or("unnamed"),
        variant_block_result.contig,
        variant_block_result.variant_type,
        variant_block_result.groups.len()
    );
    if let Err(e) = sender.send(variant_block_result) {
        log::error!(
            "Worker [{}]: Failed to send result to writer thread: {}",
            std::thread::current().name().unwrap_or("unnamed"),
            e
        );
    }
}
