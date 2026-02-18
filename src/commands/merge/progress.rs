use crossbeam_channel::{Receiver, Sender};
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use log::LevelFilter;
use std::{
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
    time::Duration,
};

pub enum ProgressEvent {
    ContigCompleted { contig: String },
    ContigScanDone,
    VariantsQueued { variants: u64 },
    ReaderDone,
    VariantsMerged { variants: u64 },
    WriterAdvanced { variants: u64, records: u64 },
    WriterSortFinalizeStart,
    WriterSortFinalizeDone,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct QueueDepthSnapshot {
    pub current: usize,
    pub peak: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PipelineQueueSnapshot {
    pub blob: QueueDepthSnapshot,
    pub result: QueueDepthSnapshot,
}

#[derive(Debug, Default)]
pub struct QueueDepthTracker {
    current: AtomicUsize,
    peak: AtomicUsize,
}

impl QueueDepthTracker {
    pub fn increment(&self) -> usize {
        let current = self.current.fetch_add(1, Ordering::Relaxed) + 1;
        self.update_peak(current);
        current
    }

    pub fn decrement(&self) -> usize {
        let mut observed = self.current.load(Ordering::Relaxed);
        loop {
            if observed == 0 {
                return 0;
            }
            match self.current.compare_exchange_weak(
                observed,
                observed - 1,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => return observed - 1,
                Err(next) => observed = next,
            }
        }
    }

    pub fn snapshot(&self) -> QueueDepthSnapshot {
        QueueDepthSnapshot {
            current: self.current.load(Ordering::Relaxed),
            peak: self.peak.load(Ordering::Relaxed),
        }
    }

    fn update_peak(&self, current: usize) {
        let mut peak = self.peak.load(Ordering::Relaxed);
        while current > peak {
            match self.peak.compare_exchange_weak(
                peak,
                current,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(next) => peak = next,
            }
        }
    }
}

#[derive(Debug, Default)]
pub struct PipelineQueueMetrics {
    pub blob: QueueDepthTracker,
    pub result: QueueDepthTracker,
}

impl PipelineQueueMetrics {
    pub fn snapshot(&self) -> PipelineQueueSnapshot {
        PipelineQueueSnapshot {
            blob: self.blob.snapshot(),
            result: self.result.snapshot(),
        }
    }
}

pub fn compute_progress_enabled(
    force_progress: bool,
    disable_progress: bool,
    is_stderr_tty: bool,
    max_level: LevelFilter,
) -> bool {
    if disable_progress || !is_stderr_tty {
        return false;
    }
    if matches!(max_level, LevelFilter::Debug | LevelFilter::Trace) {
        return false;
    }
    if force_progress {
        return true;
    }
    matches!(
        max_level,
        LevelFilter::Info | LevelFilter::Warn | LevelFilter::Error
    )
}

pub fn compute_pipeline_queue_capacities(num_threads: usize) -> (usize, usize) {
    let blob_queue_capacity = num_threads.saturating_mul(2).clamp(4, 64);
    let result_queue_capacity = num_threads.saturating_mul(4).clamp(8, 128);
    (blob_queue_capacity, result_queue_capacity)
}

pub fn resolve_pipeline_queue_capacities(
    num_threads: usize,
    blob_override: Option<usize>,
    result_override: Option<usize>,
) -> (usize, usize) {
    let (default_blob, default_result) = compute_pipeline_queue_capacities(num_threads);
    (
        blob_override.unwrap_or(default_blob),
        result_override.unwrap_or(default_result),
    )
}

pub fn send_progress_event(progress_sender: Option<&Sender<ProgressEvent>>, event: ProgressEvent) {
    if let Some(sender) = progress_sender {
        if sender.send(event).is_err() {
            log::debug!("Progress channel receiver closed unexpectedly.");
        }
    }
}

fn contig_progress_style() -> ProgressStyle {
    ProgressStyle::with_template("{prefix:>8} [{bar:40.cyan/blue}] {pos:>4}/{len:4} {msg}")
        .expect("contig progress style must be valid")
        .progress_chars("=> ")
}

fn merge_spinner_style() -> ProgressStyle {
    ProgressStyle::with_template("{prefix:>8} {spinner} {msg}")
        .expect("merge spinner style must be valid")
}

fn merge_progress_style() -> ProgressStyle {
    ProgressStyle::with_template("{prefix:>8} [{bar:40.green/blue}] {pos:>8}/{len:8} {msg}")
        .expect("merge progress style must be valid")
        .progress_chars("=> ")
}

pub fn build_merge_progress_message(
    stage: &str,
    queued_variants: u64,
    merged_variants: u64,
    written_variants: u64,
    written_records: u64,
    queue_snapshot: Option<PipelineQueueSnapshot>,
) -> String {
    let base = format!(
        "stage={stage} read={queued_variants} merged={merged_variants} written={written_variants} records={written_records}"
    );
    if let Some(snapshot) = queue_snapshot {
        format!(
            "{base} blob_q={}/{} result_q={}/{}",
            snapshot.blob.current,
            snapshot.blob.peak,
            snapshot.result.current,
            snapshot.result.peak
        )
    } else {
        base
    }
}

fn update_merge_progress_message(
    merge_bar: &ProgressBar,
    stage: &str,
    queued_variants: u64,
    merged_variants: u64,
    written_variants: u64,
    written_records: u64,
    queue_snapshot: Option<PipelineQueueSnapshot>,
) {
    merge_bar.set_message(build_merge_progress_message(
        stage,
        queued_variants,
        merged_variants,
        written_variants,
        written_records,
        queue_snapshot,
    ));
}

fn merge_stage_label(
    contig_scan_done: bool,
    reader_done: bool,
    sort_finalize_active: bool,
    queued_variants: u64,
    merged_variants: u64,
    written_variants: u64,
) -> &'static str {
    if sort_finalize_active {
        return "final_sort_flush";
    }
    if !contig_scan_done {
        return "processing";
    }
    if !reader_done {
        return "finalizing_bnd_workload";
    }
    if merged_variants < queued_variants {
        return "processing_backlog";
    }
    if written_variants < queued_variants {
        return "flushing_output";
    }
    "final_drain"
}

pub fn run_progress_ui(
    progress_receiver: Receiver<ProgressEvent>,
    total_contigs: u64,
    queue_metrics: Option<Arc<PipelineQueueMetrics>>,
) {
    let multi = MultiProgress::with_draw_target(ProgressDrawTarget::stderr_with_hz(10));

    let contig_bar = multi.add(ProgressBar::new(total_contigs));
    contig_bar.set_prefix("contigs");
    contig_bar.set_style(contig_progress_style());
    contig_bar.set_message("reading");

    let merge_bar = multi.add(ProgressBar::new_spinner());
    merge_bar.set_prefix("variants");
    merge_bar.set_style(merge_spinner_style());
    merge_bar.enable_steady_tick(Duration::from_millis(100));

    let mut queued_variants = 0u64;
    let mut merged_variants = 0u64;
    let mut written_variants = 0u64;
    let mut written_records = 0u64;
    let mut contig_scan_done = total_contigs == 0;
    let mut reader_done = false;
    let mut sort_finalize_active = false;
    let mut merge_bar_in_progress_mode = false;

    if contig_scan_done {
        merge_bar.set_style(merge_progress_style());
        merge_bar.set_length(queued_variants);
        merge_bar.set_position(0);
        merge_bar_in_progress_mode = true;
    }

    let queue_snapshot = queue_metrics.as_ref().map(|metrics| metrics.snapshot());
    let stage = merge_stage_label(
        contig_scan_done,
        reader_done,
        sort_finalize_active,
        queued_variants,
        merged_variants,
        written_variants,
    );
    update_merge_progress_message(
        &merge_bar,
        stage,
        queued_variants,
        merged_variants,
        written_variants,
        written_records,
        queue_snapshot,
    );

    for event in progress_receiver {
        match event {
            ProgressEvent::ContigCompleted { contig } => {
                contig_bar.inc(1);
                contig_bar.set_message(contig);
            }
            ProgressEvent::ContigScanDone => {
                contig_scan_done = true;
                contig_bar.set_message("contig scan complete");
                if !merge_bar_in_progress_mode && !sort_finalize_active {
                    merge_bar.set_style(merge_progress_style());
                    merge_bar.set_length(queued_variants);
                    merge_bar.set_position(merged_variants.min(queued_variants));
                    merge_bar_in_progress_mode = true;
                }
            }
            ProgressEvent::VariantsQueued { variants } => {
                queued_variants += variants;
                if merge_bar_in_progress_mode {
                    merge_bar.set_length(queued_variants);
                }
            }
            ProgressEvent::ReaderDone => {
                reader_done = true;
                if !merge_bar_in_progress_mode && !sort_finalize_active {
                    merge_bar.set_style(merge_progress_style());
                    merge_bar.set_length(queued_variants);
                    merge_bar.set_position(merged_variants.min(queued_variants));
                    merge_bar_in_progress_mode = true;
                } else if merge_bar_in_progress_mode {
                    merge_bar.set_length(queued_variants);
                }
            }
            ProgressEvent::VariantsMerged { variants } => {
                merged_variants += variants;
            }
            ProgressEvent::WriterAdvanced { variants, records } => {
                written_variants += variants;
                written_records += records;
            }
            ProgressEvent::WriterSortFinalizeStart => {
                sort_finalize_active = true;
                merge_bar.set_style(merge_spinner_style());
                merge_bar_in_progress_mode = false;
            }
            ProgressEvent::WriterSortFinalizeDone => {
                sort_finalize_active = false;
                merge_bar.set_style(merge_progress_style());
                merge_bar.set_length(queued_variants);
                merge_bar.set_position(merged_variants.min(queued_variants));
                merge_bar_in_progress_mode = true;
            }
        }

        if merge_bar_in_progress_mode {
            merge_bar.set_position(merged_variants.min(queued_variants));
        }

        let queue_snapshot = queue_metrics.as_ref().map(|metrics| metrics.snapshot());
        let stage = merge_stage_label(
            contig_scan_done,
            reader_done,
            sort_finalize_active,
            queued_variants,
            merged_variants,
            written_variants,
        );
        update_merge_progress_message(
            &merge_bar,
            stage,
            queued_variants,
            merged_variants,
            written_variants,
            written_records,
            queue_snapshot,
        );
    }

    if !merge_bar_in_progress_mode {
        merge_bar.set_style(merge_progress_style());
        merge_bar.set_length(queued_variants);
        merge_bar.set_position(merged_variants.min(queued_variants));
    }

    let queue_snapshot = queue_metrics.as_ref().map(|metrics| metrics.snapshot());
    let stage = merge_stage_label(
        contig_scan_done,
        reader_done,
        sort_finalize_active,
        queued_variants,
        merged_variants,
        written_variants,
    );
    contig_bar.finish_with_message("done");
    merge_bar.finish_with_message(build_merge_progress_message(
        stage,
        queued_variants,
        merged_variants,
        written_variants,
        written_records,
        queue_snapshot,
    ));
    let _ = multi.clear();
}
