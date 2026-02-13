use crate::{error::SvxError, utils::util::Result};
use crossbeam_channel::Sender;
use std::{any::Any, thread};

use super::progress::ProgressEvent;

fn panic_payload_message(panic_payload: &(dyn Any + Send + 'static)) -> String {
    if let Some(message) = panic_payload.downcast_ref::<&str>() {
        return (*message).to_owned();
    }
    if let Some(message) = panic_payload.downcast_ref::<String>() {
        return message.clone();
    }
    "unknown panic payload".to_owned()
}

fn join_thread_result(thread_name: &str, handle: thread::JoinHandle<Result<()>>) -> Result<()> {
    match handle.join() {
        Ok(result) => result.map_err(|e| crate::svx_error!("{thread_name} thread failed: {e}")),
        Err(panic_payload) => Err(crate::svx_error!(
            "{thread_name} thread panicked: {}",
            panic_payload_message(panic_payload.as_ref())
        )),
    }
}

fn join_progress_thread_result(handle: thread::JoinHandle<()>) -> Result<()> {
    handle.join().map_err(|panic_payload| {
        crate::svx_error!(
            "Progress thread panicked: {}",
            panic_payload_message(panic_payload.as_ref())
        )
    })
}

fn aggregate_shutdown_errors(errors: Vec<SvxError>) -> Result<()> {
    if errors.is_empty() {
        return Ok(());
    }
    if errors.len() == 1 {
        return Err(errors.into_iter().next().expect("single error expected"));
    }
    let summary = errors
        .into_iter()
        .enumerate()
        .map(|(index, error)| format!("{}. {}", index + 1, error))
        .collect::<Vec<_>>()
        .join("; ");
    Err(crate::svx_error!(
        "Multiple thread shutdown errors: {summary}"
    ))
}

pub(crate) fn finalize_merge_threads(
    reader_thread: thread::JoinHandle<Result<()>>,
    writer_thread: thread::JoinHandle<Result<()>>,
    progress_sender: Option<Sender<ProgressEvent>>,
    progress_thread: Option<thread::JoinHandle<()>>,
) -> Result<()> {
    let mut errors = Vec::new();

    match join_thread_result("Reader", reader_thread) {
        Ok(()) => log::debug!("Reader thread joined successfully."),
        Err(error) => errors.push(error),
    }
    match join_thread_result("Writer", writer_thread) {
        Ok(()) => log::debug!("Writer thread joined successfully."),
        Err(error) => errors.push(error),
    }

    drop(progress_sender);

    if let Some(handle) = progress_thread {
        match join_progress_thread_result(handle) {
            Ok(()) => {}
            Err(error) => {
                log::debug!("Progress thread terminated unexpectedly.");
                errors.push(error);
            }
        }
    }

    aggregate_shutdown_errors(errors)
}
