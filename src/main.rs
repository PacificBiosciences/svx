use anyhow::Result;
use std::time;
use svx::{
    cli::{Command, FULL_VERSION, init_verbose, parse_cli_with_config},
    commands::merge,
    utils::util::handle_error_and_exit,
};

#[allow(unused_imports)]
use svx::utils::util_intern::readable_size;

#[cfg(any(target_os = "linux", target_os = "macos"))]
use svx::utils::util_intern::peak_memory_usage;

fn runner() -> Result<()> {
    let cli = parse_cli_with_config();
    init_verbose(&cli);
    log::trace!("CLI options set: {:?}", cli);

    log::info!(
        "Running {}-{} [{}]",
        env!("CARGO_PKG_NAME"),
        FULL_VERSION,
        cli.command.name()
    );

    let start_timer = time::Instant::now();
    match cli.command {
        Command::Merge(args) => {
            log::trace!("Merge arguments: {:#?}", args);
            merge(args)?
        }
    }
    log::info!("Total execution time: {:.2?}", start_timer.elapsed());

    #[cfg(any(target_os = "linux", target_os = "macos"))]
    {
        let (size, unit) = readable_size(peak_memory_usage()?);
        log::info!("Peak memory use: {:.2} {}", size, unit);
    }

    log::info!("{} end", env!("CARGO_PKG_NAME"));
    Ok(())
}

fn main() {
    if let Err(e) = runner() {
        handle_error_and_exit(e);
    }
}
