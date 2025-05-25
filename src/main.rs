use clap::Parser;
use std::time;
use svx::{
    cli::{init_verbose, Cli, Command, FULL_VERSION},
    commands::merge,
    utils::util::{handle_error_and_exit, Result},
};

#[allow(unused_imports)]
use svx::utils::util_intern::readable_size;

#[cfg(target_os = "macos")]
use svx::utils::util_intern::peak_memory_usage;

fn runner() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    log::trace!("CLI options set: {:?}", cli);

    log::info!(
        "Running {}-{} [{}]",
        env!("CARGO_PKG_NAME"),
        &**FULL_VERSION,
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

    #[cfg(target_os = "macos")]
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
