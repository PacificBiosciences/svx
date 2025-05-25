use super::util::Result;
#[allow(unused_imports)]
use anyhow::anyhow;

pub fn readable_size(bytes: usize) -> (f64, &'static str) {
    let units: Vec<(f64, &'static str)> = vec![
        (1.0, "B"),
        (1024.0, "KiB"),
        (1024.0 * 1024.0, "MiB"),
        (1024.0 * 1024.0 * 1024.0, "GiB"),
        (1024.0 * 1024.0 * 1024.0 * 1024.0, "TiB"),
    ];

    let value = bytes as f64;
    let mut unit = units[0];
    for next in units.iter().skip(1) {
        if value >= next.0 {
            unit = *next;
        } else {
            break;
        }
    }

    (value / unit.0, unit.1)
}

#[cfg(target_os = "linux")]
pub fn peak_memory_usage() -> Result<usize> {
    // unsafe {
    //     let mut rusage: libc::rusage = std::mem::zeroed();
    //     let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
    //     match retval {
    //         0 => Ok(rusage.ru_maxrss as usize * 1024),
    //         _ => Err("libc::getrusage call failed"),
    //     }
    // }
    Ok(0)
}

#[cfg(target_os = "macos")]
pub fn peak_memory_usage() -> Result<usize> {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        match retval {
            0 => Ok(rusage.ru_maxrss as usize),
            _ => Err(anyhow!("libc::getrusage call failed")),
        }
    }
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
pub fn peak_memory_usage() -> Result<usize, &'static str> {
    Err("No peak_memory_usage implementation for this OS")
}
