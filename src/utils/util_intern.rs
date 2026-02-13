use super::util::Result;

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
    let status_contents = std::fs::read_to_string("/proc/self/status")
        .map_err(|error| crate::svx_error!("Failed to read /proc/self/status: {error}"))?;
    parse_linux_peak_rss_bytes_from_status(&status_contents)
}

#[cfg(target_os = "macos")]
pub fn peak_memory_usage() -> Result<usize> {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        match retval {
            0 => Ok(rusage.ru_maxrss as usize),
            _ => Err(crate::svx_error!("libc::getrusage call failed")),
        }
    }
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
pub fn peak_memory_usage() -> Result<usize> {
    Err(crate::svx_error!(
        "No peak_memory_usage implementation for this OS"
    ))
}

#[cfg(any(target_os = "linux", test))]
fn parse_linux_peak_rss_bytes_from_status(status: &str) -> Result<usize> {
    let vmhwm_line = status
        .lines()
        .find(|line| line.trim_start().starts_with("VmHWM:"))
        .ok_or_else(|| crate::svx_error!("VmHWM missing in /proc/self/status"))?;
    parse_linux_peak_rss_bytes_from_vmhwm_line(vmhwm_line)
}

#[cfg(any(target_os = "linux", test))]
fn parse_linux_peak_rss_bytes_from_vmhwm_line(vmhwm_line: &str) -> Result<usize> {
    let vmhwm_value_and_unit = vmhwm_line
        .trim_start()
        .strip_prefix("VmHWM:")
        .ok_or_else(|| crate::svx_error!("Invalid VmHWM line in /proc/self/status"))?;
    let mut fields = vmhwm_value_and_unit.split_whitespace();
    let value_kib = fields
        .next()
        .ok_or_else(|| crate::svx_error!("VmHWM value missing in /proc/self/status"))?
        .parse::<usize>()
        .map_err(|error| crate::svx_error!("Invalid VmHWM value in /proc/self/status: {error}"))?;
    let unit = fields
        .next()
        .ok_or_else(|| crate::svx_error!("VmHWM unit missing in /proc/self/status"))?;
    if unit != "kB" {
        return Err(crate::svx_error!(
            "Unexpected VmHWM unit in /proc/self/status: {unit}"
        ));
    }
    value_kib
        .checked_mul(1024)
        .ok_or_else(|| crate::svx_error!("VmHWM value in /proc/self/status overflows usize"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_linux_peak_rss_bytes_from_status_reads_vmhwm() {
        let status = "\
Name:\tsvx
State:\tR (running)
VmRSS:\t   512 kB
VmHWM:\t 12345 kB
Threads:\t8
";
        assert_eq!(
            parse_linux_peak_rss_bytes_from_status(status).unwrap(),
            12_641_280
        );
    }

    #[test]
    fn parse_linux_peak_rss_bytes_from_status_errors_when_missing_vmhwm() {
        let status = "\
Name:\tsvx
State:\tR (running)
VmRSS:\t   512 kB
Threads:\t8
";
        assert!(
            parse_linux_peak_rss_bytes_from_status(status)
                .unwrap_err()
                .to_string()
                .contains("VmHWM")
        );
    }

    #[test]
    fn parse_linux_peak_rss_bytes_from_status_errors_on_invalid_vmhwm_value() {
        let status = "\
Name:\tsvx
VmHWM:\tnot_a_number kB
";
        assert!(
            parse_linux_peak_rss_bytes_from_status(status)
                .unwrap_err()
                .to_string()
                .contains("Invalid VmHWM")
        );
    }

    #[cfg(target_os = "macos")]
    #[test]
    fn peak_memory_usage_macos_is_non_zero() {
        assert!(peak_memory_usage().unwrap() > 0);
    }

    #[cfg(target_os = "macos")]
    #[test]
    fn peak_memory_usage_macos_is_monotonic() {
        let before = peak_memory_usage().unwrap();
        let mut buffer = vec![0u8; 8 * 1024 * 1024];
        buffer.iter_mut().step_by(4096).for_each(|byte| *byte = 1);
        std::hint::black_box(buffer);
        let after = peak_memory_usage().unwrap();
        assert!(after >= before);
    }
}
