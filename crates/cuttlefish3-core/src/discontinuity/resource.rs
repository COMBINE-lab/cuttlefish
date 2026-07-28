//! Process resource discovery used by adaptive external-memory planning.
//!
//! This module is intentionally independent of graph data structures. The
//! public profiling helpers are re-exported from [`super`], preserving the
//! original API while keeping platform-specific code out of the algorithm.

#[cfg(all(target_os = "linux", target_env = "gnu"))]
unsafe extern "C" {
    fn malloc_trim(pad: usize) -> i32;
    fn getrlimit(resource: i32, limits: *mut RLimit) -> i32;
    fn setrlimit(resource: i32, limits: *const RLimit) -> i32;
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
#[repr(C)]
struct RLimit {
    current: u64,
    maximum: u64,
}

/// Returns the process soft file-descriptor limit.
///
/// The conservative fallback is used on unsupported platforms and if querying
/// the process limit fails.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
pub(crate) fn open_file_limit() -> usize {
    const RLIMIT_NOFILE: i32 = 7;
    unsafe {
        let mut limits = RLimit {
            current: 0,
            maximum: 0,
        };
        if getrlimit(RLIMIT_NOFILE, &mut limits) != 0 {
            return 1024;
        }
        usize::try_from(limits.current).unwrap_or(usize::MAX)
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
pub(crate) fn open_file_limit() -> usize {
    1024
}

/// Raises the soft descriptor limit to the hard limit, returning the result.
///
/// A process may always raise its soft limit up to its hard limit without
/// privileges, and there is no reason for this one not to: every descriptor it
/// opens is a file it is actively using, and the phases that fan out already
/// adapt their width to whatever budget they are given. Doing it here means a
/// user never has to discover `ulimit -n` to build a large graph -- the stock
/// soft limit of 1024 on most distributions is enough to run but tight enough
/// to force the fanout planners to narrow, and on a systemd host the hard
/// limit is usually 524288.
///
/// Returns `(soft_before, soft_after)`.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
pub fn raise_open_file_limit() -> (usize, usize) {
    const RLIMIT_NOFILE: i32 = 7;
    unsafe {
        let mut limits = RLimit {
            current: 0,
            maximum: 0,
        };
        if getrlimit(RLIMIT_NOFILE, &mut limits) != 0 {
            return (0, 0);
        }
        let before = usize::try_from(limits.current).unwrap_or(usize::MAX);
        if limits.current >= limits.maximum {
            return (before, before);
        }
        let raised = RLimit {
            current: limits.maximum,
            maximum: limits.maximum,
        };
        if setrlimit(RLIMIT_NOFILE, &raised) != 0 {
            return (before, before);
        }
        (before, usize::try_from(limits.maximum).unwrap_or(usize::MAX))
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
pub fn raise_open_file_limit() -> (usize, usize) {
    (0, 0)
}

/// Returns the number of descriptors currently owned by this process.
#[cfg(target_os = "linux")]
pub(crate) fn current_open_file_count() -> usize {
    std::fs::read_dir("/proc/self/fd")
        .map(|entries| entries.count())
        .unwrap_or(3)
}

#[cfg(not(target_os = "linux"))]
pub(crate) fn current_open_file_count() -> usize {
    3
}

/// Asks glibc to return free allocator arenas to the operating system.
///
/// This is a no-op on platforms that do not expose `malloc_trim`. It is used at
/// phase boundaries only and is never called from graph hot loops.
pub fn trim_process_allocations() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        let _ = malloc_trim(0);
    }
}

/// Reports current and peak RSS when `CF3_RS_PROFILE_RSS` is set.
///
/// Profiling is best-effort and has no effect when process status information
/// is unavailable.
pub fn report_process_memory(label: &str) {
    if std::env::var_os("CF3_RS_PROFILE_RSS").is_none() {
        return;
    }
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return;
    };
    let mut rss_kib = None;
    let mut hwm_kib = None;
    for line in status.lines() {
        if let Some(value) = line.strip_prefix("VmRSS:") {
            rss_kib = parse_status_kib(value);
        } else if let Some(value) = line.strip_prefix("VmHWM:") {
            hwm_kib = parse_status_kib(value);
        }
    }
    if let (Some(rss), Some(hwm)) = (rss_kib, hwm_kib) {
        eprintln!("cuttlefish3-rs: rss {label}: current {rss} KiB, peak {hwm} KiB");
    }
}

fn parse_status_kib(value: &str) -> Option<u64> {
    value
        .split_whitespace()
        .next()
        .and_then(|value| value.parse::<u64>().ok())
}

#[cfg(test)]
mod tests {
    use super::parse_status_kib;

    #[test]
    fn parses_linux_status_values() {
        assert_eq!(parse_status_kib("  12345 kB"), Some(12_345));
        assert_eq!(parse_status_kib(" unavailable"), None);
    }
}
