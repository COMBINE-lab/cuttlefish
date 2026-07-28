//! Process resource discovery used by adaptive external-memory planning.
//!
//! This module is intentionally independent of graph data structures. The
//! public profiling helpers are re-exported from [`super`], preserving the
//! original API while keeping platform-specific code out of the algorithm.

/// Returns the process soft file-descriptor limit.
///
/// `libc::RLIMIT_NOFILE` rather than a literal: the resource number is 7 on
/// most Linux architectures, 6 on sparc and 8 on Apple, so a hardcoded value
/// silently queries the wrong limit.
#[cfg(unix)]
pub(crate) fn open_file_limit() -> usize {
    // SAFETY: `limits` is a valid, fully initialized destination.
    unsafe {
        let mut limits = std::mem::zeroed::<libc::rlimit>();
        if libc::getrlimit(libc::RLIMIT_NOFILE, &mut limits) != 0 {
            return 1024;
        }
        usize::try_from(limits.rlim_cur).unwrap_or(usize::MAX)
    }
}

#[cfg(not(unix))]
pub(crate) fn open_file_limit() -> usize {
    1024
}

/// Raises the soft descriptor limit to the hard limit, returning the result.
///
/// A process may always raise its soft limit up to its hard limit without
/// privileges, and there is no reason for this one not to: every descriptor it
/// opens is a file it is actively using, and the phases that fan out already
/// adapt their width to whatever budget they are given. Doing it here means a
/// user never has to discover `ulimit -n` to build a large graph. It matters
/// more on macOS than on Linux, where the stock soft limit is 256 rather than
/// 1024 -- close enough to this build's floor to force heavy narrowing.
///
/// Returns `(soft_before, soft_after)`.
#[cfg(unix)]
pub fn raise_open_file_limit() -> (usize, usize) {
    // SAFETY: both calls take valid, fully initialized `rlimit` pointers.
    unsafe {
        let mut limits = std::mem::zeroed::<libc::rlimit>();
        if libc::getrlimit(libc::RLIMIT_NOFILE, &mut limits) != 0 {
            return (0, 0);
        }
        let before = usize::try_from(limits.rlim_cur).unwrap_or(usize::MAX);
        // Darwin reports an effectively infinite hard limit but refuses a soft
        // limit above `sysconf(_SC_OPEN_MAX)`, so asking for the hard limit
        // directly fails with EINVAL and leaves the process at 256.
        let ceiling = darwin_open_max().unwrap_or(limits.rlim_max);
        let target = limits.rlim_max.min(ceiling);
        if target <= limits.rlim_cur {
            return (before, before);
        }
        let raised = libc::rlimit {
            rlim_cur: target,
            rlim_max: limits.rlim_max,
        };
        if libc::setrlimit(libc::RLIMIT_NOFILE, &raised) != 0 {
            return (before, before);
        }
        (before, usize::try_from(target).unwrap_or(usize::MAX))
    }
}

#[cfg(target_vendor = "apple")]
fn darwin_open_max() -> Option<libc::rlim_t> {
    // SAFETY: `sysconf` takes a name and returns a long; -1 signals no limit.
    let value = unsafe { libc::sysconf(libc::_SC_OPEN_MAX) };
    (value > 0).then(|| value as libc::rlim_t)
}

#[cfg(not(target_vendor = "apple"))]
fn darwin_open_max() -> Option<u64> {
    None
}

#[cfg(not(unix))]
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
        let _ = libc::malloc_trim(0);
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
