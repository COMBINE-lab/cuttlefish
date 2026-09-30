//! How fast does a directory's storage absorb and return a stream, and what
//! kind of storage is it?
//!
//! Programs with large temporary files -- external sorts, index builders,
//! graph constructions -- can trade CPU for I/O, for instance by compressing
//! what they spill. Whether that pays depends on the storage under the work
//! directory, which the program usually cannot know in advance. This crate
//! answers two questions about a directory, cheaply enough to ask at startup:
//!
//! * [`classify`] reads metadata only (no I/O): the filesystem, whether it is
//!   local, networked or memory-backed, and whether the device under it
//!   rotates.
//! * [`probe`] also times a short sequential write and read-back of random
//!   bytes, capped in both bytes and time, so a fast device costs a fraction
//!   of a second and a slow one stops early. It bypasses the page cache
//!   (`O_DIRECT` on Linux, `F_NOCACHE` on macOS), so the numbers describe the
//!   device rather than memory.
//!
//! Linux and macOS are supported. Elsewhere [`classify`] answers
//! [`StorageKind::Unknown`] and [`probe`] measures buffered I/O.
//!
//! The crate reports; policy stays with the caller. [`Probe::sustains`] is a
//! small helper for the usual question -- can this storage keep up with a
//! given stream rate, with a safety margin? -- but deciding what rate a
//! program produces is the program's business.
//!
//! A probe is a snapshot of idle storage. Shared or networked storage can be
//! busier later, and a probe cannot see what other processes will do.

use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// Caps on the timed part of a [`probe`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Limits {
    /// Most bytes written, then read back.
    pub max_bytes: u64,
    /// Writing, and then reading, each stop once this much time has passed,
    /// after at least one chunk; so a probe takes about twice this, plus a
    /// chunk each way on slow storage.
    pub max_time: Duration,
    /// Bytes per write and read call; rounded up to 4 KiB.
    pub chunk_bytes: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_bytes: 256 * 1024 * 1024,
            max_time: Duration::from_millis(500),
            chunk_bytes: 4 * 1024 * 1024,
        }
    }
}

/// Where a filesystem keeps its data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StorageKind {
    /// A local block device.
    Local,
    /// Across a network: NFS, SMB, Lustre, GPFS, CephFS, BeeGFS and the like.
    Network,
    /// Memory: tmpfs or ramfs.
    Memory,
    /// Not recognised, including FUSE filesystems, which may be either.
    Unknown,
}

/// What [`classify`] learns from metadata alone.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Storage {
    /// The filesystem's name, or `"unknown"`.
    pub filesystem: &'static str,
    pub kind: StorageKind,
    /// Whether the block device under the directory rotates; for a
    /// device-mapper or md device, whether any device under it does. `None`
    /// when there is no local block device to ask about, and always on
    /// macOS, where the answer lives in IOKit.
    pub rotational: Option<bool>,
}

/// Bytes moved in a measured time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Throughput {
    pub bytes: u64,
    pub elapsed: Duration,
}

impl Throughput {
    pub fn bytes_per_second(&self) -> f64 {
        self.bytes as f64 / self.elapsed.as_secs_f64().max(1e-9)
    }
}

/// The result of a [`probe`].
#[derive(Debug, Clone, PartialEq)]
pub struct Probe {
    pub dir: PathBuf,
    pub storage: Storage,
    /// Writing, including the final `fdatasync`.
    pub write: Throughput,
    /// Reading the same bytes back.
    pub read: Throughput,
    /// Whether both passes bypassed the page cache. Without direct I/O the
    /// read may be served from cache and overstate the device.
    pub direct_io: bool,
    /// `MemAvailable` from `/proc/meminfo`, in bytes, where there is one.
    pub mem_available: Option<u64>,
}

impl Probe {
    /// The slower of the write and read rates.
    pub fn slowest_bytes_per_second(&self) -> f64 {
        self.write
            .bytes_per_second()
            .min(self.read.bytes_per_second())
    }

    /// Whether both directions keep up with `bytes_per_second` times
    /// `margin`.
    pub fn sustains(&self, bytes_per_second: f64, margin: f64) -> bool {
        self.slowest_bytes_per_second() >= bytes_per_second * margin
    }
}

/// Classifies the storage under `dir` from metadata, without I/O.
pub fn classify(dir: &Path) -> Storage {
    #[cfg(target_os = "linux")]
    {
        linux::classify(dir)
    }
    #[cfg(target_os = "macos")]
    {
        macos::classify(dir)
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        let _ = dir;
        Storage {
            filesystem: "unknown",
            kind: StorageKind::Unknown,
            rotational: None,
        }
    }
}

/// `MemAvailable` from `/proc/meminfo`, in bytes, where there is one.
pub fn mem_available() -> Option<u64> {
    let meminfo = fs::read_to_string("/proc/meminfo").ok()?;
    let line = meminfo
        .lines()
        .find(|line| line.starts_with("MemAvailable:"))?;
    let kib: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    kib.checked_mul(1024)
}

/// Classifies `dir` and times a capped write and read-back in it.
///
/// Writes a hidden file in `dir` and removes it before returning, including
/// on error.
pub fn probe(dir: &Path, limits: &Limits) -> io::Result<Probe> {
    let chunk = limits.chunk_bytes.max(1).div_ceil(ALIGN) * ALIGN;
    let mut buffer = AlignedBuffer::new(chunk);
    fill_random(buffer.as_mut_slice());
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.subsec_nanos());
    let path = dir.join(format!(".scratch-probe-{}-{nanos}", std::process::id()));
    let _remove = RemoveOnDrop(&path);

    let (write, direct_io) = match timed_write(&path, &buffer, limits, true) {
        Ok(write) => (write, true),
        // Direct I/O is refused by some filesystems (tmpfs, several network
        // ones) at open or at the first write; measure buffered instead.
        Err(error) if direct_io_refused(&error) => {
            (timed_write(&path, &buffer, limits, false)?, false)
        }
        Err(error) => return Err(error),
    };
    let read = timed_read(&path, &mut buffer, limits, direct_io)?;
    Ok(Probe {
        dir: dir.to_path_buf(),
        storage: classify(dir),
        write,
        read,
        direct_io,
        mem_available: mem_available(),
    })
}

const ALIGN: usize = 4096;

fn timed_write(
    path: &Path,
    buffer: &AlignedBuffer,
    limits: &Limits,
    direct: bool,
) -> io::Result<Throughput> {
    let mut file = open(path, true, direct)?;
    let started = Instant::now();
    let mut written = 0u64;
    while written < limits.max_bytes && (written == 0 || started.elapsed() < limits.max_time) {
        file.write_all(buffer.as_slice())?;
        written += buffer.len() as u64;
        // Buffered writes land in the page cache at memory speed, so the
        // clock would stop long before the device had done any work, and
        // the final sync alone could take seconds. Syncing each chunk keeps
        // the time cap honest.
        if !direct {
            file.sync_data()?;
        }
    }
    file.sync_data()?;
    Ok(Throughput {
        bytes: written,
        elapsed: started.elapsed(),
    })
}

fn timed_read(
    path: &Path,
    buffer: &mut AlignedBuffer,
    limits: &Limits,
    direct: bool,
) -> io::Result<Throughput> {
    let mut file = open(path, false, direct)?;
    if !direct {
        drop_cached_pages(&file);
    }
    let started = Instant::now();
    let mut read = 0u64;
    while read == 0 || started.elapsed() < limits.max_time {
        match file.read(buffer.as_mut_slice())? {
            0 => break,
            n => read += n as u64,
        }
    }
    Ok(Throughput {
        bytes: read,
        elapsed: started.elapsed(),
    })
}

fn open(path: &Path, write: bool, direct: bool) -> io::Result<File> {
    let mut options = OpenOptions::new();
    if write {
        options.write(true).create(true).truncate(true);
    } else {
        options.read(true);
    }
    #[cfg(target_os = "linux")]
    if direct {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_DIRECT);
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    if direct {
        return Err(io::Error::from(io::ErrorKind::Unsupported));
    }
    let file = options.open(path)?;
    // macOS has no `O_DIRECT`; `F_NOCACHE` keeps this descriptor's reads and
    // writes out of the unified buffer cache, which is what the probe needs.
    #[cfg(target_os = "macos")]
    if direct {
        use std::os::unix::io::AsRawFd;
        // SAFETY: the descriptor is open for the duration of the call.
        if unsafe { libc::fcntl(file.as_raw_fd(), libc::F_NOCACHE, 1) } == -1 {
            return Err(io::Error::last_os_error());
        }
    }
    Ok(file)
}

fn direct_io_refused(error: &io::Error) -> bool {
    error.kind() == io::ErrorKind::Unsupported
        || error.raw_os_error() == Some(libc::EINVAL)
        || error.raw_os_error() == Some(libc::EOPNOTSUPP)
}

/// Asks the kernel to forget the file's cached pages, so a buffered read
/// comes from the device. Advisory: it may keep some.
fn drop_cached_pages(file: &File) {
    #[cfg(target_os = "linux")]
    {
        use std::os::unix::io::AsRawFd;
        // SAFETY: the descriptor is open for the duration of the call.
        unsafe {
            libc::posix_fadvise(file.as_raw_fd(), 0, 0, libc::POSIX_FADV_DONTNEED);
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = file;
}

/// Random bytes, so compressing or deduplicating storage cannot shrink them.
fn fill_random(bytes: &mut [u8]) {
    let mut state = 0x9E37_79B9_7F4A_7C15u64
        ^ std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |elapsed| elapsed.as_nanos() as u64);
    for chunk in bytes.chunks_mut(8) {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        chunk.copy_from_slice(&state.to_le_bytes()[..chunk.len()]);
    }
}

struct RemoveOnDrop<'a>(&'a Path);

impl Drop for RemoveOnDrop<'_> {
    fn drop(&mut self) {
        let _ = fs::remove_file(self.0);
    }
}

/// A heap buffer aligned for direct I/O.
struct AlignedBuffer {
    ptr: std::ptr::NonNull<u8>,
    len: usize,
}

impl AlignedBuffer {
    fn new(len: usize) -> Self {
        let layout = std::alloc::Layout::from_size_align(len, ALIGN).expect("valid layout");
        // SAFETY: `len` is a non-zero multiple of the alignment.
        let ptr = unsafe { std::alloc::alloc_zeroed(layout) };
        let ptr =
            std::ptr::NonNull::new(ptr).unwrap_or_else(|| std::alloc::handle_alloc_error(layout));
        Self { ptr, len }
    }

    fn len(&self) -> usize {
        self.len
    }

    fn as_slice(&self) -> &[u8] {
        // SAFETY: `ptr` owns `len` initialized bytes.
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    fn as_mut_slice(&mut self) -> &mut [u8] {
        // SAFETY: `ptr` owns `len` initialized bytes, borrowed mutably here.
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for AlignedBuffer {
    fn drop(&mut self) {
        let layout = std::alloc::Layout::from_size_align(self.len, ALIGN).expect("valid layout");
        // SAFETY: allocated in `new` with this layout.
        unsafe { std::alloc::dealloc(self.ptr.as_ptr(), layout) };
    }
}

#[cfg(target_os = "linux")]
mod linux {
    use super::{Storage, StorageKind};
    use std::path::{Path, PathBuf};

    /// Filesystem magic numbers (`statfs::f_type`) and what they mean.
    const FILESYSTEMS: &[(u32, &str, StorageKind)] = &[
        (0xEF53, "ext2/3/4", StorageKind::Local),
        (0x5846_5342, "xfs", StorageKind::Local),
        (0x9123_683E, "btrfs", StorageKind::Local),
        (0x2FC1_2FC1, "zfs", StorageKind::Local),
        (0xF2F5_2010, "f2fs", StorageKind::Local),
        (0x4D44, "vfat", StorageKind::Local),
        (0x5346_544E, "ntfs", StorageKind::Local),
        (0x794C_7630, "overlayfs", StorageKind::Unknown),
        (0x0102_1994, "tmpfs", StorageKind::Memory),
        (0x8584_58F6, "ramfs", StorageKind::Memory),
        (0x6969, "nfs", StorageKind::Network),
        (0x517B, "smb", StorageKind::Network),
        (0xFF53_4D42, "cifs", StorageKind::Network),
        (0xFE53_4D42, "smb2", StorageKind::Network),
        (0x0BD0_0BD0, "lustre", StorageKind::Network),
        (0x4750_4653, "gpfs", StorageKind::Network),
        (0x00C3_6400, "cephfs", StorageKind::Network),
        (0x1983_0326, "beegfs", StorageKind::Network),
        (0x5346_414F, "afs", StorageKind::Network),
        (0xAAD7_AAEA, "panfs", StorageKind::Network),
        (0x6573_5546, "fuse", StorageKind::Unknown),
    ];

    pub(super) fn classify(dir: &Path) -> Storage {
        let (filesystem, kind) = filesystem(dir);
        let rotational = match kind {
            StorageKind::Local | StorageKind::Unknown => rotational(dir),
            StorageKind::Network | StorageKind::Memory => None,
        };
        Storage {
            filesystem,
            kind,
            rotational,
        }
    }

    fn filesystem(dir: &Path) -> (&'static str, StorageKind) {
        use std::os::unix::ffi::OsStrExt;
        let Ok(path) = std::ffi::CString::new(dir.as_os_str().as_bytes()) else {
            return ("unknown", StorageKind::Unknown);
        };
        let mut stats = std::mem::MaybeUninit::<libc::statfs>::zeroed();
        // SAFETY: `path` is NUL-terminated and `stats` is a valid out-pointer.
        if unsafe { libc::statfs(path.as_ptr(), stats.as_mut_ptr()) } != 0 {
            return ("unknown", StorageKind::Unknown);
        }
        // SAFETY: `statfs` succeeded and initialized `stats`.
        let magic = unsafe { stats.assume_init() }.f_type as u32;
        FILESYSTEMS
            .iter()
            .find(|&&(known, _, _)| known == magic)
            .map_or(("unknown", StorageKind::Unknown), |&(_, name, kind)| {
                (name, kind)
            })
    }

    /// Follows the directory's device to its sysfs node and asks whether it
    /// rotates; stacked devices (device-mapper, md) rotate if any device under
    /// them does, and a partition answers for its disk.
    fn rotational(dir: &Path) -> Option<bool> {
        use std::os::unix::fs::MetadataExt;
        let device = std::fs::metadata(dir).ok()?.dev();
        let (major, minor) = (libc::major(device), libc::minor(device));
        if major == 0 {
            return None;
        }
        let node = std::fs::canonicalize(format!("/sys/dev/block/{major}:{minor}")).ok()?;
        device_rotational(&node, 0)
    }

    fn device_rotational(node: &Path, depth: usize) -> Option<bool> {
        if depth > 8 {
            return None;
        }
        if let Ok(slaves) = std::fs::read_dir(node.join("slaves")) {
            let answers: Vec<Option<bool>> = slaves
                .flatten()
                .filter_map(|slave| std::fs::canonicalize(slave.path()).ok())
                .map(|slave| device_rotational(&slave, depth + 1))
                .collect();
            if !answers.is_empty() {
                return if answers.contains(&Some(true)) {
                    Some(true)
                } else if answers.iter().all(|answer| *answer == Some(false)) {
                    Some(false)
                } else {
                    None
                };
            }
        }
        let queue = |dir: &Path| -> Option<bool> {
            let flag = std::fs::read_to_string(dir.join("queue/rotational")).ok()?;
            Some(flag.trim() == "1")
        };
        queue(node).or_else(|| {
            // A partition has no queue of its own; its disk is the parent.
            let parent: PathBuf = node.parent()?.to_path_buf();
            queue(&parent)
        })
    }
}

#[cfg(target_os = "macos")]
mod macos {
    use super::{Storage, StorageKind};
    use std::path::Path;

    /// `statfs::f_fstypename` values and what they mean.
    const FILESYSTEMS: &[(&str, StorageKind)] = &[
        ("apfs", StorageKind::Local),
        ("hfs", StorageKind::Local),
        ("msdos", StorageKind::Local),
        ("exfat", StorageKind::Local),
        ("nfs", StorageKind::Network),
        ("smbfs", StorageKind::Network),
        ("afpfs", StorageKind::Network),
        ("webdav", StorageKind::Network),
        ("tmpfs", StorageKind::Memory),
    ];

    pub(super) fn classify(dir: &Path) -> Storage {
        use std::os::unix::ffi::OsStrExt;
        let unknown = Storage {
            filesystem: "unknown",
            kind: StorageKind::Unknown,
            rotational: None,
        };
        let Ok(path) = std::ffi::CString::new(dir.as_os_str().as_bytes()) else {
            return unknown;
        };
        let mut stats = std::mem::MaybeUninit::<libc::statfs>::zeroed();
        // SAFETY: `path` is NUL-terminated and `stats` is a valid out-pointer.
        if unsafe { libc::statfs(path.as_ptr(), stats.as_mut_ptr()) } != 0 {
            return unknown;
        }
        // SAFETY: `statfs` succeeded and initialized `stats`.
        let stats = unsafe { stats.assume_init() };
        let name: Vec<u8> = stats
            .f_fstypename
            .iter()
            .take_while(|&&byte| byte != 0)
            .map(|&byte| byte as u8)
            .collect();
        let (filesystem, kind) = FILESYSTEMS
            .iter()
            .find(|&&(known, _)| known.as_bytes() == name.as_slice())
            .map_or_else(
                // An unlisted filesystem is still local or not; the kernel
                // says which. FUSE mounts (macFUSE) report themselves local
                // whatever they front, so they stay unknown.
                || {
                    let kind = if name.windows(4).any(|w| w == b"fuse") {
                        StorageKind::Unknown
                    } else if stats.f_flags & libc::MNT_LOCAL as u32 != 0 {
                        StorageKind::Local
                    } else {
                        StorageKind::Network
                    };
                    ("unknown", kind)
                },
                |&(name, kind)| (name, kind),
            );
        Storage {
            filesystem,
            kind,
            rotational: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probes_a_directory_and_cleans_up() {
        let dir = std::env::temp_dir().join(format!("scratch-probe-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let limits = Limits {
            max_bytes: 8 * 1024 * 1024,
            max_time: Duration::from_millis(200),
            chunk_bytes: 1024 * 1024,
        };
        let probe = probe(&dir, &limits).unwrap();
        assert!(probe.write.bytes >= 1024 * 1024 && probe.write.bytes <= limits.max_bytes);
        // Reading stops at the time cap too, so it may not reach the end.
        assert!(probe.read.bytes >= 1024 * 1024 && probe.read.bytes <= probe.write.bytes);
        assert!(probe.slowest_bytes_per_second() > 0.0);
        assert!(probe.sustains(1.0, 1.0));
        assert!(!probe.sustains(f64::MAX, 1.0));
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 0, "probe file removed");
        fs::remove_dir(&dir).unwrap();
    }

    /// The buffered fallback, driven directly: which filesystems refuse
    /// direct I/O depends on the kernel (tmpfs accepts it since Linux 6.6).
    #[test]
    fn buffered_passes_measure_within_the_limits() {
        let dir = std::env::temp_dir().join(format!("scratch-probe-buf-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("probe");
        let limits = Limits {
            max_bytes: 8 * 1024 * 1024,
            max_time: Duration::from_millis(100),
            chunk_bytes: 1024 * 1024,
        };
        let mut buffer = AlignedBuffer::new(limits.chunk_bytes);
        fill_random(buffer.as_mut_slice());
        let write = timed_write(&path, &buffer, &limits, false).unwrap();
        assert!(write.bytes >= 1024 * 1024 && write.bytes <= limits.max_bytes);
        let read = timed_read(&path, &mut buffer, &limits, false).unwrap();
        assert!(read.bytes >= 1024 * 1024 && read.bytes <= write.bytes);
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn chunk_sizes_round_to_the_alignment() {
        let dir = std::env::temp_dir().join(format!("scratch-probe-odd-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let limits = Limits {
            max_bytes: 1,
            max_time: Duration::ZERO,
            chunk_bytes: 5000,
        };
        let probe = probe(&dir, &limits).unwrap();
        assert_eq!(probe.write.bytes, 8192);
        fs::remove_dir_all(&dir).unwrap();
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn classifies_memory_and_virtual_filesystems() {
        if Path::new("/dev/shm").is_dir() {
            let shm = classify(Path::new("/dev/shm"));
            assert_eq!(shm.kind, StorageKind::Memory);
            assert_eq!(shm.rotational, None);
        }
        assert!(mem_available().is_some_and(|bytes| bytes > 0));
    }

    #[cfg(target_os = "macos")]
    #[test]
    fn macos_bypasses_the_cache_on_local_disks() {
        let dir = std::env::temp_dir().join(format!("scratch-probe-mac-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let storage = classify(&dir);
        assert_eq!(storage.kind, StorageKind::Local, "{storage:?}");
        assert_eq!(storage.rotational, None);
        let limits = Limits {
            max_bytes: 4 * 1024 * 1024,
            max_time: Duration::from_millis(100),
            chunk_bytes: 1024 * 1024,
        };
        assert!(probe(&dir, &limits).unwrap().direct_io);
        fs::remove_dir_all(&dir).unwrap();
    }
}
