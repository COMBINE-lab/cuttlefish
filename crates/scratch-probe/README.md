# scratch-probe

How fast does a directory's storage absorb and return a stream, and what kind
of storage is it?

Programs with large temporary files can trade CPU for I/O, for example by
compressing what they spill. Whether that pays depends on the storage under
the work directory. `scratch-probe` answers cheaply enough to ask at startup:

- `classify(dir)` uses metadata only (no I/O). It reports the filesystem;
  whether it is local, networked (NFS, SMB, Lustre, GPFS, CephFS, BeeGFS, ...)
  or memory-backed; and whether the block device under it rotates.
  Device-mapper and md stacks are followed down to their disks.
- `probe(dir, &Limits::default())` also times a sequential write (with
  `fdatasync`) and read-back of random bytes.
  - It is capped at 256 MiB and about half a second, so fast storage costs a
    fraction of a second and slow storage stops early.
  - On Linux it uses `O_DIRECT`, so it measures the device rather than the
    page cache. It falls back to buffered I/O where direct I/O is refused.
  - The probe file is removed before returning.

```rust
let probe = scratch_probe::probe(work_dir, &scratch_probe::Limits::default())?;
// Compress spills unless the storage keeps up with our stream rate, twice over.
let compress = !probe.sustains(expected_bytes_per_second, 2.0);
```

The crate reports; policy (what rate a program produces, and what margin is
safe) stays with the caller. A probe is a snapshot of idle storage: shared or
network storage can be busier later.

Classification is Linux-specific. Elsewhere `classify` reports `Unknown`, and
`probe` times buffered I/O.
