use super::*;

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "cf3-offset-{}-{:?}",
            std::process::id(),
            std::thread::current().id(),
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn record(offset: u64, reverse: bool, is_cycle: bool) -> MaterializedStitchedCoordRecord {
    MaterializedStitchedCoordRecord {
        path_id: (1 << 61) + 17,
        rank: 23,
        label_offset: offset,
        label_len: 31,
        reverse,
        is_cycle,
        color_index: u32::MAX,
        color_count: 0,
    }
}

#[test]
fn materialized_offsets_round_trip_across_u32_and_u46_boundaries() {
    let path = Path::new("offset-round-trip");
    for offset in [
        0,
        u32::MAX as u64,
        1 << 32,
        (1 << 32) + 41,
        4_561_650_435,
        LoadedMaterializedStitchedCoordRecord::MAX_LABEL_OFFSET,
    ] {
        for flags in 0..4 {
            for color_index in [17, u32::MAX] {
                let mut input = record(offset, flags & 1 != 0, flags & 2 != 0);
                input.color_index = color_index;
                input.color_count = u32::from(color_index != u32::MAX);
                let encoded = encoded_materialized_stitched_coord_record(input, path).unwrap();
                assert_eq!(encoded.len(), 24);
                assert_eq!(
                    decoded_materialized_stitched_coord_record(&encoded, path).unwrap(),
                    input
                );
                // Production loads the native POD layout rather than the test decoder.
                let loaded = unsafe {
                    std::ptr::read_unaligned(
                        encoded
                            .as_ptr()
                            .cast::<LoadedMaterializedStitchedCoordRecord>(),
                    )
                };
                assert_eq!(loaded.label_offset(), offset);
                assert_eq!(loaded.reverse(), input.reverse);
                assert_eq!(loaded.is_cycle(), input.is_cycle);
                assert_eq!(loaded.path_id, input.path_id);
                assert_eq!(loaded.color_start, color_index);
                assert_eq!(loaded.color_count(), input.color_count);
            }
        }
    }
    for offset in [1 << 46, u64::MAX] {
        let error = encoded_materialized_stitched_coord_record(record(offset, false, false), path)
            .unwrap_err()
            .to_string();
        assert!(error.contains("offset-round-trip"));
        assert!(error.contains(&offset.to_string()));
        assert!(error.contains("46-bit maximum"));
    }
}

#[test]
fn materialized_offset_rebasing_carries_without_changing_flags() {
    let path = Path::new("offset-rebase");
    for flags in 0..4 {
        let mut loaded = LoadedMaterializedStitchedCoordRecord::new(
            17,
            23,
            u32::MAX - 2,
            31,
            flags & 1 != 0,
            flags & 2 != 0,
            u32::MAX,
            0,
        );
        loaded.rebase_label_offset(7, path).unwrap();
        assert_eq!(loaded.label_offset(), (1 << 32) + 4);
        loaded.rebase_label_offset(1 << 32, path).unwrap();
        assert_eq!(loaded.label_offset(), (2 << 32) + 4);
        assert_eq!(loaded.reverse(), flags & 1 != 0);
        assert_eq!(loaded.is_cycle(), flags & 2 != 0);
        assert_eq!(loaded.color_start, u32::MAX);
        let before = loaded;
        assert!(loaded.rebase_label_offset(1 << 46, path).is_err());
        assert_eq!(loaded, before);
        assert!(loaded.rebase_label_offset(u64::MAX, path).is_err());
        assert_eq!(loaded, before);
    }
    assert_eq!(
        checked_materialized_label_bytes((1 << 46) - 1, 1, path).unwrap(),
        1 << 46
    );
    assert!(checked_materialized_label_bytes(1 << 46, 1, path).is_err());
    assert!(checked_materialized_label_bytes(u64::MAX, 1, path).is_err());
}

fn coord(rank: u64) -> StitchedCoordRecord {
    StitchedCoordRecord {
        path_id: 17,
        rank,
        unitig_index: 0,
        reverse: false,
        is_cycle: false,
    }
}

fn tail(rank: u64, label: &[u8], colors: &[UnitigColor]) -> PendingMaterializedBucket {
    PendingMaterializedBucket {
        records: vec![LoadedMaterializedStitchedCoordRecord::new(
            17,
            rank,
            0,
            label.len() as u32,
            false,
            false,
            if colors.is_empty() { u32::MAX } else { 0 },
            colors.len() as u32,
        )],
        labels: label.to_vec(),
        colors: colors.to_vec(),
    }
}

#[test]
fn materialized_loader_rejects_old_formats_and_out_of_bounds_labels() {
    let directory = TestDirectory::new();
    let mut writer = MaterializedStitchedCoordShardWriter::create(&directory.0, 0, 7).unwrap();
    writer.write_record(&coord(1), b"AACCGG").unwrap();
    let entry = writer.finish().unwrap();
    let file = OpenOptions::new()
        .write(true)
        .open(&entry.coord_path)
        .unwrap();
    file.write_all_at(b"CF3MCB2\0", 0).unwrap();
    assert!(read_materialized_stitched_coord_bucket_file(&entry).is_err());
    file.write_all_at(MATERIALIZED_STITCH_COORD_MAGIC, 0)
        .unwrap();
    assert!(read_materialized_stitched_coord_bucket_file(&entry).is_ok());
    // A well-formed offset whose label would extend beyond the file.
    file.write_all_at(&1u32.to_le_bytes(), STITCH_COORD_HEADER_LEN + 8)
        .unwrap();
    assert!(read_materialized_stitched_coord_bucket_file(&entry).is_err());
}

#[test]
fn materialized_failed_tail_keeps_input_files() {
    let directory = TestDirectory::new();
    let mut writer = MaterializedStitchedCoordShardWriter::create(&directory.0, 0, 7).unwrap();
    writer.write_record(&coord(1), b"AACCGG").unwrap();
    let entry = writer.finish().unwrap();
    let mut invalid = tail(2, b"GGTTAA", &[]);
    invalid.records[0]
        .set_label_offset(1 << 32, &entry.label_path)
        .unwrap();
    let error = load_materialized_stitched_coord_bucket_file_group_with_tails(
        std::slice::from_ref(&entry),
        &[invalid],
    )
    .err()
    .unwrap()
    .to_string();
    assert!(error.contains(entry.label_path.to_str().unwrap()));
    assert!(entry.coord_path.exists());
    assert!(entry.label_path.exists());
}

#[test]
fn materialized_writer_checks_limit_before_appending_labels() {
    let directory = TestDirectory::new();
    let mut writer = MaterializedStitchedCoordShardWriter::create(&directory.0, 0, 7).unwrap();
    writer.label_bytes = (1 << 46) - 2;
    assert!(writer.write_record(&coord(1), b"AACCGG").is_err());
    assert!(
        writer
            .write_pending_batch(&mut tail(2, b"GGTTAA", &[]))
            .is_err()
    );
    assert_eq!(writer.records, 0);
    assert!(writer.record_buffer.is_empty());
    writer.label_out.as_mut().unwrap().flush().unwrap();
    assert_eq!(fs::metadata(&writer.label_path).unwrap().len(), 0);
}

// This exercises the real writer, native bulk loader, multiple-shard and
// retained-tail rebasing, and assembly with offsets above 4 GiB. Sparse file
// padding keeps disk use tiny, but each iteration needs approximately 4 GiB
// of RAM. Run explicitly in release as part of scalability validation.
#[test]
#[ignore = "allocates a 4 GiB label buffer; run with --ignored --test-threads=1"]
fn materialized_large_bucket_round_trip() {
    let directory = TestDirectory::new();
    for colored in [false, true] {
        let colors = |coordinate| {
            if colored {
                // Neighboring fragments agree on their shared endpoint k-mer.
                vec![
                    UnitigColor::new(0, crate::state::ColorCoordinate::from_u40(coordinate)),
                    UnitigColor::new(3, crate::state::ColorCoordinate::from_u40(coordinate + 1)),
                ]
            } else {
                Vec::new()
            }
        };
        let mut writer = MaterializedStitchedCoordShardWriter::create(&directory.0, 0, 7).unwrap();
        let base = (1u64 << 32) - 3;
        let labels = writer.label_out.as_mut().unwrap();
        labels.get_mut().set_len(base).unwrap();
        labels.seek(SeekFrom::Start(base)).unwrap();
        writer.label_bytes = base;
        if colored {
            writer
                .write_colored_record(&coord(1), b"AACCGG", &colors(42))
                .unwrap();
        } else {
            writer.write_record(&coord(1), b"AACCGG").unwrap();
        }
        writer
            .write_pending_batch(&mut tail(2, b"CGGTTA", &colors(43)))
            .unwrap();
        let first = writer.finish().unwrap();
        let mut second = MaterializedStitchedCoordShardWriter::create(&directory.0, 1, 7).unwrap();
        second
            .write_pending_batch(&mut tail(3, b"TTAACC", &colors(44)))
            .unwrap();
        let second = second.finish().unwrap();
        let mut loaded = load_materialized_stitched_coord_bucket_file_group_with_tails(
            &[first, second],
            &[tail(4, b"ACCGAT", &colors(45))],
        )
        .unwrap();
        assert_eq!(loaded.labels.len() as u64, base + 24);
        assert_eq!(
            loaded
                .records
                .iter()
                .map(|r| r.label_offset())
                .collect::<Vec<_>>(),
            vec![base, base + 6, base + 12, base + 18]
        );
        let unitigs = reduce_materialized_stitched_coord_bucket::<3>(
            &mut loaded.records,
            &loaded.labels,
            &loaded.colors,
        );
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].label, b"AACCGGTTAACCGAT");
        let expected_colors = if colored {
            (0..5)
                .map(|i| {
                    UnitigColor::new(
                        i * 3,
                        crate::state::ColorCoordinate::from_u40(42 + u64::from(i)),
                    )
                })
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        assert_eq!(unitigs[0].colors, expected_colors);
    }
}

// Also copied unchanged into the 3.0.1 snapshot for an interleaved comparison
// of the production writer/load/sort/assembly path with u32-sized buffers.
#[test]
#[ignore = "manual materialization performance benchmark"]
fn benchmark_materialized_coordinate_round_trip() {
    let directory = std::env::temp_dir().join(format!("cf3-offset-bench-{}", std::process::id()));
    fs::create_dir_all(&directory).unwrap();
    const RECORDS: u64 = 2_000_000;
    let label = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT";
    for colored in [false, true] {
        let colors = [UnitigColor::new(
            0,
            crate::state::ColorCoordinate::from_u40(42),
        )];
        let shared = SharedMaterializedWriters::new(&directory, 1, 1);
        let mut batch = SharedMaterializedBatch::new(&shared, 1);
        let started = Instant::now();
        for index in 0..RECORDS {
            let record = StitchedCoordRecord {
                path_id: index.wrapping_mul(0x9e3779b97f4a7c15),
                rank: 1,
                unitig_index: 0,
                reverse: false,
                is_cycle: false,
            };
            if colored {
                batch
                    .write_materialized_colored_record(0, &record, label, &colors)
                    .unwrap();
            } else {
                batch.write_materialized_record(0, &record, label).unwrap();
            }
        }
        batch.finish().unwrap();
        let buckets = shared.finish().unwrap();
        let write = started.elapsed();
        let started = Instant::now();
        let mut loaded = load_materialized_stitched_coord_bucket_file_group_with_tails(
            &buckets.manifest,
            &buckets.retained[0],
        )
        .unwrap();
        let load = started.elapsed();
        let started = Instant::now();
        loaded.records.sort_unstable_by_key(|r| (r.path_id, r.rank));
        let sort = started.elapsed();
        let started = Instant::now();
        let mut bases = 0;
        let mut count = 0;
        reduce_sorted_materialized_stitched_coord_bucket_with::<31, _>(
            &mut loaded.records,
            &loaded.labels,
            &loaded.colors,
            |label, colors| {
                bases += std::hint::black_box(label).len() as u64;
                std::hint::black_box(colors);
                count += 1;
            },
        );
        let assemble = started.elapsed();
        assert_eq!(count, RECORDS);
        assert_eq!(bases, RECORDS * label.len() as u64);
        println!(
            "MCOORD_BENCH colored={colored} records={RECORDS} write={:.6} load={:.6} sort={:.6} assemble={:.6}",
            write.as_secs_f64(),
            load.as_secs_f64(),
            sort.as_secs_f64(),
            assemble.as_secs_f64()
        );
    }
    fs::remove_dir_all(directory).unwrap();
}
