#!/usr/bin/env python3
import argparse
import collections
import gzip
import heapq
import itertools
import tempfile
from pathlib import Path


COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(seq):
    rc = reverse_complement(seq)
    return min(seq, rc)


def least_rotation_start(seq):
    n = len(seq)
    if n <= 1:
        return 0
    i, j, matched = 0, 1, 0
    while i < n and j < n and matched < n:
        left = seq[(i + matched) % n]
        right = seq[(j + matched) % n]
        if left == right:
            matched += 1
        elif left > right:
            i += matched + 1
            if i <= j:
                i = j + 1
            matched = 0
        else:
            j += matched + 1
            if j <= i:
                j = i + 1
            matched = 0
    return min(i, j)


def minimal_rotation(seq):
    start = least_rotation_start(seq)
    return seq[start:] + seq[:start]


def canonical_unitig(seq, k):
    if len(seq) > k and seq[: k - 1] == seq[-(k - 1) :]:
        body = seq[: len(seq) - k + 1]
        body = min(minimal_rotation(body), minimal_rotation(reverse_complement(body)))
        overlap = (body * ((k - 2) // len(body) + 1))[: k - 1]
        return body + overlap
    return canonical(seq)


def fasta_records(path):
    header = None
    seq = bytearray()
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as source:
        for raw in source:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if header is not None:
                    yield header, bytes(seq)
                header = line[1:]
                seq.clear()
            else:
                seq.extend(line.upper())
    if header is not None:
        yield header, bytes(seq)


def read_varint(source):
    value = 0
    shift = 0
    while True:
        byte = source.read(1)
        if not byte:
            raise ValueError("truncated color repository")
        value |= (byte[0] & 0x7f) << shift
        if byte[0] < 0x80:
            return value
        shift += 7
        if shift >= 35:
            raise ValueError("invalid color varint")


def load_rust_repository(path):
    path = Path(path)
    colors = {}
    with (path / "manifest.tsv").open() as manifest:
        next(manifest)
        for line in manifest:
            worker_text, records_text, file_text = line.rstrip("\n").split("\t")
            worker = int(worker_text)
            file_path = Path(file_text)
            if not file_path.is_absolute():
                file_path = path / file_path
            with file_path.open("rb") as source:
                for index in range(int(records_text)):
                    count = read_varint(source)
                    values = []
                    previous = 0
                    for _ in range(count):
                        previous += read_varint(source)
                        values.append(previous)
                    colors[worker | (index << 8)] = tuple(values)
    return colors


def rust_colored_records(fasta, repository, k):
    records = []
    for header, seq in fasta_records(fasta):
        fields = header.split()
        runs = []
        for raw_text in fields[1:]:
            raw = int(raw_text)
            offset = raw & 0xffffff
            coordinate = raw >> 24
            runs.append((offset, repository[coordinate]))
        vertex_count = len(seq) - k + 1
        expanded = []
        for index, (offset, color) in enumerate(runs):
            end = runs[index + 1][0] if index + 1 < len(runs) else vertex_count
            expanded.extend([color] * (end - offset))
        if len(expanded) != vertex_count:
            raise ValueError(f"malformed color runs for {seq!r}")
        rc = reverse_complement(seq)
        if rc < seq:
            seq = rc
            expanded.reverse()
        records.append((seq, tuple(expanded)))
    return sorted(records)


def topology(fasta, k):
    return sorted(canonical_unitig(seq, k) for _, seq in fasta_records(fasta))


def external_topology(fasta, k, directory, chunk_records):
    directory.mkdir(parents=True, exist_ok=True)
    chunks = []
    pending = []
    for _, seq in fasta_records(fasta):
        pending.append(canonical_unitig(seq, k))
        if len(pending) >= chunk_records:
            chunks.append(write_topology_chunk(pending, directory, len(chunks)))
            pending.clear()
    if pending:
        chunks.append(write_topology_chunk(pending, directory, len(chunks)))
    sources = [path.open("rb") for path in chunks]
    return heapq.merge(*sources), sources


def write_topology_chunk(records, directory, index):
    records.sort()
    path = directory / f"topology-{index:06}.txt"
    with path.open("wb") as output:
        for seq in records:
            output.write(seq)
            output.write(b"\n")
    return path


def compare_external_topology(cpp_fasta, rust_fasta, k, work_dir, chunk_records):
    with tempfile.TemporaryDirectory(dir=work_dir, prefix="cf3-topology-") as temp:
        temp = Path(temp)
        cpp, cpp_sources = external_topology(cpp_fasta, k, temp / "cpp", chunk_records)
        rust, rust_sources = external_topology(rust_fasta, k, temp / "rust", chunk_records)
        count = 0
        missing = object()
        try:
            for count, (cpp_seq, rust_seq) in enumerate(
                itertools.zip_longest(cpp, rust, fillvalue=missing), 1
            ):
                if cpp_seq is missing or rust_seq is missing:
                    raise ValueError("different topology counts")
                if cpp_seq != rust_seq:
                    raise SystemExit(
                        "FAIL: C++ and Rust topology multisets differ at record "
                        f"{count}: {cpp_seq.strip()!r} != {rust_seq.strip()!r}"
                    )
        except ValueError as error:
            raise SystemExit(
                f"FAIL: C++ and Rust topology counts differ after {count} records"
            ) from error
        finally:
            for source in cpp_sources + rust_sources:
                source.close()
    return count


def expected_colors(records, source_paths, k):
    wanted = set()
    for seq, _ in records:
        for offset in range(len(seq) - k + 1):
            wanted.add(canonical(seq[offset : offset + k]))
    observed = collections.defaultdict(list)
    for source_id, path in enumerate(source_paths, 1):
        found = set()
        for _, seq in fasta_records(path):
            for offset in range(len(seq) - k + 1):
                kmer = canonical(seq[offset : offset + k])
                if kmer in wanted:
                    found.add(kmer)
        for kmer in found:
            observed[kmer].append(source_id)
    expected = []
    for seq, _ in records:
        colors = tuple(
            tuple(observed[canonical(seq[offset : offset + k])])
            for offset in range(len(seq) - k + 1)
        )
        expected.append((seq, colors))
    return sorted(expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpp-fasta")
    parser.add_argument("--rust-fasta", required=True)
    parser.add_argument("--rust-repository")
    parser.add_argument("--source-list")
    parser.add_argument("--topology-only", action="store_true")
    parser.add_argument("--source-only", action="store_true")
    parser.add_argument("--external-work-dir")
    parser.add_argument("--chunk-records", type=int, default=1_000_000)
    parser.add_argument("-k", type=int, required=True)
    args = parser.parse_args()

    if not args.source_only and not args.cpp_fasta:
        parser.error("--cpp-fasta is required unless --source-only is used")
    if args.external_work_dir:
        if args.source_only or not args.topology_only:
            parser.error("external comparison currently requires --topology-only")
        Path(args.external_work_dir).mkdir(parents=True, exist_ok=True)
        count = compare_external_topology(
            args.cpp_fasta,
            args.rust_fasta,
            args.k,
            args.external_work_dir,
            args.chunk_records,
        )
        print(f"OK: {count} unitigs match C++ topology")
        return
    cpp_topology = None if args.source_only else topology(args.cpp_fasta, args.k)
    if args.topology_only:
        rust_topology = topology(args.rust_fasta, args.k)
    else:
        if not args.rust_repository or not args.source_list:
            parser.error("--rust-repository and --source-list are required for color validation")
        repository = load_rust_repository(args.rust_repository)
        rust = rust_colored_records(args.rust_fasta, repository, args.k)
        # Circular unitigs may differ by both strand and rotation. Keep color
        # records in their emitted order, but use the same circular
        # canonicalization as the topology-only comparison.
        rust_topology = topology(args.rust_fasta, args.k)
    if cpp_topology is not None and cpp_topology != rust_topology:
        cpp_only = collections.Counter(cpp_topology) - collections.Counter(rust_topology)
        rust_only = collections.Counter(rust_topology) - collections.Counter(cpp_topology)
        print(
            f"topology differences: {sum(cpp_only.values())} C++-only, "
            f"{sum(rust_only.values())} Rust-only"
        )
        for label, differences in (("C++", cpp_only), ("Rust", rust_only)):
            for seq, count in differences.most_common(3):
                print(f"{label}-only x{count}: {seq.decode()}" )
        raise SystemExit("FAIL: C++ and Rust topology multisets differ")
    if args.topology_only:
        print(f"OK: {len(rust_topology)} unitigs match C++ topology")
        return
    sources = [line.strip() for line in Path(args.source_list).read_text().splitlines() if line.strip()]
    expected = expected_colors(rust, sources, args.k)
    if rust != expected:
        for actual, truth in zip(rust, expected):
            if actual != truth:
                print("first color mismatch:", actual, truth)
                break
        raise SystemExit("FAIL: Rust colors differ from source-derived colors")
    if args.source_only:
        print(f"OK: {len(rust)} colored unitigs have source-derived colors")
    else:
        print(
            f"OK: {len(rust)} colored unitigs match C++ topology "
            "and source-derived colors"
        )


if __name__ == "__main__":
    main()
