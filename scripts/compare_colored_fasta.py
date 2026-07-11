#!/usr/bin/env python3
import argparse
import collections
from pathlib import Path


COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(seq):
    rc = reverse_complement(seq)
    return min(seq, rc)


def fasta_records(path):
    header = None
    seq = bytearray()
    with Path(path).open("rb") as source:
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


def topology(fasta):
    return sorted(canonical(seq) for _, seq in fasta_records(fasta))


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
    parser.add_argument("--cpp-fasta", required=True)
    parser.add_argument("--rust-fasta", required=True)
    parser.add_argument("--rust-repository")
    parser.add_argument("--source-list")
    parser.add_argument("--topology-only", action="store_true")
    parser.add_argument("-k", type=int, required=True)
    args = parser.parse_args()

    cpp_topology = topology(args.cpp_fasta)
    if args.topology_only:
        rust_topology = topology(args.rust_fasta)
    else:
        if not args.rust_repository or not args.source_list:
            parser.error("--rust-repository and --source-list are required for color validation")
        repository = load_rust_repository(args.rust_repository)
        rust = rust_colored_records(args.rust_fasta, repository, args.k)
        rust_topology = [seq for seq, _ in rust]
    if cpp_topology != rust_topology:
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
    print(f"OK: {len(rust)} colored unitigs match C++ topology and source-derived colors")


if __name__ == "__main__":
    main()
