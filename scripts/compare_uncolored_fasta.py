#!/usr/bin/env python3
"""Compare uncolored Cuttlefish C++ and Rust FASTA outputs.

The comparator ignores record identifiers, output ordering, and unitig strand by
canonicalizing each sequence against its reverse complement.
"""

import argparse
import collections
import os
import resource
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run C++ and Rust uncolored builds and diff normalized FASTA unitigs."
    )
    parser.add_argument("--cpp-bin", required=True, help="path to the C++ cuttlefish binary")
    parser.add_argument(
        "--rust-bin",
        help="path to cuttlefish3-rs; defaults to `cargo run -p cuttlefish-rs-cli --`",
    )
    parser.add_argument(
        "-s",
        "--seq",
        action="append",
        required=True,
        help="input FASTA/FASTQ; repeat or use comma-separated values for multiple inputs",
    )
    parser.add_argument("-k", "--kmer-len", type=int, required=True, help="odd k-mer length")
    parser.add_argument("--min-len", type=int, required=True, help="minimizer length")
    parser.add_argument("-c", "--cutoff", type=int, help="frequency cutoff")
    parser.add_argument(
        "--read",
        action="store_true",
        help="build a read-input graph instead of reference-input graph",
    )
    parser.add_argument(
        "--work-root",
        default=None,
        help="directory for temporary build outputs; defaults to a temporary directory",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="keep temporary output directories after the comparison",
    )
    parser.add_argument(
        "--threads",
        default="1",
        help="PARLAY_NUM_THREADS value for the C++ run; default: 1",
    )
    parser.add_argument(
        "--ulimit-n",
        type=int,
        default=32768,
        help="soft file descriptor limit to request before running builds; default: 32768",
    )
    parser.add_argument(
        "--max-diff",
        type=int,
        default=20,
        help="maximum missing/extra sequences to print",
    )
    parser.add_argument(
        "--metrics-tsv",
        help="optional path to append benchmark metrics for the C++ and Rust runs",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    raise_nofile_limit(args.ulimit_n)
    root = Path(args.work_root) if args.work_root else Path(tempfile.mkdtemp(prefix="cf3-compare-"))
    root.mkdir(parents=True, exist_ok=True)
    cxx_work = root / "cxx-work"
    rust_work = root / "rust-work"
    cxx_prefix = root / "cxx"
    rust_prefix = root / "rust"
    cxx_work.mkdir(parents=True, exist_ok=True)
    rust_work.mkdir(parents=True, exist_ok=True)

    try:
        seqs = expand_seq_args(args.seq)
        input_mode = "--read" if args.read else "--ref"
        cxx_cmd = [
            args.cpp_bin,
            "build",
            input_mode,
            "-k",
            str(args.kmer_len),
            "--min-len",
            str(args.min_len),
            "-w",
            str(cxx_work),
            "-o",
            str(cxx_prefix),
        ]
        add_seq_args(cxx_cmd, seqs)
        if args.cutoff is not None:
            cxx_cmd.extend(["-c", str(args.cutoff)])

        rust_cmd = rust_command(args, rust_work, rust_prefix, seqs, input_mode)

        cxx_env = os.environ.copy()
        cxx_env["PARLAY_NUM_THREADS"] = args.threads
        cxx_run = run(cxx_cmd, cxx_env)
        rust_run = run(rust_cmd, os.environ.copy())

        cxx_fasta = Path(f"{cxx_prefix}.fa")
        rust_fasta = Path(f"{rust_prefix}.fa")
        cxx = normalized_fasta(cxx_fasta)
        rust = normalized_fasta(rust_fasta)

        print_summary("C++", cxx_fasta, cxx)
        print_run_metrics("C++", cxx_run)
        print_summary("Rust", rust_fasta, rust)
        print_run_metrics("Rust", rust_run)
        if args.metrics_tsv:
            write_metrics_tsv(
                Path(args.metrics_tsv),
                [
                    metrics_row("C++", cxx_cmd, cxx_run, cxx, cxx_fasta, cxx_work),
                    metrics_row("Rust", rust_cmd, rust_run, rust, rust_fasta, rust_work),
                ],
            )

        if cxx == rust:
            print("OK: normalized FASTA multisets match")
            return 0

        print("FAIL: normalized FASTA multisets differ")
        print_diff("missing from Rust", cxx - rust, args.max_diff)
        print_diff("extra in Rust", rust - cxx, args.max_diff)
        return 1
    finally:
        if not args.keep and args.work_root is None:
            shutil.rmtree(root, ignore_errors=True)
        elif args.keep:
            print(f"kept comparison outputs under {root}")


def expand_seq_args(raw_values):
    seqs = []
    for raw_value in raw_values:
        seqs.extend(value for value in raw_value.split(",") if value)
    return seqs


def add_seq_args(cmd, seqs):
    for seq in seqs:
        cmd.extend(["-s", seq])


def rust_command(args, rust_work, rust_prefix, seqs, input_mode):
    common = [
        "build",
        input_mode,
        "-k",
        str(args.kmer_len),
        "--min-len",
        str(args.min_len),
        "-t",
        str(args.threads),
        "-w",
        str(rust_work),
        "-o",
        str(rust_prefix),
    ]
    add_seq_args(common, seqs)
    if args.cutoff is not None:
        common.extend(["-c", str(args.cutoff)])
    if args.rust_bin:
        return [args.rust_bin, *common]
    return ["cargo", "run", "-q", "-p", "cuttlefish-rs-cli", "--", *common]


def run(cmd, env):
    print("+ " + " ".join(cmd))
    started = time.monotonic()
    proc = subprocess.Popen(cmd, env=env)
    _, status, usage = os.wait4(proc.pid, 0)
    elapsed = time.monotonic() - started
    proc.returncode = wait_status_to_returncode(status)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return {
        "elapsed_s": elapsed,
        "max_rss_kb": usage.ru_maxrss,
    }


def wait_status_to_returncode(status):
    if hasattr(os, "waitstatus_to_exitcode"):
        return os.waitstatus_to_exitcode(status)
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    if os.WIFSIGNALED(status):
        return -os.WTERMSIG(status)
    return 1


def raise_nofile_limit(limit):
    if limit is None:
        return

    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    requested = min(limit, hard)
    if soft < requested:
        resource.setrlimit(resource.RLIMIT_NOFILE, (requested, hard))


def normalized_fasta(path):
    records = collections.Counter()
    for seq in read_fasta(path):
        normalized = canonical(seq.upper())
        records[normalized] += 1
    return records


def read_fasta(path):
    records = []
    current = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    records.append("".join(current))
                    current.clear()
                continue
            current.append(line)
    if current:
        records.append("".join(current))
    return records


def canonical(seq):
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(seq, rc)


def print_summary(label, path, records):
    total_bases = sum(len(seq) * count for seq, count in records.items())
    total_records = sum(records.values())
    print(f"{label}: {total_records} unitig(s), {total_bases} base(s), {path}")


def print_run_metrics(label, run_result):
    print(
        f"{label}: elapsed {run_result['elapsed_s']:.3f}s, "
        f"max RSS {run_result['max_rss_kb']} KiB"
    )


def metrics_row(label, cmd, run_result, records, fasta, work_dir):
    total_bases = sum(len(seq) * count for seq, count in records.items())
    total_records = sum(records.values())
    return {
        "label": label,
        "elapsed_s": f"{run_result['elapsed_s']:.6f}",
        "max_rss_kb": str(run_result["max_rss_kb"]),
        "unitigs": str(total_records),
        "bases": str(total_bases),
        "work_bytes": str(directory_bytes(work_dir)),
        "fasta": str(fasta),
        "command": shell_join(cmd),
    }


def shell_join(cmd):
    return " ".join(shlex.quote(str(arg)) for arg in cmd)


def write_metrics_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "label",
        "elapsed_s",
        "max_rss_kb",
        "unitigs",
        "bases",
        "work_bytes",
        "fasta",
        "command",
    ]
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a") as handle:
        if write_header:
            handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(row[field] for field in fields) + "\n")


def directory_bytes(path):
    total = 0
    for root, _, files in os.walk(path):
        root = Path(root)
        for name in files:
            try:
                total += (root / name).stat().st_size
            except FileNotFoundError:
                pass
    return total


def print_diff(label, diff, max_diff):
    total = sum(diff.values())
    print(f"{label}: {total}")
    for idx, (seq, count) in enumerate(sorted(diff.items())):
        if idx >= max_diff:
            print(f"... {total - max_diff} more")
            break
        suffix = f" x{count}" if count != 1 else ""
        print(f"  {seq}{suffix}")


if __name__ == "__main__":
    sys.exit(main())
