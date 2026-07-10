#!/usr/bin/env python3
"""Generate deterministic small uncolored FASTA fixtures."""

from pathlib import Path


RECORDS = [
    (
        "linear_repeat",
        "ACGTTGCATGTCAGTACGTTGCATGTCAGTA",
    ),
    (
        "branch_left",
        "GATTACAGATTACCGTAAGGCTA",
    ),
    (
        "branch_right",
        "GATTACAGATTACCGTCCGGCTA",
    ),
    (
        "cycle_like",
        "ATCGATCGATCGATCGATCGATC",
    ),
    (
        "repeat_path",
        "AACCGGTTAACCGGTTAACCGGTT",
    ),
    (
        "n_split",
        "AACCGGTTAANNTTACCGGTTAA",
    ),
    (
        "cutoff_copy_1",
        "TTAACCGGTTAACTTAACCGGTTAA",
    ),
    (
        "cutoff_copy_2",
        "TTAACCGGTTAACTTAACCGGTTAA",
    ),
]


def main():
    repo = Path(__file__).resolve().parents[1]
    out = repo / "data" / "generated" / "uncolored-synthetic.fa"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        for name, seq in RECORDS:
            handle.write(f">{name}\n{seq}\n")
    print(out)


if __name__ == "__main__":
    main()
