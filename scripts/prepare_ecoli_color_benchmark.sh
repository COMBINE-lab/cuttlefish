#!/usr/bin/env bash
set -euo pipefail

root=${1:-data/ecoli-refseq-150}
accessions="$root/accessions.txt"
archive="$root/ncbi_dataset.zip"
package="$root/package"
genomes="$root/genomes.list"

if ! command -v datasets >/dev/null 2>&1; then
    module load ncbi_tools/22.0.0
fi

test -s "$accessions"
mkdir -p "$root"

if [[ ! -s "$archive" ]]; then
    datasets download genome accession \
        --inputfile "$accessions" \
        --include genome \
        --filename "$archive" \
        --no-progressbar
fi

if [[ ! -d "$package/ncbi_dataset/data" ]]; then
    rm -rf "$package"
    unzip -q "$archive" -d "$package"
fi

find "$package/ncbi_dataset/data" -type f -name '*_genomic.fna' \
    | LC_ALL=C sort > "$genomes"

expected=$(wc -l < "$accessions")
observed=$(wc -l < "$genomes")
if [[ "$observed" -ne "$expected" ]]; then
    printf 'expected %d genomes, found %d\n' "$expected" "$observed" >&2
    exit 1
fi

printf 'prepared %d genome colors under %s\n' "$observed" "$package"
