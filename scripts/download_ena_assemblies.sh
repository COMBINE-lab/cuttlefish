#!/usr/bin/env bash
# Download the first N assemblies from the ENA 661k bacterial assembly index.

set -euo pipefail

readonly INDEX_URL="https://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt"
readonly URL_PREFIX="https://ftp.ebi.ac.uk"
readonly PATH_PREFIX="/ebi/ftp/pub/databases/ENA2018-bacteria-661k/Assemblies/"

usage() {
    cat <<'EOF'
Usage: download_ena_assemblies.sh N OUTPUT_DIR

Download the first N genomes in the ENA 661k bacterial assembly index into
OUTPUT_DIR. The directory is created when necessary. Existing partial files
are resumed, and completed files are retained.

Environment:
  ENA_DOWNLOAD_JOBS  concurrent file downloads (default: 32)

Outputs:
  *.contigs.fa.gz                 downloaded assemblies
  genomes.list                    absolute paths for the selected assemblies
  .sampleid_assembly_paths.txt    cached ENA source index
EOF
}

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 2
fi

readonly genome_count="$1"
readonly requested_output_dir="$2"
readonly jobs="${ENA_DOWNLOAD_JOBS:-32}"

if [[ ! "$genome_count" =~ ^[1-9][0-9]*$ ]]; then
    printf 'error: N must be a positive integer, got %q\n' "$genome_count" >&2
    exit 2
fi
if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'error: ENA_DOWNLOAD_JOBS must be a positive integer, got %q\n' "$jobs" >&2
    exit 2
fi
if ! command -v aria2c >/dev/null 2>&1; then
    printf 'error: aria2c is required; install the aria2 package\n' >&2
    exit 127
fi

mkdir -p "$requested_output_dir"
readonly output_dir="$(cd "$requested_output_dir" && pwd -P)"
readonly index_path="$output_dir/.sampleid_assembly_paths.txt"
download_input="$(mktemp "$output_dir/.ena-download.XXXXXX")"
list_tmp="$(mktemp "$output_dir/.genomes-list.XXXXXX")"
cleanup() {
    rm -f "$download_input" "$list_tmp"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

printf 'Refreshing ENA assembly index...\n' >&2
aria2c \
    --allow-overwrite=true \
    --auto-file-renaming=false \
    --conditional-get=true \
    --console-log-level=warn \
    --dir="$output_dir" \
    --download-result=hide \
    --file-allocation=none \
    --max-tries=5 \
    --out="$(basename "$index_path")" \
    --remote-time=true \
    --retry-wait=2 \
    --summary-interval=0 \
    "$INDEX_URL"

# aria2's input format places per-file options on indented lines following a
# URL. Keeping all transfers in one process allows connection and DNS reuse.
awk \
    -v count="$genome_count" \
    -v path_prefix="$PATH_PREFIX" \
    -v url_prefix="$URL_PREFIX" \
    '
    BEGIN { FS = "\t" }
    NR <= count {
        accession = $1
        path = $2
        if (accession == "" || index(path, path_prefix) != 1) {
            printf "invalid ENA index entry at line %d\n", NR > "/dev/stderr"
            invalid = 1
            exit 4
        }
        parts = split(path, component, "/")
        filename = component[parts]
        expected = accession ".contigs.fa.gz"
        if (filename != expected) {
            printf "unexpected assembly filename at line %d: %s\n", NR, filename > "/dev/stderr"
            invalid = 1
            exit 4
        }
        sub(/^\/ebi\/ftp/, "", path)
        print url_prefix path
        print "  out=" filename
    }
    NR == count { exit }
    END {
        if (!invalid && NR < count) {
            printf "ENA index contains only %d entries; requested %d\n", NR, count > "/dev/stderr"
            exit 3
        }
    }
    ' "$index_path" >"$download_input"

printf 'Downloading %s assemblies with %s concurrent transfers into %s...\n' \
    "$genome_count" "$jobs" "$output_dir" >&2
aria2c \
    --allow-overwrite=true \
    --auto-file-renaming=false \
    --connect-timeout=20 \
    --console-log-level=warn \
    --continue=true \
    --dir="$output_dir" \
    --disk-cache=64M \
    --download-result=hide \
    --enable-http-keep-alive=true \
    --enable-http-pipelining=true \
    --file-allocation=none \
    --input-file="$download_input" \
    --max-concurrent-downloads="$jobs" \
    --max-connection-per-server=1 \
    --max-file-not-found=3 \
    --max-tries=8 \
    --min-split-size=20M \
    --remote-time=true \
    --retry-wait=2 \
    --split=1 \
    --summary-interval=10 \
    --timeout=60

awk \
    -v count="$genome_count" \
    -v output_dir="$output_dir" \
    '
    BEGIN { FS = "\t" }
    NR <= count {
        parts = split($2, component, "/")
        print output_dir "/" component[parts]
    }
    NR == count { exit }
    ' "$index_path" >"$list_tmp"
mv -f "$list_tmp" "$output_dir/genomes.list"

printf 'Downloaded %s assemblies. Cuttlefish input list: %s/genomes.list\n' \
    "$genome_count" "$output_dir" >&2
