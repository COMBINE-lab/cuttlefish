#!/usr/bin/env bash
# Run one cuttlefish3-rs (or C++ cuttlefish) build, capture wall time, peak RSS,
# peak open-descriptor count, per-phase timers, and output counts, and append one
# row to a TSV.
#
# Usage:
#   scripts/bench_variant.sh --variant NAME --list FILE --mode ref|read [options]
#
# Options:
#   --variant NAME     label for this binary/configuration (required)
#   --list FILE        input file list                     (required)
#   --mode ref|read    graph input class                   (required)
#   --bin PATH         binary (default target/release/cuttlefish3-rs)
#   --impl rust|cpp    which CLI dialect to emit           (default rust)
#   --color            build a colored graph
#   --threads N        default 256
#   -k N               default 31
#   --min-len N        default 12 for ref, 15 for read
#   --cutoff N         passed through when set
#   --ulimit-n N       set RLIMIT_NOFILE for the run
#   --mem-limit SIZE   run under a systemd scope with MemoryMax=SIZE (e.g. 64G).
#                      Caps page cache so bucket traffic must reach the device,
#                      which is the only way to see I/O volume on a host whose
#                      RAM would otherwise absorb it.
#   --root DIR         bench root (default /scratch4/rob/cf3-bench)
#   --tsv FILE         results TSV (default $ROOT/results.tsv)
#   --keep             keep the work dir and output FASTA
#   --extra "ARGS"     extra arguments appended verbatim
set -u -o pipefail

ROOT=/scratch4/rob/cf3-bench
BIN=target/release/cuttlefish3-rs
IMPL=rust
THREADS=256
K=31
MINLEN=
CUTOFF=
COLOR=0
KEEP=0
ULIMIT_N=
MEM_LIMIT=
EXTRA=
VARIANT=
LIST=
MODE=
TSV=

while [[ $# -gt 0 ]]; do
  case "$1" in
    --variant) VARIANT=$2; shift 2;;
    --list) LIST=$2; shift 2;;
    --mode) MODE=$2; shift 2;;
    --bin) BIN=$2; shift 2;;
    --impl) IMPL=$2; shift 2;;
    --color) COLOR=1; shift;;
    --threads) THREADS=$2; shift 2;;
    -k|--kmer-len) K=$2; shift 2;;
    --min-len) MINLEN=$2; shift 2;;
    --cutoff) CUTOFF=$2; shift 2;;
    --ulimit-n) ULIMIT_N=$2; shift 2;;
    --mem-limit) MEM_LIMIT=$2; shift 2;;
    --root) ROOT=$2; shift 2;;
    --tsv) TSV=$2; shift 2;;
    --keep) KEEP=1; shift;;
    --extra) EXTRA=$2; shift 2;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done

[[ -n $VARIANT && -n $LIST && -n $MODE ]] || { echo "need --variant, --list, --mode" >&2; exit 2; }
[[ -r $LIST ]] || { echo "cannot read list: $LIST" >&2; exit 2; }
: "${TSV:=$ROOT/results.tsv}"
if [[ -z $MINLEN ]]; then [[ $MODE == read ]] && MINLEN=15 || MINLEN=12; fi

DATASET=$(basename "$LIST" .list)
COLORTAG=$([[ $COLOR == 1 ]] && echo colored || echo uncolored)
TAG="${VARIANT}.${DATASET}.${COLORTAG}.t${THREADS}.k${K}"
WORK="$ROOT/work/$TAG"
OUTPREFIX="$ROOT/out/$TAG"
LOG="$ROOT/out/$TAG.log"

rm -rf "$WORK" "$OUTPREFIX".fa "$OUTPREFIX".fa.fa
mkdir -p "$WORK" "$ROOT/out"

if [[ $IMPL == rust ]]; then
  ARGS=(build "--$MODE" --list "$LIST" --kmer-len "$K" --min-len "$MINLEN"
        --threads "$THREADS" --work-dir "$WORK" --output "$OUTPREFIX")
  [[ $COLOR == 1 ]] && ARGS+=(--color)
  [[ -n $CUTOFF ]] && ARGS+=(--cutoff "$CUTOFF")
else
  ARGS=(build "--$MODE" -l "$LIST" -k "$K" --min-len "$MINLEN"
        -w "$WORK" -o "$OUTPREFIX")
  [[ $COLOR == 1 ]] && ARGS+=(--color)
  [[ -n $CUTOFF ]] && ARGS+=(-c "$CUTOFF")
  export PARLAY_NUM_THREADS=$THREADS
fi
[[ -n $EXTRA ]] && read -r -a EXTRA_ARR <<< "$EXTRA" && ARGS+=("${EXTRA_ARR[@]}")

echo "== $TAG" | tee "$LOG"
echo "== $(date -Is)  loadavg $(cut -d' ' -f1-3 /proc/loadavg)" | tee -a "$LOG"
{ printf '== cmd:'; printf ' %q' "$BIN" "${ARGS[@]}"; echo; } | tee -a "$LOG"

# Peak-descriptor sampler: polls the child's /proc/<pid>/fd while it runs.
FDPEAK_FILE=$(mktemp "$ROOT/out/.fdpeak.XXXXXX")
echo 0 > "$FDPEAK_FILE"
# Last-seen /proc/<pid>/io counters. They are monotonic per process and vanish
# when it exits, so the sampler keeps the highest it saw rather than reading
# once at the end. `wchar` is the volume the algorithm asks for; `write_bytes`
# is what actually reached the block layer, and on a large-memory host the two
# diverge sharply because page cache absorbs the difference.
IOFILE=$(mktemp "$ROOT/out/.io.XXXXXX")
echo "0 0 0 0" > "$IOFILE"

TIMEFILE=$(mktemp "$ROOT/out/.time.XXXXXX")
[[ -n $ULIMIT_N ]] && ulimit -n "$ULIMIT_N"

# Device counters bracket the run; they include unrelated system activity, so
# they only ever corroborate the per-process numbers.
# /sys/block is keyed by kernel name, so a device-mapper path has to be resolved
# through its symlink target (/dev/mapper/vg-lv -> ../dm-N) first.
DEV=$(lsblk -no KNAME "$(df --output=source "$WORK" | tail -1)" 2>/dev/null | head -1)
devstat() { cat "/sys/block/$1/stat" 2>/dev/null || echo; }
DEVSTAT_BEFORE=$(devstat "$DEV")

if [[ -n $MEM_LIMIT ]]; then
  RUNNER=(systemd-run --user --scope -q -p "MemoryMax=$MEM_LIMIT" -p "MemorySwapMax=0")
else
  RUNNER=()
fi

CF3_RS_PROFILE_RSS=1 "${RUNNER[@]}" /usr/bin/time -v -o "$TIMEFILE" "$BIN" "${ARGS[@]}" >>"$LOG" 2>&1 &
RUNPID=$!

(
  peak=0
  rchar=0; wchar=0; rbytes=0; wbytes=0
  # /usr/bin/time forks the real binary; watch the whole process subtree.
  while kill -0 $RUNPID 2>/dev/null; do
    n=0
    for p in $(pgrep -P $RUNPID 2>/dev/null) $RUNPID; do
      c=$(ls /proc/$p/fd 2>/dev/null | wc -l)
      n=$((n + c))
      while read -r key value; do
        case $key in
          rchar:)       (( value > rchar ))  && rchar=$value ;;
          wchar:)       (( value > wchar ))  && wchar=$value ;;
          read_bytes:)  (( value > rbytes )) && rbytes=$value ;;
          write_bytes:) (( value > wbytes )) && wbytes=$value ;;
        esac
      done < "/proc/$p/io" 2>/dev/null
    done
    (( n > peak )) && { peak=$n; echo "$peak" > "$FDPEAK_FILE"; }
    echo "$rchar $wchar $rbytes $wbytes" > "$IOFILE"
    sleep 0.2
  done
) &
SAMPLER=$!

wait $RUNPID; STATUS=$?
wait $SAMPLER 2>/dev/null

DEVSTAT_AFTER=$(devstat "$DEV")

WALL=$(awk -F': ' '/Elapsed \(wall clock\)/ {print $2}' "$TIMEFILE")
RSS=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$TIMEFILE")
FDPEAK=$(cat "$FDPEAK_FILE")
read -r RCHAR WCHAR RBYTES WBYTES < "$IOFILE"
# Fields 3 and 7 of /sys/block/<dev>/stat are sectors read and written.
DEVREAD=0; DEVWRITE=0
if [[ -n $DEVSTAT_BEFORE && -n $DEVSTAT_AFTER ]]; then
  DEVREAD=$(( ( $(awk '{print $3}' <<<"$DEVSTAT_AFTER") - $(awk '{print $3}' <<<"$DEVSTAT_BEFORE") ) * 512 ))
  DEVWRITE=$(( ( $(awk '{print $7}' <<<"$DEVSTAT_AFTER") - $(awk '{print $7}' <<<"$DEVSTAT_BEFORE") ) * 512 ))
fi

# Phase timers and counts from the tool's own stderr.
PART=$(grep -oP 'partition and bucket emission completed in \K[0-9.]+' "$LOG" | tail -1)
BUILD=$(grep -oP 'graph build completed in \K[0-9.]+' "$LOG" | tail -1)
LOCAL=$(grep -oP 'local contraction phase completed in \K[0-9.]+' "$LOG" | tail -1)
CONTRACT=$(grep -oP 'discontinuity graph contracted in \K[0-9.]+' "$LOG" | tail -1)
EXPAND=$(grep -oP 'contracted graph expanded in \K[0-9.]+' "$LOG" | tail -1)
PARTWORKERS=$(grep -oP '(?:(?:uncolored|colored) partition using |partition streaming [0-9]+ reader\(s\) into )\K[0-9]+' "$LOG" | tail -1)
UNITIGS=$(grep -oP 'wrote \K[0-9]+(?= (?:colored )?unitig)' "$LOG" | tail -1)
BASES=$(grep -oP 'unitig\(s\), \K[0-9]+(?= base)' "$LOG" | tail -1)
# The C++ CLI reports neither total; count them from the FASTA it wrote.
if [[ -z ${UNITIGS:-} ]]; then
  FA="$OUTPREFIX.fa"; [[ -s $FA ]] || FA="$OUTPREFIX.fa.fa"
  if [[ -s $FA ]]; then
    read -r UNITIGS BASES < <(awk '/^>/{u++;next}{b+=length($0)}END{print u, b}' "$FA")
  fi
fi
WORKBYTES=$(du -sb "$WORK" 2>/dev/null | cut -f1)
FASTA=$(stat -c %s "$OUTPREFIX".fa 2>/dev/null || stat -c %s "$OUTPREFIX".fa.fa 2>/dev/null || echo 0)

if [[ ! -s $TSV ]]; then
  printf 'variant\timpl\tdataset\tmode\tcolor\tthreads\tk\tmin_len\tstatus\twall\tmax_rss_kb\tfd_peak\tpart_workers\tpartition_s\tbuild_s\tlocal_s\tcontract_s\texpand_s\tunitigs\tbases\twork_bytes\tfasta_bytes\trchar\twchar\tread_bytes\twrite_bytes\tdev_read\tdev_write\tmem_limit\ttag\n' > "$TSV"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$VARIANT" "$IMPL" "$DATASET" "$MODE" "$COLORTAG" "$THREADS" "$K" "$MINLEN" "$STATUS" \
  "$WALL" "$RSS" "$FDPEAK" "${PARTWORKERS:-}" "${PART:-}" "${BUILD:-}" "${LOCAL:-}" "${CONTRACT:-}" "${EXPAND:-}" \
  "${UNITIGS:-}" "${BASES:-}" "${WORKBYTES:-}" "$FASTA" \
  "$RCHAR" "$WCHAR" "$RBYTES" "$WBYTES" "$DEVREAD" "$DEVWRITE" "${MEM_LIMIT:-none}" "$TAG" >> "$TSV"

rm -f "$TIMEFILE" "$FDPEAK_FILE" "$IOFILE"
if [[ $KEEP == 0 ]]; then rm -rf "$WORK"; fi

echo "== status=$STATUS wall=$WALL rss=${RSS}KiB fd_peak=$FDPEAK unitigs=${UNITIGS:-?} bases=${BASES:-?}"
echo "== io: wchar=$(numfmt --to=iec ${WCHAR:-0}) write_bytes=$(numfmt --to=iec ${WBYTES:-0}) rchar=$(numfmt --to=iec ${RCHAR:-0}) read_bytes=$(numfmt --to=iec ${RBYTES:-0}) dev_write=$(numfmt --to=iec ${DEVWRITE:-0}) mem_limit=${MEM_LIMIT:-none}"
exit $STATUS
