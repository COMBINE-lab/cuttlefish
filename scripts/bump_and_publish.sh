#!/usr/bin/env bash
set -euo pipefail

# Bump the shared workspace version and publish the cuttlefish-rs crates to
# crates.io in dependency order. Modeled on salmon's bump_and_publish.sh,
# adapted for this workspace: two published crates whose version lives in a
# single [workspace.package] block.
#
# `cuttlefish-rs-compat-tests` is `publish = false` and never appears here.

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./scripts/bump_and_publish.sh <version> [--publish] [--dry-run] [--skip-tests]

Bumps [workspace.package] version (both published crates inherit it), commits,
tags, pushes, then optionally publishes to crates.io in dependency order.

Options:
  --publish     Publish to crates.io after bumping, committing, tagging, pushing
  --dry-run     Show what would be done; do not modify tracked files, commit,
                tag, push, or publish. Runs `cargo publish --dry-run` per crate
                so packaging is validated without uploading.
  --skip-tests  Skip the `cargo test --workspace` preflight. The fixture tests
                are this project's actual correctness gate, so skip them only
                when you have just run them yourself.
  -h, --help    Show this help message

Publish order (dependency-topological):
  cuttlefish-rs -> cuttlefish-rs-cli
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false
SKIP_TESTS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish) PUBLISH=true ;;
        --dry-run) DRY_RUN=true ;;
        --skip-tests) SKIP_TESTS=true ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1" ;;
        *)
            [[ -z "$VERSION" ]] || die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || { usage; exit 1; }

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$WORKSPACE_ROOT"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
TAG="v${VERSION}"

# Dependency-topological publish order. crates.io will not accept a crate until
# the versions of its path+version deps are already in the index; `cargo
# publish` blocks (wait-for-publish) until each appears, so a sequential loop is
# enough.
CRATES=(
    cuttlefish-rs
    cuttlefish-rs-cli
)

MANIFEST_BACKUP=""
LOCKFILE_BACKUP=""
MANIFEST_UPDATED=false
COMMIT_CREATED=false

cleanup() {
    local status=$?
    if [[ "$status" -ne 0 && "$DRY_RUN" == false && "$MANIFEST_UPDATED" == true && "$COMMIT_CREATED" == false ]]; then
        [[ -n "$MANIFEST_BACKUP" && -f "$MANIFEST_BACKUP" ]] && cp "$MANIFEST_BACKUP" "$ROOT_CARGO"
        [[ -n "$LOCKFILE_BACKUP" && -f "$LOCKFILE_BACKUP" ]] && cp "$LOCKFILE_BACKUP" "$LOCKFILE"
        echo "restored $ROOT_CARGO and $LOCKFILE after failure" >&2
    fi
    [[ -n "$MANIFEST_BACKUP" && -f "$MANIFEST_BACKUP" ]] && rm -f "$MANIFEST_BACKUP"
    [[ -n "$LOCKFILE_BACKUP" && -f "$LOCKFILE_BACKUP" ]] && rm -f "$LOCKFILE_BACKUP"
    return "$status"
}
trap cleanup EXIT

[[ -f "$ROOT_CARGO" ]] || die "not found: $ROOT_CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE (run 'cargo build' once to generate it)"

# Current workspace version lives under [workspace.package].
CURRENT_VERSION="$(sed -n '/^\[workspace.package\]/,/^\[/{s/^version = "\(.*\)"/\1/p}' "$ROOT_CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine [workspace.package] version from $ROOT_CARGO"

[[ "$CURRENT_VERSION" != "$VERSION" ]] || die "workspace version is already $VERSION"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi
if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi
git remote get-url origin >/dev/null 2>&1 || die "git remote 'origin' is not configured"

echo "Current workspace version : $CURRENT_VERSION"
echo "New workspace version     : $VERSION"
echo "Tag                       : $TAG"
echo "Publish                   : $([[ "$PUBLISH" == true ]] && echo yes || echo no)"
echo "Dry-run                   : $([[ "$DRY_RUN" == true ]] && echo yes || echo no)"
echo "Crates (in order)         : ${CRATES[*]}"
echo

echo "Preflight: cargo check"
cargo check -q

# The fixture tests assert exact unitig and base counts. They are the check that
# distinguishes "it compiles" from "it builds the right graph", and a release is
# exactly when that distinction matters.
if [[ "$SKIP_TESTS" == false ]]; then
    echo "Preflight: cargo test --workspace"
    cargo test --workspace -q
else
    echo "Preflight: skipping tests (--skip-tests)"
fi

echo "Updating [workspace.package] version: $CURRENT_VERSION -> $VERSION"
if [[ "$DRY_RUN" == false ]]; then
    MANIFEST_BACKUP="$(mktemp "${TMPDIR:-/tmp}/cuttlefish-Cargo.toml.XXXXXX")"
    LOCKFILE_BACKUP="$(mktemp "${TMPDIR:-/tmp}/cuttlefish-Cargo.lock.XXXXXX")"
    cp "$ROOT_CARGO" "$MANIFEST_BACKUP"
    cp "$LOCKFILE" "$LOCKFILE_BACKUP"

    # Bump the version line within the [workspace.package] block, and the
    # version="X" requirement on the internal cuttlefish-rs path dep in
    # [workspace.dependencies] (so the published CLI resolves to the new
    # library version rather than the previous release).
    sed -i.bak "/^\[workspace.package\]/,/^\[/{s/^version = \".*\"/version = \"${VERSION}\"/}" "$ROOT_CARGO"
    sed -i.bak "/^\[workspace.dependencies\]/,/^\[profile/{s|^\(cuttlefish-rs = { path = \"crates/cuttlefish-rs\", version = \"\).*\(\" }\)|\1${VERSION}\2|}" "$ROOT_CARGO"
    rm -f "${ROOT_CARGO}.bak"

    MANIFEST_UPDATED=true

    # Refresh Cargo.lock for the new versions.
    cargo check -q

    UPDATED_VERSION="$(sed -n '/^\[workspace.package\]/,/^\[/{s/^version = "\(.*\)"/\1/p}' "$ROOT_CARGO" | head -1)"
    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "workspace version update failed"

    DEP_VERSION="$(sed -n 's|^cuttlefish-rs = { path = "crates/cuttlefish-rs", version = "\(.*\)" }|\1|p' "$ROOT_CARGO" | head -1)"
    [[ "$DEP_VERSION" == "$VERSION" ]] || die "[workspace.dependencies] cuttlefish-rs version update failed"
else
    echo "Dry-run: would rewrite [workspace.package] version and the cuttlefish-rs dep version in $ROOT_CARGO"
fi

if [[ "$DRY_RUN" == true ]]; then
    # Per-crate packaging validation. Only meaningful in --dry-run mode, where
    # the version is NOT bumped, so cuttlefish-rs-cli's library requirement
    # still resolves against the already-published version in the index. In
    # --publish mode the version IS bumped, so the CLI's `cuttlefish-rs =
    # "^<new>"` cannot resolve until the library is actually published; the real
    # publish loop below handles that ordering via cargo's wait-for-publish.
    #
    # Until the library has been published at the version the CLI requires,
    # the CLI's dependency cannot resolve against the index -- whether the
    # index holds nothing at all or only older versions (e.g. the 0.0.x name
    # placeholders). Both are facts about the index and not about this
    # workspace, and the real publish loop below resolves them by publishing
    # the library first; so a CLI resolution failure of that shape is reported
    # and tolerated rather than treated as a validation failure.
    echo
    echo "Per-crate package validation (cargo publish --dry-run, in order)"
    validation_failed=false
    for crate in "${CRATES[@]}"; do
        echo "--- $crate"
        status=0
        output=$(cargo publish -p "$crate" --dry-run --allow-dirty 2>&1) || status=$?
        printf '%s\n' "$output"
        if [[ $status -ne 0 ]]; then
            if [[ "$crate" != "cuttlefish-rs" ]] && \
               grep -q 'failed to select a version for the requirement `cuttlefish-rs' <<<"$output"; then
                echo ":: $crate cannot be validated until cuttlefish-rs is published" >&2
                echo "::    at the required version; the index cannot resolve its" >&2
                echo "::    dependency yet. Not treated as a failure." >&2
            else
                validation_failed=true
                echo "::  $crate failed packaging validation" >&2
            fi
        fi
    done
    [[ "$validation_failed" == false ]] || die "packaging validation failed"
    echo
    echo "Dry-run complete (no commit, tag, push, or publish performed)"
    exit 0
fi

run git add "$ROOT_CARGO" "$LOCKFILE"
run git commit -m "chore(release): bump cuttlefish workspace to v${VERSION}"
COMMIT_CREATED=true

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$PUBLISH" == true ]]; then
    for crate in "${CRATES[@]}"; do
        echo "Publishing $crate ..."
        run cargo publish -p "$crate"
    done
    echo "All cuttlefish-rs crates published for v${VERSION}"
else
    echo "Skipping crates.io publish; re-run with --publish to publish v${VERSION}"
fi

echo
echo "Release bump complete for v${VERSION}"
echo "The pushed tag triggers .github/workflows/release.yml (dist), which"
echo "builds the platform binaries and creates the GitHub release itself,"
echo "with the CHANGELOG.md section for ${VERSION} as its body."
