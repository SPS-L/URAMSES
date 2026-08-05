#!/usr/bin/env bash
# Replace a module kit directory with the contents of a RAMSES kit zip.
#
#   refresh_kit.sh <kit-zip> <mod-dir> <expected-lib> <expected-tag>
#
# Validates that the zip carries the expected library, at least one .mod, and
# a BUILDINFO.txt whose ramses_tag matches <expected-tag>. Only then does it
# empty <mod-dir> and unpack into it.
#
# Emptying first is the point: it evicts modules a RAMSES release has dropped.
# sparse_matrix_mod.mod became sparse_matrix_optimized_mod.mod between v3.40
# and v3.50, and unpacking over the top would have left the stale file behind
# for the next build to find.
#
# Validation happens entirely in a staging directory, so a failure of any kind
# leaves <mod-dir> untouched.
set -u

if [ $# -ne 4 ]; then
    echo "usage: $0 <kit-zip> <mod-dir> <expected-lib> <expected-tag>" >&2
    exit 2
fi
ZIP="$1"
MOD_DIR="$2"
LIB="$3"
TAG="$4"

[ -f "$ZIP" ]     || { echo "FAIL: kit zip not found: $ZIP"; exit 1; }
[ -d "$MOD_DIR" ] || { echo "FAIL: module directory not found: $MOD_DIR"; exit 1; }

# A failed mktemp leaves STAGE empty, and `unzip -d ""` silently unpacks into
# the caller's cwd instead of a staging dir -- which, run from the repo root
# as this script normally is, would scatter kit files across the working
# tree. Fail loudly instead of letting that happen.
STAGE="$(mktemp -d)" || { echo "FAIL: could not create a staging directory (mktemp -d failed)"; exit 1; }
[ -n "$STAGE" ] && [ -d "$STAGE" ] || {
    echo "FAIL: could not create a staging directory (mktemp -d returned no usable path)"
    exit 1
}
trap 'rm -rf "$STAGE"' EXIT

unzip -q -j "$ZIP" -d "$STAGE" || { echo "FAIL: could not unzip $ZIP"; exit 1; }

[ -f "$STAGE/$LIB" ] || { echo "FAIL: $(basename "$ZIP") contains no $LIB"; exit 1; }

MOD_COUNT="$(ls -1 "$STAGE"/*.mod 2>/dev/null | wc -l | tr -d ' ')"
[ "$MOD_COUNT" -gt 0 ] || { echo "FAIL: $(basename "$ZIP") contains no .mod files"; exit 1; }

[ -f "$STAGE/BUILDINFO.txt" ] || {
    echo "FAIL: $(basename "$ZIP") contains no BUILDINFO.txt"
    echo "The kit predates toolchain recording, or the RAMSES release workflow did not run it."
    exit 1
}

ZIP_TAG="$(sed -n 's/^ramses_tag[[:space:]]\{1,\}//p' "$STAGE/BUILDINFO.txt" | head -n1)"
if [ "$ZIP_TAG" != "$TAG" ]; then
    echo "FAIL: $(basename "$ZIP") BUILDINFO says ramses_tag '$ZIP_TAG', expected '$TAG'"
    exit 1
fi

# The release notes publish mod_abi_version straight out of BUILDINFO, so it
# has to agree with every module actually in the zip -- otherwise a release
# would advertise the wrong compiler requirement. A one-file sample can miss
# a skew that only shows up further down the list, so walk every .mod in a
# single pass, reading its banner (they are gzipped text) as we go. Stays
# lenient when a banner is unreadable, or when BUILDINFO makes no claim at
# all, matching how the Makefiles treat an undeterminable version -- but a
# banner that *can* be read and disagrees is never lenient.
CLAIMED_ABI="$(sed -n 's/^mod_abi_version[[:space:]]\{1,\}//p' "$STAGE/BUILDINFO.txt" | head -n1)"
ACTUAL_ABI=""
for MOD in "$STAGE"/*.mod; do
    [ -e "$MOD" ] || continue
    MOD_ABI="$(gzip -dc "$MOD" 2>/dev/null | head -c 256 \
        | sed -n "s/.*module version '\([0-9]*\)'.*/\1/p")"
    [ -n "$MOD_ABI" ] || continue
    [ -n "$ACTUAL_ABI" ] || ACTUAL_ABI="$MOD_ABI"
    if [ -n "$CLAIMED_ABI" ] && [ "$MOD_ABI" != "$CLAIMED_ABI" ]; then
        echo "FAIL: $(basename "$ZIP") BUILDINFO claims .mod ABI $CLAIMED_ABI but $(basename "$MOD") is ABI $MOD_ABI"
        exit 1
    fi
done

# ${MOD_DIR:?} refuses to expand if the variable is somehow empty, so a bug
# upstream cannot turn this into rm -rf /*. The brace pattern below empties
# dotfiles too (but never "." or ".."), so a stray .DS_Store or editor swap
# file cannot survive a refresh; rm never follows a symlink argument into its
# target, so a symlink inside <mod-dir> is unlinked rather than descended into.
rm -rf "${MOD_DIR:?}"/{*,.[!.]*,..?*}
cp "$STAGE"/* "$MOD_DIR"/

echo "OK: $MOD_DIR refreshed from $(basename "$ZIP") ($MOD_COUNT .mod files, ABI ${ACTUAL_ABI:-unknown}, $LIB, RAMSES $TAG)"
