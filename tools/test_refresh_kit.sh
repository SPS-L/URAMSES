#!/usr/bin/env bash
# Local test for tools/refresh_kit.sh.
#
#   test_refresh_kit.sh
#
# Builds synthetic kit zips and asserts refresh_kit.sh replaces a kit
# directory wholesale, evicts files the new kit dropped, and leaves the
# target untouched whenever validation fails.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/tools/refresh_kit.sh"
FAILURES=0

fail() { echo "FAIL: $*"; FAILURES=$((FAILURES + 1)); }
ok()   { echo "ok: $*"; }

TMPD="$(mktemp -d)"
trap 'rm -rf "$TMPD"' EXIT
cd "$TMPD" || exit 1

buildinfo() {   # buildinfo <tag> > file
    cat <<EOF
kit_dir          modules_lin
ramses_tag       $1
ramses_commit    deadbeef
built_utc        2026-09-01T10:12:44Z
runner_image     ubuntu-22.04
compiler         GNU Fortran (Test) 11.4.0
target           x86_64-pc-linux-gnu
mod_abi_version  15
fflags           -O2
ldflags          -fopenmp
blas             OpenBLAS 0.3.20
runtime_floor    glibc 2.35
consume_with     Makefile.linux
EOF
}

# A .mod file is gzipped text opening with a version banner; fabricating one
# lets these tests assert on ABI handling without invoking a compiler.
fakemod() { printf "GFORTRAN module version '%s' created from %s.f90\n" "$1" "$2" | gzip; }

# --- a well-formed kit zip: lib + two mods + BUILDINFO -------------------
mkdir -p src_good
printf 'archive' > src_good/libramses.a
fakemod 15 alpha > src_good/alpha.mod
fakemod 15 beta  > src_good/beta.mod
buildinfo v9.99 > src_good/BUILDINFO.txt
( cd src_good && zip -q -j "$TMPD/good.zip" . -r )

# --- zip whose BUILDINFO disagrees with its own modules ------------------
mkdir -p src_skew
printf 'archive' > src_skew/libramses.a
fakemod 16 alpha > src_skew/alpha.mod
buildinfo v9.99 > src_skew/BUILDINFO.txt   # buildinfo() claims ABI 15
( cd src_skew && zip -q -j "$TMPD/skew.zip" . -r )

# --- zip missing BUILDINFO ----------------------------------------------
mkdir -p src_nobi
printf 'archive' > src_nobi/libramses.a
fakemod 15 alpha > src_nobi/alpha.mod
( cd src_nobi && zip -q -j "$TMPD/nobi.zip" . -r )

# --- zip missing the library --------------------------------------------
mkdir -p src_nolib
fakemod 15 alpha > src_nolib/alpha.mod
buildinfo v9.99 > src_nolib/BUILDINFO.txt
( cd src_nolib && zip -q -j "$TMPD/nolib.zip" . -r )

# --- zip with no .mod files ---------------------------------------------
mkdir -p src_nomod
printf 'archive' > src_nomod/libramses.a
buildinfo v9.99 > src_nomod/BUILDINFO.txt
( cd src_nomod && zip -q -j "$TMPD/nomod.zip" . -r )

# A target kit dir holding a stale module the new kit no longer ships.
fresh_target() {
    rm -rf "$TMPD/modules_lin"
    mkdir -p "$TMPD/modules_lin"
    printf 'old' > "$TMPD/modules_lin/libramses.a"
    fakemod 15 alpha             > "$TMPD/modules_lin/alpha.mod"
    fakemod 15 sparse_matrix_mod > "$TMPD/modules_lin/sparse_matrix_mod.mod"
}

# --- usage ---------------------------------------------------------------
"$SCRIPT" >/dev/null 2>&1
[ $? -eq 2 ] && ok "no args exits 2" || fail "no args should exit 2"

# --- happy path ----------------------------------------------------------
fresh_target
OUT="$("$SCRIPT" "$TMPD/good.zip" "$TMPD/modules_lin" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 0 ] && ok "valid kit exits 0" || fail "valid kit exited $RC"
[ -f "$TMPD/modules_lin/libramses.a" ]  && ok "library present"      || fail "library missing"
[ -f "$TMPD/modules_lin/beta.mod" ]     && ok "new module unpacked"  || fail "beta.mod missing"
[ -f "$TMPD/modules_lin/BUILDINFO.txt" ] && ok "BUILDINFO unpacked"  || fail "BUILDINFO missing"
[ "$(cat "$TMPD/modules_lin/libramses.a")" = "archive" ] \
    && ok "library overwritten with new content" || fail "library not replaced"

# The reason this script exists rather than unzip-over-the-top:
[ ! -f "$TMPD/modules_lin/sparse_matrix_mod.mod" ] \
    && ok "stale module evicted" || fail "stale sparse_matrix_mod.mod survived"

# --- tag mismatch --------------------------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
"$SCRIPT" "$TMPD/good.zip" "$TMPD/modules_lin" libramses.a v1.00 >/dev/null 2>&1
[ $? -eq 1 ] && ok "tag mismatch exits 1" || fail "tag mismatch should exit 1"
AFTER="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "tag mismatch left the kit untouched" \
    || fail "kit was modified despite tag mismatch: '$BEFORE' -> '$AFTER'"

# --- BUILDINFO disagreeing with the shipped modules ----------------------
# The release table publishes mod_abi_version straight from BUILDINFO, so a
# kit whose metadata contradicts its own .mod files must never be installed.
fresh_target
BEFORE="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
OUT="$("$SCRIPT" "$TMPD/skew.zip" "$TMPD/modules_lin" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "ABI skew exits 1" || fail "ABI skew should exit 1"
echo "$OUT" | grep -q "claims .mod ABI 15" && ok "names the claimed ABI" || fail "claimed ABI missing"
echo "$OUT" | grep -q "ABI 16"             && ok "names the actual ABI"  || fail "actual ABI missing"
AFTER="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "ABI skew left the kit untouched" || fail "kit was modified"

# --- missing BUILDINFO ---------------------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
"$SCRIPT" "$TMPD/nobi.zip" "$TMPD/modules_lin" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing BUILDINFO exits 1" || fail "should exit 1"
AFTER="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "missing BUILDINFO left the kit untouched" \
    || fail "kit was modified"

# --- missing library -----------------------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
"$SCRIPT" "$TMPD/nolib.zip" "$TMPD/modules_lin" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing library exits 1" || fail "should exit 1"
AFTER="$(ls "$TMPD/modules_lin" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "missing library left the kit untouched" || fail "kit was modified"

# --- no .mod files -------------------------------------------------------
fresh_target
"$SCRIPT" "$TMPD/nomod.zip" "$TMPD/modules_lin" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "no .mod files exits 1" || fail "should exit 1"

# --- absent inputs -------------------------------------------------------
"$SCRIPT" "$TMPD/nope.zip" "$TMPD/modules_lin" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing zip exits 1" || fail "missing zip should exit 1"

"$SCRIPT" "$TMPD/good.zip" "$TMPD/no-such-dir" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing target dir exits 1" || fail "missing dir should exit 1"

echo ""
if [ "$FAILURES" -eq 0 ]; then
    echo "PASS: all assertions held"
    exit 0
fi
echo "$FAILURES failure(s)"
exit 1
