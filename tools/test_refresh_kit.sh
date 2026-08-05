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
kit_dir          modules_l
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

# --- zip where an early .mod agrees but a later one doesn't --------------
# A one-file sample would miss this: alpha.mod matches the claimed ABI,
# zeta.mod (later alphabetically) does not.
mkdir -p src_skew_multi
printf 'archive' > src_skew_multi/libramses.a
fakemod 15 alpha > src_skew_multi/alpha.mod   # agrees with claimed ABI 15
fakemod 16 zeta  > src_skew_multi/zeta.mod    # disagrees
buildinfo v9.99 > src_skew_multi/BUILDINFO.txt
( cd src_skew_multi && zip -q -j "$TMPD/skew_multi.zip" . -r )

# --- zip missing BUILDINFO ----------------------------------------------
mkdir -p src_nobi
printf 'archive' > src_nobi/libramses.a
fakemod 15 alpha > src_nobi/alpha.mod
( cd src_nobi && zip -q -j "$TMPD/nobi.zip" . -r )

# --- zip whose BUILDINFO is missing a required key ------------------------
# A valid ramses_tag is not enough: the other twelve keys are cells in the
# published release table, and a silently missing one ships a blank column.
mkdir -p src_missingkey
printf 'archive' > src_missingkey/libramses.a
fakemod 15 alpha > src_missingkey/alpha.mod
buildinfo v9.99 | grep -v '^mod_abi_version' > src_missingkey/BUILDINFO.txt
( cd src_missingkey && zip -q -j "$TMPD/missingkey.zip" . -r )

# --- zip whose BUILDINFO names the wrong kit_dir --------------------------
# kit_dir is written by RAMSES and printed into the release table, but was
# never checked against the directory it is actually being installed into.
mkdir -p src_wrongkitdir
printf 'archive' > src_wrongkitdir/libramses.a
fakemod 15 alpha > src_wrongkitdir/alpha.mod
buildinfo v9.99 | sed 's/^kit_dir .*/kit_dir          modules_m/' > src_wrongkitdir/BUILDINFO.txt
( cd src_wrongkitdir && zip -q -j "$TMPD/wrongkitdir.zip" . -r )

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
    rm -rf "$TMPD/modules_l"
    mkdir -p "$TMPD/modules_l"
    printf 'old' > "$TMPD/modules_l/libramses.a"
    fakemod 15 alpha             > "$TMPD/modules_l/alpha.mod"
    fakemod 15 sparse_matrix_mod > "$TMPD/modules_l/sparse_matrix_mod.mod"
}

# --- usage ---------------------------------------------------------------
"$SCRIPT" >/dev/null 2>&1
[ $? -eq 2 ] && ok "no args exits 2" || fail "no args should exit 2"

# --- happy path ----------------------------------------------------------
fresh_target
OUT="$("$SCRIPT" "$TMPD/good.zip" "$TMPD/modules_l" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 0 ] && ok "valid kit exits 0" || fail "valid kit exited $RC"
[ -f "$TMPD/modules_l/libramses.a" ]  && ok "library present"      || fail "library missing"
[ -f "$TMPD/modules_l/beta.mod" ]     && ok "new module unpacked"  || fail "beta.mod missing"
[ -f "$TMPD/modules_l/BUILDINFO.txt" ] && ok "BUILDINFO unpacked"  || fail "BUILDINFO missing"
[ "$(cat "$TMPD/modules_l/libramses.a")" = "archive" ] \
    && ok "library overwritten with new content" || fail "library not replaced"

# The reason this script exists rather than unzip-over-the-top:
[ ! -f "$TMPD/modules_l/sparse_matrix_mod.mod" ] \
    && ok "stale module evicted" || fail "stale sparse_matrix_mod.mod survived"

# --- dotfiles evicted alongside everything else ---------------------------
# rm -rf DIR/* alone skips names beginning with "." -- a stray .DS_Store or
# editor swap file would then survive every future refresh forever.
fresh_target
mkdir -p "$TMPD/modules_l/.stale_hidden_dir"
printf 'ds_store' > "$TMPD/modules_l/.stale_hidden"
fakemod 15 leftover > "$TMPD/modules_l/.stale_hidden_dir/leftover.mod"
"$SCRIPT" "$TMPD/good.zip" "$TMPD/modules_l" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 0 ] && ok "refresh with planted dotfiles exits 0" \
    || fail "refresh with planted dotfiles should exit 0"
[ ! -e "$TMPD/modules_l/.stale_hidden" ] \
    && ok "hidden file evicted" || fail "hidden file .stale_hidden survived"
[ ! -e "$TMPD/modules_l/.stale_hidden_dir" ] \
    && ok "hidden directory evicted" || fail "hidden directory .stale_hidden_dir survived"

# --- tag mismatch --------------------------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
"$SCRIPT" "$TMPD/good.zip" "$TMPD/modules_l" libramses.a v1.00 >/dev/null 2>&1
[ $? -eq 1 ] && ok "tag mismatch exits 1" || fail "tag mismatch should exit 1"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "tag mismatch left the kit untouched" \
    || fail "kit was modified despite tag mismatch: '$BEFORE' -> '$AFTER'"

# --- BUILDINFO disagreeing with the shipped modules ----------------------
# The release table publishes mod_abi_version straight from BUILDINFO, so a
# kit whose metadata contradicts its own .mod files must never be installed.
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
OUT="$("$SCRIPT" "$TMPD/skew.zip" "$TMPD/modules_l" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "ABI skew exits 1" || fail "ABI skew should exit 1"
echo "$OUT" | grep -q "claims .mod ABI 15" && ok "names the claimed ABI" || fail "claimed ABI missing"
echo "$OUT" | grep -q "ABI 16"             && ok "names the actual ABI"  || fail "actual ABI missing"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "ABI skew left the kit untouched" || fail "kit was modified"

# --- ABI skew hiding behind an early agreeing .mod ------------------------
# Sampling only the first .mod (alpha.mod, which agrees) would miss zeta.mod
# disagreeing further down the list.
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
OUT="$("$SCRIPT" "$TMPD/skew_multi.zip" "$TMPD/modules_l" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "multi-mod ABI skew exits 1" || fail "multi-mod ABI skew should exit 1"
echo "$OUT" | grep -q "zeta.mod" && ok "names the offending module" || fail "offending module name missing"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "multi-mod ABI skew left the kit untouched" || fail "kit was modified"

# --- missing BUILDINFO ---------------------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
"$SCRIPT" "$TMPD/nobi.zip" "$TMPD/modules_l" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing BUILDINFO exits 1" || fail "should exit 1"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "missing BUILDINFO left the kit untouched" \
    || fail "kit was modified"

# --- BUILDINFO missing a required key -------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
OUT="$("$SCRIPT" "$TMPD/missingkey.zip" "$TMPD/modules_l" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "missing key exits 1" || fail "missing key should exit 1"
echo "$OUT" | grep -q "mod_abi_version" && ok "names the missing key" || fail "missing key not named"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "missing key left the kit untouched" || fail "kit was modified"

# --- BUILDINFO kit_dir does not match the target directory ----------------
# Feeding one platform's zip at another platform's directory must be caught,
# not silently installed.
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
OUT="$("$SCRIPT" "$TMPD/wrongkitdir.zip" "$TMPD/modules_l" libramses.a v9.99 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "wrong kit_dir exits 1" || fail "wrong kit_dir should exit 1"
echo "$OUT" | grep -q "modules_m" && ok "names the BUILDINFO's kit_dir" || fail "BUILDINFO kit_dir not named"
echo "$OUT" | grep -q "modules_l" && ok "names the target directory" || fail "target directory not named"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "wrong kit_dir left the kit untouched" || fail "kit was modified"

# --- missing library -----------------------------------------------------
fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
"$SCRIPT" "$TMPD/nolib.zip" "$TMPD/modules_l" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing library exits 1" || fail "should exit 1"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "missing library left the kit untouched" || fail "kit was modified"

# --- no .mod files -------------------------------------------------------
fresh_target
"$SCRIPT" "$TMPD/nomod.zip" "$TMPD/modules_l" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "no .mod files exits 1" || fail "should exit 1"

# --- absent inputs -------------------------------------------------------
"$SCRIPT" "$TMPD/nope.zip" "$TMPD/modules_l" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing zip exits 1" || fail "missing zip should exit 1"

"$SCRIPT" "$TMPD/good.zip" "$TMPD/no-such-dir" libramses.a v9.99 >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing target dir exits 1" || fail "missing dir should exit 1"

# --- mktemp -d failure must abort, not silently unpack into cwd ----------
# A failed mktemp leaves STAGE empty, and `unzip -d ""` unpacks into the
# caller's cwd instead. We can't make the real mktemp fail cleanly on
# demand, so shadow it on PATH with one that always fails.
FAKEBIN="$TMPD/fakebin"
mkdir -p "$FAKEBIN"
cat > "$FAKEBIN/mktemp" <<'EOF'
#!/usr/bin/env bash
echo "mktemp: simulated failure" >&2
exit 1
EOF
chmod +x "$FAKEBIN/mktemp"

fresh_target
BEFORE="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
mkdir -p "$TMPD/cwd_probe"
(
    cd "$TMPD/cwd_probe" || exit 1
    PATH="$FAKEBIN:$PATH" "$SCRIPT" "$TMPD/good.zip" "$TMPD/modules_l" libramses.a v9.99 >/dev/null 2>&1
)
RC=$?
[ $RC -ne 0 ] && ok "mktemp failure exits non-zero" || fail "mktemp failure should exit non-zero"
LEFTOVER="$(ls -A "$TMPD/cwd_probe" 2>/dev/null)"
[ -z "$LEFTOVER" ] && ok "mktemp failure did not scatter files into cwd" \
    || fail "mktemp failure left files in cwd: $LEFTOVER"
AFTER="$(ls "$TMPD/modules_l" | sort | tr '\n' ' ')"
[ "$BEFORE" = "$AFTER" ] && ok "mktemp failure left the kit untouched" || fail "kit was modified"

echo ""
if [ "$FAILURES" -eq 0 ]; then
    echo "PASS: all assertions held"
    exit 0
fi
echo "$FAILURES failure(s)"
exit 1
