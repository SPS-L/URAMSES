#!/usr/bin/env bash
# Local test for tools/check_kit.sh.
#
#   test_check_kit.sh
#
# Builds synthetic kit directories in a temp dir and asserts check_kit.sh
# accepts a matching kit, rejects a mismatched one, and degrades gracefully
# when the kit has no BUILDINFO.txt or the compiler cannot be probed.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/tools/check_kit.sh"
FAILURES=0

fail() { echo "FAIL: $*"; FAILURES=$((FAILURES + 1)); }
ok()   { echo "ok: $*"; }

TMPD="$(mktemp -d)"
trap 'rm -rf "$TMPD"' EXIT

# A .mod written by the local gfortran: matches by construction.
mkdir -p "$TMPD/good"
printf 'module probe_mod\nend module probe_mod\n' > "$TMPD/probe.f90"
gfortran -fsyntax-only -J"$TMPD/good" "$TMPD/probe.f90"
LOCAL_ABI="$(gzip -dc "$TMPD/good"/probe_mod.mod | head -c 256 \
    | sed -n "s/.*module version '\([0-9]*\)'.*/\1/p")"
echo "local gfortran emits .mod ABI $LOCAL_ABI"

cat > "$TMPD/good/BUILDINFO.txt" <<EOF
kit_dir          modules_l
ramses_tag       v9.99
ramses_commit    deadbeef
built_utc        2026-09-01T10:12:44Z
runner_image     ubuntu-24.04
compiler         GNU Fortran (Test) 13.3.0
target           x86_64-pc-linux-gnu
mod_abi_version  $LOCAL_ABI
fflags           -O2
ldflags          -fopenmp
blas             OpenBLAS 0.3.26
runtime_floor    glibc 2.39
consume_with     Makefile.linux
EOF

# A kit whose .mod claims an ABI no compiler emits.
mkdir -p "$TMPD/bad"
printf "GFORTRAN module version '99' created from FAKE.f90\n" | gzip > "$TMPD/bad/aaa.mod"
sed "s/^mod_abi_version.*/mod_abi_version  99/" "$TMPD/good/BUILDINFO.txt" > "$TMPD/bad/BUILDINFO.txt"

# A kit predating provenance tracking.
mkdir -p "$TMPD/nobi"
cp "$TMPD/good"/probe_mod.mod "$TMPD/nobi/"

# --- usage ---------------------------------------------------------------
"$SCRIPT" >/dev/null 2>&1
[ $? -eq 2 ] && ok "no args exits 2" || fail "no args should exit 2"

"$SCRIPT" "$TMPD/good" gfortran >/dev/null 2>&1
[ $? -eq 2 ] && ok "two args exits 2" || fail "two args should exit 2"

# --- matching kit --------------------------------------------------------
OUT="$("$SCRIPT" "$TMPD/good" gfortran linux 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 0 ] && ok "matching kit exits 0" || fail "matching kit exited $RC"
echo "$OUT" | grep -q "RAMSES v9.99"       && ok "prints kit tag"      || fail "no kit tag in output"
echo "$OUT" | grep -q "GNU Fortran (Test)" && ok "prints kit compiler" || fail "no kit compiler in output"
echo "$OUT" | grep -q "OK: module version $LOCAL_ABI" \
    && ok "reports the matching ABI" || fail "no OK line"

# --- mismatched kit ------------------------------------------------------
OUT="$("$SCRIPT" "$TMPD/bad" gfortran linux 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "mismatched kit exits 1" || fail "mismatched kit exited $RC"
echo "$OUT" | grep -q "version 99"        && ok "names the kit ABI"   || fail "kit ABI missing"
echo "$OUT" | grep -q "version $LOCAL_ABI" && ok "names the local ABI" || fail "local ABI missing"

# platform-specific remedies
"$SCRIPT" "$TMPD/bad" gfortran macos   2>&1 | grep -q "brew install gcc" \
    && ok "macos remedy"   || fail "macos remedy missing"
"$SCRIPT" "$TMPD/bad" gfortran windows 2>&1 | grep -q "pacman" \
    && ok "windows remedy" || fail "windows remedy missing"
"$SCRIPT" "$TMPD/bad" gfortran linux   2>&1 | grep -q "apt install" \
    && ok "linux remedy"   || fail "linux remedy missing"

# --- kit with no BUILDINFO ----------------------------------------------
OUT="$("$SCRIPT" "$TMPD/nobi" gfortran linux 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 0 ] && ok "kit without BUILDINFO still exits 0" || fail "exited $RC"
echo "$OUT" | grep -q "no BUILDINFO.txt" && ok "says BUILDINFO is absent" || fail "no notice"

# --- unprobeable compiler ------------------------------------------------
OUT="$("$SCRIPT" "$TMPD/good" definitely-not-a-compiler linux 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 0 ] && ok "unprobeable compiler exits 0 (lenient)" || fail "exited $RC"
echo "$OUT" | grep -q "Skipped" && ok "reports skip" || fail "no skip notice"

# --- missing kit dir -----------------------------------------------------
"$SCRIPT" "$TMPD/does-not-exist" gfortran linux >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing kit dir exits 1" || fail "missing kit dir should exit 1"

# --- mktemp -d failure must fail clearly, not fall through to "Skipped" --
# A failed mktemp used to leave TMPD empty, produce a raw shell error from
# the probe compile, and then report "Skipped" and exit 0 -- silently
# passing the ABI tripwire. We can't make the real mktemp fail cleanly on
# demand, so shadow it on PATH with one that always fails, the same
# technique test_refresh_kit.sh uses.
FAKEBIN="$TMPD/fakebin"
mkdir -p "$FAKEBIN"
cat > "$FAKEBIN/mktemp" <<'EOF'
#!/usr/bin/env bash
echo "mktemp: simulated failure" >&2
exit 1
EOF
chmod +x "$FAKEBIN/mktemp"

OUT="$(PATH="$FAKEBIN:$PATH" "$SCRIPT" "$TMPD/good" gfortran linux 2>&1)"; RC=$?
echo "$OUT"
[ $RC -eq 1 ] && ok "mktemp failure exits 1" || fail "mktemp failure exited $RC"
echo "$OUT" | grep -qi "Skipped" && fail "mktemp failure should not fall through to the lenient Skipped path" \
    || ok "mktemp failure did not silently skip"

echo ""
if [ "$FAILURES" -eq 0 ]; then
    echo "PASS: all assertions held"
    exit 0
fi
echo "$FAILURES failure(s)"
exit 1
