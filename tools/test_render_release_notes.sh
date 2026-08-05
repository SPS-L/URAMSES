#!/usr/bin/env bash
# Local test for tools/render_release_notes.sh.
#
#   test_render_release_notes.sh
#
# Renders notes from three synthetic BUILDINFO files plus a RAMSES body and
# asserts every platform is represented, the RAMSES prose survives verbatim,
# and the Intel-kit caveat is present.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/tools/render_release_notes.sh"
FAILURES=0

fail() { echo "FAIL: $*"; FAILURES=$((FAILURES + 1)); }
ok()   { echo "ok: $*"; }

TMPD="$(mktemp -d)"
trap 'rm -rf "$TMPD"' EXIT

mk() {  # mk <file> <kit_dir> <compiler> <abi> <target> <floor> <consume_with>
    cat > "$1" <<EOF
kit_dir          $2
ramses_tag       v9.99
ramses_commit    deadbeef
built_utc        2026-09-01T10:12:44Z
runner_image     test-runner
compiler         $3
target           $5
mod_abi_version  $4
fflags           -O2 -cpp -fopenmp
ldflags          -fopenmp
blas             OpenBLAS 0.3.20
runtime_floor    $6
consume_with     $7
EOF
}

mk "$TMPD/lin.txt" modules_lin "GNU Fortran (Ubuntu) 11.4.0" 15 \
   x86_64-pc-linux-gnu "glibc 2.35" Makefile.linux
mk "$TMPD/mac.txt" modules_mac "GNU Fortran (Homebrew) 15.1.0" 16 \
   aarch64-apple-darwin24 "macOS 15 arm64" Makefile.macos
mk "$TMPD/win.txt" modules_win_gfortran "GNU Fortran (MinGW) 14.2.0" 16 \
   x86_64-w64-mingw32 "static" Makefile.windows

cat > "$TMPD/body.md" <<'EOF'
## What's new

* Faster KLU refactorisation on the decomposed solver.
* Fixed a `CHGPRM` edge case | with a pipe character in it.
EOF

# --- usage ---------------------------------------------------------------
"$SCRIPT" >/dev/null 2>&1
[ $? -eq 2 ] && ok "no args exits 2" || fail "no args should exit 2"

"$SCRIPT" "$TMPD/nope.md" "$TMPD/lin.txt" "$TMPD/mac.txt" "$TMPD/win.txt" >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing body exits 1" || fail "missing body should exit 1"

"$SCRIPT" "$TMPD/body.md" "$TMPD/nope.txt" "$TMPD/mac.txt" "$TMPD/win.txt" >/dev/null 2>&1
[ $? -eq 1 ] && ok "missing buildinfo exits 1" || fail "missing buildinfo should exit 1"

# --- happy path ----------------------------------------------------------
OUT="$("$SCRIPT" "$TMPD/body.md" "$TMPD/lin.txt" "$TMPD/mac.txt" "$TMPD/win.txt")"; RC=$?
echo "$OUT"
echo "--- end of rendered notes ---"
[ $RC -eq 0 ] && ok "renders, exit 0" || fail "exited $RC"

echo "$OUT" | grep -q "## Toolchain and kit provenance" && ok "has the provenance heading" \
    || fail "provenance heading missing"

for m in Makefile.linux Makefile.macos Makefile.windows; do
    echo "$OUT" | grep -q "$m" && ok "row for $m" || fail "no row for $m"
done
for d in modules_lin modules_mac modules_win_gfortran; do
    echo "$OUT" | grep -q "$d" && ok "names $d" || fail "does not name $d"
done
echo "$OUT" | grep -q "11.4.0" && ok "linux compiler shown" || fail "linux compiler missing"
echo "$OUT" | grep -q "15.1.0" && ok "macos compiler shown" || fail "macos compiler missing"
echo "$OUT" | grep -q "14.2.0" && ok "windows compiler shown" || fail "windows compiler missing"

# three data rows plus a header and a separator
ROWS="$(echo "$OUT" | grep -c '^| ')"
[ "$ROWS" -eq 5 ] && ok "table has 5 pipe rows" || fail "table has $ROWS pipe rows, want 5"

echo "$OUT" | grep -q "BUILDINFO.txt" && ok "points at BUILDINFO.txt" || fail "no BUILDINFO pointer"
echo "$OUT" | grep -qi "Intel" && ok "carries the Intel-kit caveat" || fail "Intel caveat missing"

# the RAMSES body must survive verbatim, pipe character and all
echo "$OUT" | grep -q "Faster KLU refactorisation on the decomposed solver." \
    && ok "RAMSES prose preserved" || fail "RAMSES prose lost"
echo "$OUT" | grep -q 'edge case | with a pipe character in it' \
    && ok "RAMSES body copied verbatim, not escaped" || fail "RAMSES body was mangled"
echo "$OUT" | grep -q '^---$' && ok "separator before the RAMSES body" || fail "no separator"

echo ""
if [ "$FAILURES" -eq 0 ]; then
    echo "PASS: all assertions held"
    exit 0
fi
echo "$FAILURES failure(s)"
exit 1
