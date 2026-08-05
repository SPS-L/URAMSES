#!/usr/bin/env bash
# Report a module kit's provenance and verify the local gfortran can read it.
#
#   check_kit.sh <mod-dir> <fortran-compiler> <linux|macos|windows>
#
# Prints the RAMSES tag and build compiler recorded in the kit's
# BUILDINFO.txt, then compares the GFORTRAN module ABI version baked into the
# kit's .mod files against the version the local compiler emits. A .mod file
# can only be read by the gfortran generation that wrote it, so a mismatch is
# a hard error carrying a platform-specific remedy.
#
# Exits 0 on match or when either version cannot be determined (the same
# lenient behaviour the Makefiles had before this script existed), 1 on
# mismatch or a missing kit, 2 on usage error.
set -u

if [ $# -ne 3 ]; then
    echo "usage: $0 <mod-dir> <fortran-compiler> <linux|macos|windows>" >&2
    exit 2
fi
MOD_DIR="$1"
FC="$2"
PLATFORM="$3"

case "$PLATFORM" in
    linux|macos|windows) ;;
    *) echo "usage: $0 <mod-dir> <fortran-compiler> <linux|macos|windows>" >&2; exit 2 ;;
esac

[ -d "$MOD_DIR" ] || { echo "ERROR: module directory '$MOD_DIR' not found."; exit 1; }

# Read one aligned key/value pair out of a BUILDINFO.txt.
bi() { sed -n "s/^$1[[:space:]]\{1,\}//p" "$MOD_DIR/BUILDINFO.txt" | head -n1; }

# Read the GFORTRAN module ABI integer out of a .mod file (they are gzipped
# text with a version banner). Empty output means unreadable.
mod_abi() {
    gzip -dc "$1" 2>/dev/null | head -c 256 \
        | sed -n "s/.*module version '\([0-9]*\)'.*/\1/p"
}

if [ -f "$MOD_DIR/BUILDINFO.txt" ]; then
    echo "  Kit:      RAMSES $(bi ramses_tag), built $(bi built_utc) on $(bi runner_image)"
    echo "  Built by: $(bi compiler)"
else
    echo "  Kit:      no BUILDINFO.txt (kit predates provenance tracking)"
fi

# The kit's ABI is read from its own .mod files rather than from BUILDINFO,
# so a hand-patched kit cannot lie its way past this check.
FIRST_MOD="$(ls -1 "$MOD_DIR"/*.mod 2>/dev/null | sort | head -n1)"
KIT_VER=""
[ -n "$FIRST_MOD" ] && KIT_VER="$(mod_abi "$FIRST_MOD")"

TMPD="$(mktemp -d)" || { echo "ERROR: could not create a temp directory (mktemp -d failed)."; exit 1; }
[ -n "$TMPD" ] && [ -d "$TMPD" ] || {
    echo "ERROR: could not create a temp directory (mktemp -d returned no usable path)."
    exit 1
}
trap 'rm -rf "$TMPD"' EXIT
printf 'module probe_mod\nend module probe_mod\n' > "$TMPD/probe.f90"
"$FC" -c "$TMPD/probe.f90" -o "$TMPD/probe.o" -J"$TMPD" >/dev/null 2>&1 || true
FC_VER="$(mod_abi "$TMPD/probe_mod.mod")"

if [ -z "$KIT_VER" ] || [ -z "$FC_VER" ]; then
    echo "  Skipped: could not determine module versions"
    exit 0
fi

if [ "$KIT_VER" = "$FC_VER" ]; then
    echo "  OK: module version $KIT_VER matches $FC"
    exit 0
fi

echo "ERROR: $MOD_DIR/ holds GFORTRAN module version $KIT_VER but $FC emits version $FC_VER."
echo "The .mod files can only be read by the gfortran release that wrote them."
case "$PLATFORM" in
    linux)
        echo "Install a matching gfortran and pass it explicitly, e.g. on Ubuntu"
        echo "(the kit's own compiler is printed above; adjust the number to match):"
        echo "  sudo apt install gfortran-13"
        echo "  make -f build/Makefile.linux FC=gfortran-13"
        ;;
    macos)
        echo "Install a matching gfortran (brew install gcc) and if it is versioned, pass it:"
        echo "  make -f build/Makefile.macos FC=gfortran-15"
        ;;
    windows)
        echo "Update the MSYS2 toolchain from an MSYS2 shell:"
        echo "  pacman -Syu mingw-w64-x86_64-gcc-fortran"
        ;;
esac
exit 1
