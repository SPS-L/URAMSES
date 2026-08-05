#!/usr/bin/env bash
# Render the body of a URAMSES release.
#
#   render_release_notes.sh <body-file> <buildinfo-lin> <buildinfo-mac> <buildinfo-win>
#
# Writes to stdout a generated toolchain table describing the three module
# kits, followed verbatim by <body-file> (the RAMSES release body). URAMSES
# ships pre-compiled kits, so the compiler behind each one is the single fact
# a user needs before adding their own Fortran models.
set -u

if [ $# -ne 4 ]; then
    echo "usage: $0 <body-file> <buildinfo-lin> <buildinfo-mac> <buildinfo-win>" >&2
    exit 2
fi
BODY="$1"
BI_LIN="$2"
BI_MAC="$3"
BI_WIN="$4"

for f in "$BODY" "$BI_LIN" "$BI_MAC" "$BI_WIN"; do
    [ -f "$f" ] || { echo "FAIL: input not found: $f" >&2; exit 1; }
done

# One aligned key/value pair out of a BUILDINFO.txt.
bi() { sed -n "s/^$2[[:space:]]\{1,\}//p" "$1" | head -n1; }

# `consume_with` is stamped into the kit by stepss-ramses CI, which now writes
# the qualified `build/Makefile.<plat>` and needs no help here. Kits stamped
# before that change carry a bare `Makefile.<plat>`, so prefix those rather
# than publish a path that no longer resolves. Qualified values pass through
# untouched; both directions are pinned in test_render_release_notes.sh.
row() {   # row <buildinfo-file>
    consume="$(bi "$1" consume_with)"
    case "$consume" in
        Makefile.*) consume="build/$consume" ;;
    esac
    printf '| `%s` | `%s` | %s | %s | `%s` | %s |\n' \
        "$consume" \
        "$(bi "$1" kit_dir)" \
        "$(bi "$1" compiler)" \
        "$(bi "$1" mod_abi_version)" \
        "$(bi "$1" target)" \
        "$(bi "$1" runtime_floor)"
}

cat <<'EOF'
## Toolchain and kit provenance

The `modules_*/` directories in this release ship pre-compiled RAMSES
libraries and `.mod` files. A `.mod` file can only be read by the gfortran
generation that wrote it, so to add your own Fortran models you need a
compiler emitting the same module ABI version as your platform's kit.
`make -f <makefile> check-deps` verifies this before building.

EOF

echo '| Makefile | Kit directory | Built with | `.mod` ABI | Target | Runtime floor |'
echo '| --- | --- | --- | --- | --- | --- |'
row "$BI_LIN"
row "$BI_MAC"
row "$BI_WIN"

cat <<EOF

Exact compile flags, BLAS build and RAMSES commit for each kit are recorded in
that kit's \`BUILDINFO.txt\`, for example \`modules_l/BUILDINFO.txt\`.

The Windows/Intel kit in \`modules_wi/\` (used by \`build/msvs/URAMSES.sln\`) is
**not** refreshed by this release: no CI produces Intel-format modules. See
\`modules_wi/BUILDINFO.txt\` for what it was last built from.

---

EOF

cat "$BODY"
