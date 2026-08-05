# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

URAMSES is a Fortran framework for integrating custom user-defined models into PyRAMSES (Python) and STEPSS (Java) power system simulators. Users write Fortran model subroutines, register them in router files, and compile everything into a shared library (`ramses.so` / `ramses.dll`) or standalone executable (`dynsim`).

## Build Commands (Linux)

```bash
make -f build/Makefile.linux            # Build shared library + executable
make -f build/Makefile.linux dll        # Shared library only (ramses.so)
make -f build/Makefile.linux exe        # Executable only (dynsim)
make -f build/Makefile.linux clean      # Remove build artifacts
make -f build/Makefile.linux check-deps # Verify dependencies (gfortran, OpenBLAS, libramses.a)
```

Output goes to `Release_l/`.

## Build Commands (macOS and Windows/MinGW)

`build/Makefile.macos` and `build/Makefile.windows` mirror
`build/Makefile.linux`, with the same targets
(`all`/`dll`/`exe`/`clean`/`check-deps`/`help`) and the same wildcard model
discovery. They differ only in the module directory, output directory, and
library naming.

Every platform carries one suffix throughout: `l` (Linux), `m` (macOS), `wg`
(Windows/MinGW gfortran) and `wi` (Windows/Intel). A kit directory and the
output directory it builds into always share it.

| Route | Makefile | Modules | Output dir | Shared lib |
|---|---|---|---|---|
| Linux | `build/Makefile.linux` | `modules_l/` (`libramses.a`) | `Release_l/` | `ramses.so` |
| macOS arm64 | `build/Makefile.macos` | `modules_m/` (`libramses.a`) | `Release_m/` | `ramses.so` |
| Windows MinGW | `build/Makefile.windows` | `modules_wg/` (`libramses.lib`) | `Release_wg/` | `ramses.dll` |
| Windows Intel | `build/msvs/URAMSES.sln` | `modules_wi/` (`libramses.lib`) | `Release_wi/` | `ramses.dll` |

Constraints worth knowing before editing these files:

- The Makefiles live in `build/` but are **always invoked from the repo root**
  (`make -f build/Makefile.linux`). `-f` does not change make's working
  directory, so every path inside them stays repo-root-relative. Do not add
  `../` prefixes, and do not `cd build` to run them.
- The `.vfproj` files live in `build/msvs/` and their paths *are* relative to
  their own directory, two levels below the root, hence the `..\..\src\`,
  `..\..\custom_models\`, `..\..\modules_wi` and `..\..\Release_wi` prefixes.
  The two routes resolve paths differently; keep them straight.
- macOS uses `ramses.so`, **not** `ramses.dylib`. PyRAMSES separates platforms
  by the `libs/lin`, `libs/mac`, `libs/win` folder, not by filename.
- `modules_m/` is arm64; macOS builds require Apple Silicon.
- Each kit's GFORTRAN module ABI version is recorded in its
  `modules_<plat>/BUILDINFO.txt` (see the README's "Which compiler do I need?"
  section) rather than hardcoded here, because it changes whenever RAMSES bumps
  a runner's compiler. Each `check-deps` compares the local compiler against
  the kit and fails early on a mismatch.
- The Windows gfortran kit ships `libramses.lib`, which MinGW's linker does not
  probe for `-lramses`, so `build/Makefile.windows` passes it by explicit path.
- `FC` is assigned with `=` not `?=` in all three: make predefines `FC` as `f77`,
  which `?=` would leave in place. A command-line `FC=` still overrides.

Kits come from the matching `stepss-ramses` release, e.g.
`uramses-modules_{l,m,wg}-v3.51.zip` for RAMSES v3.51.

## Module kits are CI-managed

`modules_l/`, `modules_m/` and `modules_wg/` are written by
`.github/workflows/sync-ramses-release.yml` and must not be hand-edited. When
stepss-ramses publishes a release it dispatches here; the workflow refreshes
all three kits from that release's `uramses-modules_*` bundles, builds every
Makefile on the same runner images RAMSES built on, and only then
fast-forwards `master` and publishes a matching release under the same tag.

Every kit carries a `BUILDINFO.txt` recording the compiler, `.mod` ABI
version, flags, BLAS and runtime floor behind it. `tools/check_kit.sh` reads it
during `check-deps` and refuses a compiler that cannot read the kit.

`modules_wi/` (Intel, consumed by `build/msvs/URAMSES.sln`) is outside all of this and is
still refreshed by hand.

The runner images in the sync workflow deliberately mirror those in
stepss-ramses' `release.yml` — currently `ubuntu-24.04` (gfortran 13.3, `.mod`
ABI 15). Do not replace that pin with `ubuntu-latest`: the kit is built by the
pinned image over in stepss-ramses, so a floating image here would pair the
next compiler GitHub promotes against a kit built by a different one. Bumping
the image is fine, but bump both repositories in the same release.

The `.mod` ABI is coarser than the compiler version — gfortran 9 through 13 all
write version 15, and 14 is the next bump — so an image change is not
automatically an ABI break. `ubuntu-22.04` → `ubuntu-24.04` (gfortran 11.4 →
13.3) kept ABI 15. What it did move is the **glibc floor**, 2.35 → 2.39: the
published Linux binaries no longer run on Ubuntu 22.04 or Debian 12.

The three kits are not on a common ABI and never have been. As of v3.55:
`modules_l` is gfortran 13.3 / ABI 15, while `modules_m` and `modules_wg` are
gfortran 16.1 / ABI 16, because `macos-15` and `windows-latest` track much
newer toolchains than any Ubuntu LTS image. Read the per-kit `BUILDINFO.txt`
rather than assuming one number covers all three; `check_kit.sh` compares
against the kit in hand, so it is already right about this.

Helper scripts in `tools/` each have a `tools/test_*.sh` alongside them. Run
them after any change:

```sh
bash tools/test_check_kit.sh
bash tools/test_refresh_kit.sh
bash tools/test_render_release_notes.sh
```

## The model routers come from RAMSES

`src/usr_*_models.f90` are the five model routers, and they are **copies of
`stepss-ramses/src/devices/usr_*_models.f90`**, not independent files. The
linker prefers an explicitly-listed object over an archive member, so these
shadow the copies inside `libramses.a`: whatever they do not register becomes
unreachable, however many models the kit actually contains.

They were once stripped to empty templates (every `case` commented out), which
left a URAMSES build able to resolve exactly one model, `VFAULT`, out of the
~57 compiled into the kit. The dispatcher tries the router first and falls back
to only `CONSTANT`, `1ST_ORDER`, `GENERIC1` and `GENERIC2`, so every named
model (ST1A, AC1A, SEXS, DEGOV1, GAST, HYGOV, HVDC_*, PQ, IBG, WT3/WT4, BESS,
GFOL/GFOR, PMU, vfd_load) was unreachable. Do not strip them again.

Notes when editing:

- Match the **prefixed** name: the routers normalise a model name by
  prepending `exc_`/`inj_`/`tor_`/`twop_`/`dctl_` when absent, then select on
  that, so the label is `case('inj_vfd_load')`.
- Register models from `custom_models/` alongside the pre-compiled ones, as
  `exc_ENTSOE_lim` is in `usr_exc_models.f90`. Do not register a model the kit
  already exports, or the link fails on a duplicate symbol.
- These routers carry a deliberately narrower model list than the RAMSES
  originals. Re-copying those files wholesale reintroduces entries that were
  removed on purpose; diff before replacing.
- `FUNCTIONS_IN_MODELS` is **not** built here. It comes from the kit, as
  `modules_<plat>/functions_in_models.mod` to compile against and inside
  `libramses.a` to link against. A local copy collides with the archive's the
  moment anything pulls that member in.

## Regression gate

`tools/nordic_gate.sh` runs the Nordic voltage-collapse case through `dynsim`
and compares the trajectory against `tests/baselines/nordic_baseline.npz`, the
baseline generated in stepss-ramses. The case does not touch `custom_models/`,
so a correct build matches it exactly. All three gfortran jobs in the sync
workflow run it after building; it needs numpy on the runner.

```sh
make -f build/Makefile.linux all
bash tools/nordic_gate.sh Release_l/dynsim . tests/baselines/nordic_baseline.npz
```

## Architecture

**Build dependency chain:**
```
main.f90 → c_interface.f90 → usr_*_models.f90 → custom_models/*.f90 + FUNCTIONS_IN_MODELS.f90
                                                        ↓
                                                 libramses.a (pre-compiled)
```

**Key directories:**
- `src/`: framework code (C interface, model routers, main entry point)
- `custom_models/`: user model implementations (auto-discovered by every Makefile build)
- `examples/Nordic/`, `tests/baselines/`: the regression case and baseline behind `tools/nordic_gate.sh`
- `build/`: the three Makefiles, plus `build/msvs/` for the Visual Studio route
- `modules_l/`: pre-compiled RAMSES library and `.mod` files (Linux/gfortran)
- `modules_m/`: pre-compiled RAMSES library and `.mod` files (macOS arm64/gfortran)
- `modules_wg/`: pre-compiled RAMSES library and `.mod` files (Windows/MinGW)
- `modules_wi/`: pre-compiled RAMSES library and `.mod` files (Windows/Intel)

**Model router pattern** (`src/usr_*_models.f90`): Each file contains an `assoc_*_ptr` subroutine that maps string model names to Fortran subroutine pointers via `select case`. Five model categories exist: `exc` (exciters), `inj` (injectors), `tor` (torque/governors), `twop` (two-port devices), `dctl` (discrete control).

**Model subroutine pattern**: Each model in `custom_models/` is a single subroutine using mode-based dispatch (`select case` on `mode`): `define_var_and_par`, `define_obs`, `diffstate`, `algstate`, `update_state`. Models use `FUNCTIONS_IN_MODELS` for shared utilities (`ppower`, `qpower`, `vrectif`, `vcomp`, etc.).

## Adding a New Model

1. Create `custom_models/<type>_<NAME>.f90` following the naming convention (`exc_`, `inj_`, `tor_`, `twop_`, `dctl_`)
2. Register in the corresponding `src/usr_<type>_models.f90` by adding a `case` to the `select case` block
3. Rebuild with the Makefile for your platform, e.g. `make -f build/Makefile.linux clean all`

All three Makefiles auto-discover `.f90` files in `custom_models/` via wildcard. Only the Intel Visual Studio route requires files to be added to the project by hand.

Do not copy a model that the linked RAMSES library already exports: defining a subroutine the library also defines produces a duplicate-symbol link error. Register those by name in the routers instead, as `VFAULT` is in `src/usr_inj_models.f90`.

## Fortran Conventions

- Free-form Fortran 90/95, compiled with `-ffree-line-length-none`
- Double precision throughout
- `#ifdef DLL` preprocessor guards for Windows DLL exports (`!DEC$ ATTRIBUTES DLLEXPORT`)
- OpenMP enabled (`-fopenmp`)
- `FUNCTIONS_IN_MODELS.f90` must compile before any model that `use`s it (enforced in Makefile)

## No Test Suite

There is no automated test framework. Models are validated through integration with PyRAMSES or STEPSS simulations.
