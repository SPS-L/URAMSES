# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

URAMSES is a Fortran framework for integrating custom user-defined models into PyRAMSES (Python) and STEPSS (Java) power system simulators. Users write Fortran model subroutines, register them in router files, and compile everything into a shared library (`ramses.so` / `ramses.dll`) or standalone executable (`dynsim`).

## Build Commands (Linux)

```bash
make -f Makefile.linux            # Build shared library + executable
make -f Makefile.linux dll        # Shared library only (ramses.so)
make -f Makefile.linux exe        # Executable only (dynsim)
make -f Makefile.linux clean      # Remove build artifacts
make -f Makefile.linux check-deps # Verify dependencies (gfortran, OpenBLAS, libramses.a)
```

Output goes to `Release_gnu_l/`.

## Build Commands (macOS and Windows/MinGW)

`Makefile.macos` and `Makefile.mingw` mirror `Makefile.linux` — same targets
(`all`/`dll`/`exe`/`clean`/`check-deps`/`help`), same wildcard model discovery —
and differ only in the module directory, output directory, and library naming:

| Route | Makefile | Modules | Output dir | Shared lib |
|---|---|---|---|---|
| Linux | `Makefile.linux` | `modules_lin/` (`libramses.a`) | `Release_gnu_l/` | `ramses.so` |
| macOS arm64 | `Makefile.macos` | `modules_mac/` (`libramses.a`) | `Release_gnu_m/` | `ramses.so` |
| Windows MinGW | `Makefile.mingw` | `modules_win_gfortran/` (`libramses.lib`) | `Release_gnu_w/` | `ramses.dll` |
| Windows Intel | `URAMSES.sln` | `modules/` (`libramses.lib`) | `Release_intel_w64/` | `ramses.dll` |

Constraints worth knowing before editing these files:

- macOS uses `ramses.so`, **not** `ramses.dylib` — PyRAMSES separates platforms
  by the `libs/lin`, `libs/mac`, `libs/win` folder, not by filename.
- `modules_mac/` is arm64; macOS builds require Apple Silicon.
- `modules_mac/` and `modules_win_gfortran/` are GFORTRAN module version 16
  (gfortran 15+); `modules_lin/` is version 15. Each `check-deps` compares the
  local compiler against the kit and fails early on a mismatch.
- The Windows gfortran kit ships `libramses.lib`, which MinGW's linker does not
  probe for `-lramses`, so `Makefile.mingw` passes it by explicit path.
- `FC` is assigned with `=` not `?=` in all three: make predefines `FC` as `f77`,
  which `?=` would leave in place. A command-line `FC=` still overrides.

Kits come from the matching `stepss-ramses` release, e.g.
`uramses-modules_{lin,mac,win_gfortran}-v3.51.zip` for RAMSES v3.51.

## Architecture

**Build dependency chain:**
```
main.f90 → c_interface.f90 → usr_*_models.f90 → my_models/*.f90 + FUNCTIONS_IN_MODELS.f90
                                                        ↓
                                                 libramses.a (pre-compiled)
```

**Key directories:**
- `src/` — Framework code (C interface, model routers, utility functions, main entry point)
- `my_models/` — User model implementations (auto-discovered by every Makefile build)
- `modules_lin/` — Pre-compiled RAMSES library and `.mod` files (Linux/gfortran)
- `modules_mac/` — Pre-compiled RAMSES library and `.mod` files (macOS arm64/gfortran)
- `modules_win_gfortran/` — Pre-compiled RAMSES library and `.mod` files (Windows/MinGW)
- `modules/` — Pre-compiled RAMSES library and `.mod` files (Windows/Intel)

**Model router pattern** (`src/usr_*_models.f90`): Each file contains an `assoc_*_ptr` subroutine that maps string model names to Fortran subroutine pointers via `select case`. Five model categories exist: `exc` (exciters), `inj` (injectors), `tor` (torque/governors), `twop` (two-port devices), `dctl` (discrete control).

**Model subroutine pattern**: Each model in `my_models/` is a single subroutine using mode-based dispatch (`select case` on `mode`): `define_var_and_par`, `define_obs`, `diffstate`, `algstate`, `update_state`. Models use `FUNCTIONS_IN_MODELS` for shared utilities (`ppower`, `qpower`, `vrectif`, `vcomp`, etc.).

## Adding a New Model

1. Create `my_models/<type>_<NAME>.f90` following the naming convention (`exc_`, `inj_`, `tor_`, `twop_`, `dctl_`)
2. Register in the corresponding `src/usr_<type>_models.f90` by adding a `case` to the `select case` block
3. Rebuild with the Makefile for your platform, e.g. `make -f Makefile.linux clean all`

All three Makefiles auto-discover `.f90` files in `my_models/` via wildcard. Only the Intel Visual Studio route requires files to be added to the project by hand.

Do not copy a model that the linked RAMSES library already exports — defining a subroutine the library also defines produces a duplicate-symbol link error. Register those by name in the routers instead, as `VFAULT` is in `src/usr_inj_models.f90`.

## Fortran Conventions

- Free-form Fortran 90/95, compiled with `-ffree-line-length-none`
- Double precision throughout
- `#ifdef DLL` preprocessor guards for Windows DLL exports (`!DEC$ ATTRIBUTES DLLEXPORT`)
- OpenMP enabled (`-fopenmp`)
- `FUNCTIONS_IN_MODELS.f90` must compile before any model that `use`s it (enforced in Makefile)

## No Test Suite

There is no automated test framework. Models are validated through integration with PyRAMSES or STEPSS simulations.
