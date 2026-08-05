# URAMSES

**User-defined device model framework for the RAMSES power system simulator**

URAMSES lets you compile your own Fortran device models and link them against a pre-compiled RAMSES library, as part of the [STEPSS](https://stepss.sps-lab.org/) power system simulation platform. You write a model subroutine, register it in a router file, and build either a shared library (`ramses.so`/`ramses.dll`) loaded by [PyRAMSES](https://stepss.sps-lab.org/pyramses/) or a standalone executable (`dynsim`).

## Features

- **Five model categories** — exciters (`exc_`), torque/governors (`tor_`), injectors (`inj_`), two-port devices (`twop_`), and discrete controls (`dctl_`)
- **No RAMSES source required** — models link against a pre-compiled RAMSES library shipped per platform: `modules_lin/` (Linux/gfortran, `libramses.a`), `modules_mac/` (macOS arm64/gfortran, `libramses.a`), `modules_win_gfortran/` (Windows/MinGW, `libramses.lib`), and `modules/` (Windows/Intel, `libramses.lib`)
- **Automatic model discovery in every Makefile build** — Linux, macOS, and Windows/MinGW pick up any `.f90` file placed in `my_models/` via wildcard; no Makefile edits needed. Only the Intel Visual Studio route needs files added by hand
- **Four build routes** — gfortran Makefiles for Linux, macOS, and Windows/MinGW; or Visual Studio with Intel oneAPI Fortran (Windows)
- **Dual outputs** — a shared library for PyRAMSES/STEPSS integration and a standalone `dynsim` executable
- **Example models included** — an ENTSO-E exciter and an air-conditioning load model with parameter files in `my_models/`

## Installation

**Requirements:** gfortran and OpenBLAS (Linux, macOS, Windows/MinGW); or Visual Studio 2019+ and the Intel oneAPI Fortran compiler (Windows/Intel). To use the compiled models from Python you also need [PyRAMSES](https://stepss.sps-lab.org/pyramses/) (`pip install pyramses`).

> **gfortran version must match the kit.** Fortran `.mod` files are readable only
> by the gfortran release that wrote them. The `modules_mac/` and
> `modules_win_gfortran/` kits are built with **GFORTRAN module version 16**
> (gfortran 15 or newer); `modules_lin/` uses **module version 15** (gfortran 13
> or 14). Each Makefile's `check-deps` target compares your compiler against the
> kit and reports a mismatch before the build starts.

### Linux

```bash
# Ubuntu/Debian
sudo apt install gfortran libopenblas-dev

# Fedora/RHEL
sudo dnf install gcc-gfortran openblas-devel

# Arch Linux
sudo pacman -S gcc-fortran openblas
```

### macOS (Apple Silicon)

macOS ships no Fortran compiler, so gfortran comes from Homebrew:

```bash
brew install gcc openblas
```

macOS builds require **Apple Silicon** — the `modules_mac/` kit is arm64. If
Homebrew installs a versioned binary only, pass it explicitly with
`FC=gfortran-15`.

### Windows (MinGW / gfortran)

Install [MSYS2](https://www.msys2.org/), then from an MSYS2 shell:

```bash
pacman -S mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-openblas
```

Build from the **"MSYS2 MinGW 64-bit"** shell, not the plain MSYS shell — the
latter links against the `msys-2.0` runtime and produces a DLL that CPython
cannot load. The Makefile refuses to run there.

### Windows (Intel oneAPI)

- **Microsoft Visual Studio** (2019 or later recommended)
- **Intel oneAPI Fortran Compiler** (formerly Intel Fortran)
- **PyRAMSES** (Python package) or the **STEPSS** Java interface

For detailed installation instructions of the Intel oneAPI Fortran compiler, refer to the included PDF:
[Installing the Intel oneAPI Fortran compiler.pdf](Installing%20the%20Intel%20oneAPI%20Fortran%20compiler.pdf)

## Quick Start

### Linux

```bash
# Build everything (shared library + executable)
make -f Makefile.linux

# Run standalone simulation
./Release_gnu_l/dynsim
```

### macOS (Apple Silicon)

```bash
# Build everything (shared library + executable)
make -f Makefile.macos

# Run standalone simulation
./Release_gnu_m/dynsim
```

To link against Apple's Accelerate framework instead of OpenBLAS:

```bash
make -f Makefile.macos BLAS=accelerate
```

### Windows (MinGW / gfortran)

From the "MSYS2 MinGW 64-bit" shell:

```bash
# Build everything (shared library + executable)
make -f Makefile.windows

# Run standalone simulation
./Release_gnu_w/dynsim.exe
```

### Windows (Intel oneAPI)

1. Open `URAMSES.sln` in Visual Studio
2. Select `Release|x64` configuration
3. Build → Build Solution
4. Run `Release_intel_w64\dynsim.exe`

### Load your models in PyRAMSES

```python
import pyramses

ram = pyramses.sim('/path/to/URAMSES/Release_gnu_l')             # Linux
# ram = pyramses.sim('/path/to/URAMSES/Release_gnu_m')           # macOS
# ram = pyramses.sim(r'C:\path\to\URAMSES\Release_gnu_w')        # Windows (MinGW)
# ram = pyramses.sim(r'C:\path\to\URAMSES\Release_intel_w64')    # Windows (Intel)
```

PyRAMSES loads `ramses.so` on Linux and macOS, and `ramses.dll` on Windows.

## Model Types

URAMSES supports several types of power system models:

- **Exciters (`exc_*`)**: Generator excitation system models
- **Injectors (`inj_*`)**: Current/voltage injection models for faults/disturbances
- **Torque (`tor_*`)**: Mechanical torque models for generators
- **Two-port (`twop_*`)**: Two-port network models (e.g., SVC, STATCOM)
- **Discrete Control (`dctl_*`)**: Discrete control system models

## Project Structure

```
URAMSES/
├── src/                    # Source code files (common for Windows/Linux)
│   ├── c_interface.f90     # C interface for Python integration
│   ├── main.f90            # Main entry point (executable only)
│   ├── FUNCTIONS_IN_MODELS.f90  # Helper functions for models
│   ├── usr_exc_models.f90  # Exciter model associations
│   ├── usr_inj_models.f90  # Injector model associations
│   ├── usr_tor_models.f90  # Torque model associations
│   ├── usr_twop_models.f90 # Two-port model associations
│   └── usr_dctl_models.f90 # Discrete control model associations
├── my_models/              # Your custom models (common for Windows/Linux)
│   ├── exc_*.f90           # Exciter models
│   ├── inj_*.f90           # Injector models
│   ├── tor_*.f90           # Torque models
│   └── *.txt               # Model parameter files
├── modules/                # Pre-compiled modules (Windows/Intel Fortran)
│   ├── *.mod               # Module interface files
│   └── libramses.lib       # Pre-compiled RAMSES library
├── modules_lin/            # Pre-compiled modules (Linux/gfortran)
│   ├── *.mod               # Module interface files
│   └── libramses.a         # Pre-compiled RAMSES library
├── modules_mac/            # Pre-compiled modules (macOS arm64/gfortran)
│   ├── *.mod               # Module interface files
│   └── libramses.a         # Pre-compiled RAMSES library
├── modules_win_gfortran/   # Pre-compiled modules (Windows/MinGW gfortran)
│   ├── *.mod               # Module interface files
│   └── libramses.lib       # Pre-compiled RAMSES library
├── URAMSES.sln             # Visual Studio solution file (Windows/Intel)
├── dllramses.vfproj        # DLL project - ramses.dll (Windows/Intel)
├── exeramses.vfproj        # Executable project - dynsim.exe (Windows/Intel)
├── MDL.vfproj              # Model library project - ramsesmdl.dll (Windows/Intel)
├── Makefile.linux          # Makefile for Linux builds
├── Makefile.macos          # Makefile for macOS builds (Apple Silicon)
├── Makefile.windows        # Makefile for Windows/MinGW builds
├── Release_intel_w64/      # Compiled output (Windows/Intel)
├── Release_gnu_w/          # Compiled output (Windows/MinGW)
├── Release_gnu_m/          # Compiled output (macOS)
└── Release_gnu_l/          # Compiled output (Linux)
```

## Building

### Building on Linux

#### Build Process

1. **Check dependencies**: The Makefile automatically verifies that gfortran and OpenBLAS are installed
2. **Auto-detect sources**: Automatically finds all `.f90` files in `src/` and `my_models/` directories
3. **Compile sources**: Compiles all detected source files
4. **Link**: Links against pre-compiled `libramses.a` from `modules_lin/` and OpenBLAS
5. **Output**: Creates `ramses.so` and `dynsim` in `Release_gnu_l/`

**Note**: The Makefile uses `wildcard` to automatically detect all `.f90` files in `my_models/`. You don't need to manually add new model files to the Makefile - just place them in `my_models/` and rebuild.

#### Makefile Targets

```bash
make -f Makefile.linux            # Build both ramses.so and dynsim (default)
make -f Makefile.linux dll        # Build only ramses.so (shared library)
make -f Makefile.linux exe        # Build only dynsim (executable)
make -f Makefile.linux clean      # Remove build artifacts
make -f Makefile.linux check-deps # Verify dependencies
make -f Makefile.linux help       # Show help
```

#### Output

After successful build, the following files will be in `Release_gnu_l/`:
```
Release_gnu_l/ramses.so   # Shared library for PyRAMSES
Release_gnu_l/dynsim      # Standalone executable
```

### Building on macOS (Apple Silicon)

`Makefile.macos` mirrors the Linux build: it auto-detects sources, links against
`libramses.a` from `modules_mac/`, and writes `Release_gnu_m/`. The shared
library is named `ramses.so` — the same as on Linux, because PyRAMSES keeps the
two apart in per-platform `libs/` folders rather than by filename. Do not rename
it to `ramses.dylib`.

```bash
make -f Makefile.macos                    # Build both ramses.so and dynsim (default)
make -f Makefile.macos dll                # Build only ramses.so (shared library)
make -f Makefile.macos exe                # Build only dynsim (executable)
make -f Makefile.macos clean              # Remove build artifacts
make -f Makefile.macos check-deps         # Verify dependencies
make -f Makefile.macos help               # Show help

make -f Makefile.macos BLAS=accelerate    # Use Apple Accelerate instead of OpenBLAS
make -f Makefile.macos FC=gfortran-15     # Use a versioned Homebrew gfortran
```

Two guards stop the build early with an explanatory message: the kit is arm64,
so a non-arm64 host is rejected up front; and `check-deps` refuses to continue
if your gfortran writes a different `.mod` format version than the kit.

#### Output

```
Release_gnu_m/ramses.so   # Shared library for PyRAMSES
Release_gnu_m/dynsim      # Standalone executable
```

### Building on Windows (MinGW / gfortran)

`Makefile.windows` is the free-toolchain Windows route and is independent of the
Intel Visual Studio projects — the two can coexist, writing to `Release_gnu_w/`
and `Release_intel_w64/` respectively. Run it from the **"MSYS2 MinGW 64-bit"**
shell; the plain MSYS shell is rejected because it produces a DLL linked to the
`msys-2.0` runtime that CPython cannot load.

```bash
make -f Makefile.windows               # Build both ramses.dll and dynsim.exe (default)
make -f Makefile.windows dll           # Build only ramses.dll (shared library)
make -f Makefile.windows exe           # Build only dynsim.exe (executable)
make -f Makefile.windows clean         # Remove build artifacts
make -f Makefile.windows check-deps    # Verify dependencies
make -f Makefile.windows help          # Show help
```

The Windows gfortran kit ships its archive as `libramses.lib` rather than
`libramses.a`. MinGW's linker does not probe that name for `-lramses`, so the
Makefile passes `modules_win_gfortran/libramses.lib` by explicit path.

#### Output

```
Release_gnu_w/ramses.dll   # Shared library for PyRAMSES
Release_gnu_w/dynsim.exe   # Standalone executable
```

The resulting DLL depends on the MSYS2 runtime DLLs (`libgfortran`, `libgcc_s`,
`libgomp`, `libopenblas`). Keep the MinGW `bin` directory on `PATH`, or copy
those DLLs next to `ramses.dll`, when loading it from Python.

### Building on Windows (Intel oneAPI)

#### Visual Studio Projects

The solution contains three projects:

1. **dllramses (ramses.dll)**
   - **Purpose**: Creates the main dynamic link library for PyRAMSES integration
   - **Output**: `ramses.dll` - Used by PyRAMSES to access your custom models
   - **Usage**: Primary project for Python integration on Windows

2. **exeramses (dynsim.exe)**
   - **Purpose**: Creates a standalone executable for direct simulation
   - **Output**: `dynsim.exe` - Command-line simulation tool
   - **Usage**: Run simulations directly without Python/Java interface
   - **Features**: Includes all your custom models for standalone operation

3. **MDL (ramsesmdl.dll)**
   - **Purpose**: Auxiliary project that compiles selected model files into a separate model library
   - **Output**: `ramsesmdl.dll`

#### Step-by-Step Build Process

1. **Open Solution**: Open `URAMSES.sln` in Microsoft Visual Studio
2. **Verify Compiler**: Ensure Intel Fortran compiler is properly configured
3. **Select Configuration**: Choose `Release|x64` configuration
4. **Build**: Right-click solution → "Build Solution"

#### Output Files

All compiled files will be created in `Release_intel_w64/`:
- `ramses.dll` - For PyRAMSES/STEPSS integration
- `dynsim.exe` - Standalone executable

## Which compiler do I need?

The `modules_*/` directories ship RAMSES pre-compiled. A gfortran `.mod` file
can only be read by the compiler generation that wrote it, so building your
own models needs a gfortran matching the kit for your platform.

Each kit records exactly what built it:

```sh
cat modules_lin/BUILDINFO.txt          # Linux
cat modules_mac/BUILDINFO.txt          # macOS arm64
cat modules_win_gfortran/BUILDINFO.txt # Windows MinGW
```

These files appear once the kit-sync CI has refreshed a kit at least once;
until then `check-deps` reports the kit as predating provenance tracking
instead of printing a compiler.

`check-deps` compares your compiler against the kit and fails early on a
mismatch, naming both module versions:

```sh
make -f Makefile.linux check-deps
```

If your distribution's default gfortran is the wrong generation, install a
matching one and pass it explicitly:

```sh
make -f Makefile.linux FC=gfortran-11
```

The same information appears in the toolchain table at the top of every
release. The Intel kit in `modules/` (used by `URAMSES.sln`) is maintained by
hand — see `modules/BUILDINFO.txt`.

## Adding Custom Models

### Step-by-Step Process

#### 1. Create Your Model File

Place your generated `.f90` model files (created by CODEGEN) into the `my_models/` directory.

**Linux, macOS, Windows/MinGW**: The Makefiles automatically detect and compile any `.f90` files in this directory.
**Windows/Intel**: You'll need to add the file to the Visual Studio project (see step 3).

**Example**: If you create `my_models/exc_MYMODEL.f90`, it will be automatically included in every Makefile build.

#### 2. Register Model in Association Files

Edit the appropriate association file in `src/` to register your models. This tells RAMSES which subroutine to call for your model.

**For Exciters** (`src/usr_exc_models.f90`):
```fortran
select case (modelname)
   case('YOUR_MODEL_NAME')
      exc_ptr => your_model_subroutine
end select
```

**For Injectors** (`src/usr_inj_models.f90`):
```fortran
select case (modelname)
   case('YOUR_MODEL_NAME')
      inj_ptr => your_model_subroutine
end select
```

**For Torque Models** (`src/usr_tor_models.f90`):
```fortran
select case (modelname)
   case('YOUR_MODEL_NAME')
      tor_ptr => your_model_subroutine
end select
```

**For Two-port Models** (`src/usr_twop_models.f90`):
```fortran
select case (modelname)
   case('YOUR_MODEL_NAME')
      twop_ptr => your_model_subroutine
end select
```

**For Discrete Control Models** (`src/usr_dctl_models.f90`):
```fortran
select case (modelname)
   case('YOUR_MODEL_NAME')
      dctl_ptr => your_model_subroutine
end select
```

**Important**: The `modelname` in the `case` statement must match exactly the name used in your simulation case files.

#### 3. Add to Visual Studio Project (Windows/Intel Only)

For Windows/Intel builds, you need to manually add the model file to the Visual Studio project:

1. Right-click on the `dllramses` project in Solution Explorer
2. Select "Add" → "Existing Item"
3. Navigate to `my_models/` and select your `.f90` files
4. Click "Add"

**Note**: For the Makefile builds (Linux, macOS, Windows/MinGW) this step is **not required**. Each Makefile detects all `.f90` files in `my_models/` using `wildcard`, so your new model will be compiled automatically.

#### 4. Rebuild

- **Linux**:
  ```bash
  make -f Makefile.linux clean all
  ```

- **macOS**:
  ```bash
  make -f Makefile.macos clean all
  ```

- **Windows/MinGW**:
  ```bash
  make -f Makefile.windows clean all
  ```

  Each Makefile will automatically compile your new model file(s) from `my_models/`.

- **Windows/Intel**: Rebuild the solution in Visual Studio (Build → Rebuild Solution)

### Model File Naming Conventions

While not strictly required, following naming conventions helps organization:
- **Exciters**: `exc_*.f90` (e.g., `exc_MYMODEL.f90`)
- **Injectors**: `inj_*.f90` (e.g., `inj_MYMODEL.f90`)
- **Torque**: `tor_*.f90` (e.g., `tor_MYMODEL.f90`)
- **Two-port**: `twop_*.f90` (e.g., `twop_MYMODEL.f90`)
- **Discrete Control**: `dctl_*.f90` (e.g., `dctl_MYMODEL.f90`)

## Using Your Models

### With PyRAMSES (Python)

```python
import pyramses

# Linux
ram = pyramses.sim('/path/to/your/URAMSES/Release_gnu_l')

# Windows
ram = pyramses.sim(r'C:\path\to\your\URAMSES\Release_intel_w64')

# Your models are now available for use in simulations
```

### With STEPSS (Java) - Windows Only

Use `ramses.dll` with the STEPSS Java interface — your custom models will be available in STEPSS simulations.

### Standalone Simulation

```bash
# Linux
cd Release_gnu_l
./dynsim

# Windows
cd Release_intel_w64
dynsim.exe
```

## Examples

The `my_models/` directory contains example models:
- `exc_ENTSOE_lim.f90`: ENTSO-E exciter model with limiters
- `inj_AIR_COND1_mod.f90`: Air conditioning load model

Models already provided by the RAMSES library in `modules_lin/` must not be
duplicated here — defining a subroutine that the library also exports produces
a duplicate-symbol link error. Register those by name in the `src/usr_*_models.f90`
routers instead, as `VFAULT` is in `src/usr_inj_models.f90`.

## Troubleshooting

### Linux Issues

1. **gfortran not found**
   ```bash
   # Ubuntu/Debian
   sudo apt install gfortran
   ```

2. **OpenBLAS not found**
   ```bash
   # Ubuntu/Debian
   sudo apt install libopenblas-dev
   ```

3. **Module files not found**
   - Ensure `modules_lin/` directory exists and contains `.mod` files
   - Verify `libramses.a` is present in `modules_lin/`

4. **Undefined reference errors**
   - Check that your model subroutine names match those in association files
   - Ensure all dependencies are properly linked
   - Verify that OpenBLAS is properly installed

5. **New model not being compiled**
   - Ensure the model file has a `.f90` extension
   - Check that the file is in the `my_models/` directory
   - Run `make -f Makefile.linux clean all` to force a full rebuild

### macOS Issues

1. **`Cannot read module file ... created by a different version of GNU Fortran`**
   - The `modules_mac/` kit is built with GFORTRAN module version 16 (gfortran 15+)
   - Run `make -f Makefile.macos check-deps` to see both versions side by side
   - Install a matching compiler with `brew install gcc`, then pass it explicitly
     if Homebrew only provides a versioned binary: `make -f Makefile.macos FC=gfortran-15`

2. **`modules_mac/ ships arm64 objects`**
   - macOS builds require Apple Silicon

3. **OpenBLAS not found**
   - `brew install openblas` (Homebrew keeps it keg-only; the Makefile finds it
     via `brew --prefix`)
   - Or skip it entirely with `make -f Makefile.macos BLAS=accelerate`

4. **`gfortran: command not found`**
   - macOS ships no Fortran compiler: `brew install gcc`

### Windows Issues (MinGW / gfortran)

1. **Makefile refuses to run in the MSYS shell**
   - Use the "MSYS2 MinGW 64-bit" shell. The plain MSYS shell links against the
     `msys-2.0` runtime and produces a DLL that CPython cannot load

2. **`Cannot read module file ... created by a different version of GNU Fortran`**
   - The `modules_win_gfortran/` kit needs gfortran 15+ (module version 16)
   - Update the toolchain: `pacman -Syu mingw-w64-x86_64-gcc-fortran`

3. **`cannot find -lramses`**
   - Expected: the kit ships `libramses.lib`, which MinGW's linker does not probe
     for `-lramses`. The Makefile passes it by explicit path — do not replace that
     with `-lramses`

4. **Python cannot load `ramses.dll`**
   - The DLL needs the MSYS2 runtime (`libgfortran`, `libgcc_s`, `libgomp`,
     `libopenblas`). Add the MinGW `bin` directory to `PATH` or copy those DLLs
     alongside `ramses.dll`

### Windows Issues (Intel oneAPI)

1. **Compilation Errors**: Ensure Intel Fortran compiler is properly installed and configured
2. **Missing Models**: Verify model names match exactly in association files
3. **DLL Loading**: Check that the path to `ramses.dll` is correct in PyRAMSES
4. **Model Parameters**: Ensure parameter files are properly formatted

### Debug Tips

- Check compiler output for compilation errors
- Verify model subroutine names match exactly in association files
- Test with simple models first before complex implementations
- On Linux, use `ldd Release_gnu_l/ramses.so` to check library dependencies

## Platform Comparison

| Feature | Linux | macOS (arm64) | Windows (MinGW) | Windows (Intel) |
|---------|-------|---------------|-----------------|-----------------|
| Compiler | gfortran | gfortran (Homebrew) | gfortran (MinGW-w64) | Intel Fortran |
| BLAS Library | OpenBLAS | OpenBLAS or Accelerate | OpenBLAS | Intel MKL |
| Build System | `Makefile.linux` | `Makefile.macos` | `Makefile.windows` | Visual Studio |
| Output Library | `ramses.so` | `ramses.so` | `ramses.dll` | `ramses.dll` |
| Output Executable | `dynsim` | `dynsim` | `dynsim.exe` | `dynsim.exe` |
| Output Directory | `Release_gnu_l/` | `Release_gnu_m/` | `Release_gnu_w/` | `Release_intel_w64/` |
| Module Directory | `modules_lin/` | `modules_mac/` | `modules_win_gfortran/` | `modules/` |
| Model Auto-Detection | ✅ Automatic (wildcard) | ✅ Automatic (wildcard) | ✅ Automatic (wildcard) | ❌ Manual (VS project) |
| Host Prerequisites | gfortran, OpenBLAS | Homebrew gcc, OpenBLAS | MSYS2 toolchain, OpenBLAS | VS, Intel Fortran |

## Documentation

| Document | Description |
|----------|-------------|
| [URAMSES developer guide](https://stepss.sps-lab.org/developer/uramses/) | Building URAMSES and registering user models |
| [PyRAMSES documentation](https://stepss.sps-lab.org/pyramses/) | Python interface that loads the compiled library |
| [Installing the Intel oneAPI Fortran compiler.pdf](Installing%20the%20Intel%20oneAPI%20Fortran%20compiler.pdf) | Windows compiler setup (local PDF) |

For issues and questions, contact the Sustainable Power Systems Lab (SPS-L) at [info@sps-lab.org](mailto:info@sps-lab.org).

## License

URAMSES is distributed under the **Apache License 2.0** — see [LICENSE.rst](LICENSE.rst). Copyright © Petros Aristidou.

Note that RAMSES itself, which URAMSES links against (the pre-compiled `libramses` in `modules/` and `modules_lin/`), is **not** covered by this licence — see [NOTICE](NOTICE) for the terms applying to the RAMSES, PFC and CODEGEN components.

## Authors

Developed and maintained by the [Sustainable Power Systems Laboratory (SPS-L)](https://sps-lab.org/) at the Cyprus University of Technology, under the direction of Dr. Petros Aristidou.

RAMSES, the simulator URAMSES links against, builds on the original work of Dr. Thierry Van Cutsem (University of Liège) — see [NOTICE](NOTICE).
