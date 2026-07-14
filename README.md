# URAMSES - User Models for PyRAMSES

## About

URAMSES is a project that enables the integration of custom user models into PyRAMSES (Python interface for RAMSES power system simulator) and STEPSS. This repository provides the framework and tools needed to compile and link your own Fortran models with the simulation environment.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Project Structure](#project-structure)
- [Model Types](#model-types)
- [Quick Start](#quick-start)
- [Building](#building)
  - [Linux](#building-on-linux)
  - [Windows](#building-on-windows)
  - [Docker](#building-with-docker)
- [Adding Custom Models](#adding-custom-models)
- [Using Your Models](#using-your-models)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Platform Comparison](#platform-comparison)
- [Documentation & Support](#documentation--support)
- [License](#license)

## Prerequisites

### Windows
- **Microsoft Visual Studio** (2019 or later recommended)
- **Intel oneAPI Fortran Compiler** (formerly Intel Fortran)
- **PyRAMSES** (Python package) or **STEPSS** (Java package)

For detailed installation instructions of the Intel oneAPI Fortran compiler, refer to the included PDF:
[Installing the Intel oneAPI Fortran compiler.pdf](Installing%20the%20Intel%20oneAPI%20Fortran%20compiler.pdf)

### Linux
- **gfortran** (GNU Fortran compiler)
- **OpenBLAS** (optimized BLAS library)
- **PyRAMSES** (Python package)

Install dependencies on your Linux distribution:

```bash
# Ubuntu/Debian
sudo apt install gfortran libopenblas-dev

# Fedora/RHEL
sudo dnf install gcc-gfortran openblas-devel

# Arch Linux
sudo pacman -S gcc-fortran openblas
```

### Docker (any platform)
- **Docker** with Docker Compose (v2)
- No compiler or library installation required

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
│   ├── twop_*.f90          # Two-port models
│   └── *.txt               # Model parameter files
├── modules/                # Pre-compiled modules (Windows/Intel Fortran)
│   ├── *.mod               # Module interface files
│   └── libramses.lib       # Pre-compiled RAMSES library
├── modules_lin/            # Pre-compiled modules (Linux/gfortran)
│   ├── *.mod               # Module interface files
│   └── libramses.a         # Pre-compiled RAMSES library
├── URAMSES.sln             # Visual Studio solution file (Windows)
├── dllramses.vfproj        # DLL project - ramses.dll (Windows)
├── exeramses.vfproj        # Executable project - dynsim.exe (Windows)
├── MDL.vfproj              # Model library project (Windows)
├── Makefile.gfortran       # Makefile for Linux builds
├── docker/                 # Docker build environment
│   └── Dockerfile          # Ubuntu 24.04 + gfortran + OpenBLAS
├── docker-compose.yml      # One-command Docker build
├── build.sh                # Docker build helper script
├── output/                 # Docker build output (ramses.so)
├── Release_intel_w64/      # Compiled output (Windows)
└── Release_gnu_l/          # Compiled output (Linux)
```

## Model Types

URAMSES supports several types of power system models:

- **Exciters (`exc_*`)**: Generator excitation system models
- **Injectors (`inj_*`)**: Current/voltage injection models for faults/disturbances
- **Torque (`tor_*`)**: Mechanical torque models for generators
- **Two-port (`twop_*`)**: Two-port network models (e.g., SVC, STATCOM)
- **Discrete Control (`dctl_*`)**: Discrete control system models

## Quick Start

### Linux
```bash
# Build everything (shared library + executable)
make -f Makefile.gfortran

# Run standalone simulation
./Release_gnu_l/dynsim
```

### Windows
1. Open `URAMSES.sln` in Visual Studio
2. Select `Release|x64` configuration
3. Build → Build Solution
4. Run `Release_intel_w64\dynsim.exe`

### Docker
```bash
docker compose build              # One-time image build
docker compose run --rm uramses-build  # Build ramses.so → output/ramses.so
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
make -f Makefile.gfortran            # Build both ramses.so and dynsim (default)
make -f Makefile.gfortran dll        # Build only ramses.so (shared library)
make -f Makefile.gfortran exe        # Build only dynsim (executable)
make -f Makefile.gfortran clean      # Remove build artifacts
make -f Makefile.gfortran check-deps # Verify dependencies
make -f Makefile.gfortran help       # Show help
```

#### Output

After successful build, the following files will be in `Release_gnu_l/`:
```
Release_gnu_l/ramses.so   # Shared library for PyRAMSES
Release_gnu_l/dynsim      # Standalone executable
```

### Building on Windows

#### Visual Studio Projects

The solution contains three main projects:

1. **dllramses (ramses.dll)**
   - **Purpose**: Creates the main dynamic link library for PyRAMSES integration
   - **Output**: `ramses.dll` - Used by PyRAMSES to access your custom models
   - **Usage**: Primary project for Python integration on Windows

2. **exeramses (dynsim.exe)**
   - **Purpose**: Creates a standalone executable for direct simulation
   - **Output**: `dynsim.exe` - Command-line simulation tool
   - **Usage**: Run simulations directly without Python/Java interface
   - **Features**: Includes all your custom models for standalone operation

#### Step-by-Step Build Process

1. **Open Solution**: Open `URAMSES.sln` in Microsoft Visual Studio
2. **Verify Compiler**: Ensure Intel Fortran compiler is properly configured
3. **Select Configuration**: Choose `Release|x64` configuration
4. **Build**: Right-click solution → "Build Solution"

#### Output Files

All compiled files will be created in `Release_intel_w64/`:
- `ramses.dll` - For PyRAMSES/STEPSS integration
- `dynsim.exe` - Standalone executable

### Building with Docker

Docker provides a self-contained build environment — no compiler or library installation required on the host.

#### Project Files

```
docker/Dockerfile       # Build environment image (Ubuntu 24.04, gfortran, OpenBLAS)
docker-compose.yml      # One-command build configuration
build.sh                # Convenience wrapper script
```

#### One-Time Setup

```bash
docker compose build    # Builds the image (~2 min)
```

#### Build ramses.so

```bash
# Option 1: docker compose
docker compose run --rm uramses-build

# Option 2: helper script
./build.sh
```

The compiled `ramses.so` is written to the `output/` directory on the host.

#### How It Works

- The repository root is bind-mounted into the container at `/uramses`
- The container runs `make -f Makefile.gfortran dll` and copies the result to `/output`
- Edits to `my_models/` on the host are visible immediately — no image rebuild needed
- Only the shared library (`ramses.so`) is built, not the standalone executable

#### Adding Models with Docker

1. Create or edit model files in `my_models/` on the host
2. Register in the appropriate `src/usr_*_models.f90`
3. Run `docker compose run --rm uramses-build` (or `./build.sh`)
4. Copy `output/ramses.so` to your PyRAMSES/STEPSS environment

## Adding Custom Models

### Step-by-Step Process

#### 1. Create Your Model File

Place your generated `.f90` model files (created by CODEGEN) into the `my_models/` directory. 

**Linux**: The Makefile will automatically detect and compile any `.f90` files in this directory.  
**Windows**: You'll need to add the file to the Visual Studio project (see step 3).

**Example**: If you create `my_models/exc_MYMODEL.f90`, it will be automatically included in Linux builds.

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

#### 3. Add to Visual Studio Project (Windows Only)

For Windows builds, you need to manually add the model file to the Visual Studio project:

1. Right-click on the `dllramses` project in Solution Explorer
2. Select "Add" → "Existing Item"
3. Navigate to `my_models/` and select your `.f90` files
4. Click "Add"

**Note**: For Linux builds, this step is **not required**. The Makefile automatically detects all `.f90` files in `my_models/` using `wildcard`, so your new model will be compiled automatically.

#### 4. Rebuild

- **Linux**: 
  ```bash
  make -f Makefile.gfortran clean all
  ```
  The Makefile will automatically compile your new model file(s) from `my_models/`.

- **Windows**: Rebuild the solution in Visual Studio (Build → Rebuild Solution)

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

```java
// Use ramses.dll with STEPSS Java interface
// Your custom models will be available in STEPSS simulations
```

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

The `my_models/` directory contains several example models:
- `exc_ENTSOE_lim.f90`: ENTSO-E exciter model with limiters
- `exc_GENERIC3.f90`: Generic exciter model type 3
- `exc_GENERIC4.f90`: Generic exciter model type 4
- `exc_ST1A.f90`: IEEE ST1A exciter model
- `inj_AIR_COND1_mod.f90`: Air conditioning load model
- `tor_ENTSOE_simp.f90`: Simplified ENTSO-E torque model

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
   - Run `make -f Makefile.gfortran clean all` to force a full rebuild

### Windows Issues

1. **Compilation Errors**: Ensure Intel Fortran compiler is properly installed and configured
2. **Missing Models**: Verify model names match exactly in association files
3. **DLL Loading**: Check that the path to `ramses.dll` is correct in PyRAMSES
4. **Model Parameters**: Ensure parameter files are properly formatted

### Docker Issues

1. **`docker compose` not found**
   - Ensure Docker Desktop or Docker Engine with Compose v2 is installed
   - Older installations may need `docker-compose` (with hyphen) instead

2. **Permission denied on `output/`**
   - The `output/` directory is created by the container. If permissions are wrong:
     ```bash
     sudo chown -R $(id -u):$(id -g) output/
     ```

3. **Stale build artifacts**
   - The container bind-mounts the repo, so old object files may persist:
     ```bash
     docker compose run --rm uramses-build make -f Makefile.gfortran clean dll
     ```

### Debug Tips

- Check compiler output for compilation errors
- Verify model subroutine names match exactly in association files
- Test with simple models first before complex implementations
- On Linux, use `ldd Release_gnu_l/ramses.so` to check library dependencies

## Platform Comparison

| Feature | Linux | Windows | Docker |
|---------|-------|---------|--------|
| Compiler | gfortran | Intel Fortran | gfortran (in container) |
| BLAS Library | OpenBLAS | Intel MKL | OpenBLAS (in container) |
| Build System | Makefile | Visual Studio | Docker Compose |
| Output Library | `ramses.so` | `ramses.dll` | `ramses.so` |
| Output Executable | `dynsim` | `dynsim.exe` | N/A |
| Output Directory | `Release_gnu_l/` | `Release_intel_w64/` | `output/` |
| Module Directory | `modules_lin/` | `modules/` | `modules_lin/` |
| Model Auto-Detection | ✅ Automatic (wildcard) | ❌ Manual (VS project) | ✅ Automatic (wildcard) |
| Host Prerequisites | gfortran, OpenBLAS | VS, Intel Fortran | Docker only |

## Documentation & Support

### Documentation

For comprehensive PyRAMSES documentation, visit: [https://stepss.sps-lab.org/pyramses/](https://stepss.sps-lab.org/pyramses/)

### Support

For issues and questions, contact:
- **Sustainable Power Systems Lab (SPS-L)**
- **Website:** [https://sps-lab.org](https://sps-lab.org)
- **Email:** info@sps-lab.org

## License

This project is licensed under the Apache License 2.0. See [LICENSE.rst](LICENSE.rst) for details.

Note that RAMSES itself, which URAMSES links against, is **not** covered by this licence — see [NOTICE](NOTICE) for the terms applying to the RAMSES, PFC and CODEGEN components.

---

**Last Updated:** April 2026