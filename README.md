# URAMSES - User Models for PyRAMSES

## About

URAMSES is a project that enables the integration of custom user models into PyRAMSES (Python interface for RAMSES power system simulator) and STEPSS. This repository provides the framework and tools needed to compile and link your own Fortran models with the simulation environment.

## Prerequisites

### Required Software
- **Microsoft Visual Studio** (2019 or later recommended)
- **Intel oneAPI Fortran Compiler** (formerly Intel Fortran)
- **PyRAMSES** (Python package) or **STEPSS## (java package)

### Installation Guide
For detailed installation instructions of the Intel oneAPI Fortran compiler, refer to the included PDF:
[Installing the Intel oneAPI Fortran compiler.pdf](Installing%20the%20Intel%20oneAPI%20Fortran%20compiler.pdf)

## Project Structure

```
URAMSES/
├── src/                    # Source code files
│   ├── c_interface.f90     # C interface for Python integration
│   ├── main.f90           # Main entry point
│   ├── usr_exc_models.f90 # Exciter model associations
│   ├── usr_inj_models.f90 # Injector model associations
│   ├── usr_tor_models.f90 # Torque model associations
│   ├── usr_twop_models.f90 # Two-port model associations
│   └── usr_dctl_models.f90 # Discrete control model associations
├── my_models/             # Your custom models go here
│   ├── exc_*.f90          # Exciter models
│   ├── inj_*.f90          # Injector models
│   ├── tor_*.f90          # Torque models
│   ├── twop_*.f90         # Two-port models
│   └── *.txt              # Model parameter files
├── modules/               # Compiled Fortran modules
├── URAMSES.sln           # Visual Studio solution file
├── dllramses.vfproj      # Main DLL project (ramses.dll)
├── exeramses.vfproj      # Executable project (dynsim.exe)
├── MDL.vfproj           # Model library project (ramsesmdl.dll)
└── Release_intel_w64/    # Compiled output directory
```

## Visual Studio Projects

The solution contains three main projects:

### 1. dllramses (ramses.dll)
- **Purpose**: Creates the main dynamic link library for PyRAMSES integration
- **Output**: `ramses.dll` - Used by PyRAMSES to access your custom models
- **Usage**: Primary project for Python integration

### 2. exeramses (dynsim.exe)
- **Purpose**: Creates a standalone executable for direct simulation
- **Output**: `dynsim.exe` - Command-line simulation tool
- **Usage**: Run simulations directly without Java interface
- **Features**: Includes all your custom models for standalone operation

### 3. MDL (ramsesmdl.dll)
- **Purpose**: Creates a model library DLL
- **Output**: `ramsesmdl.dll` - Additional model library
- **Usage**: Supplementary model library for advanced use cases

## Model Types

URAMSES supports several types of power system models:

- **Exciters (`exc_*`)**: Generator excitation system models
- **Injectors (`inj_*`)**: Current/voltage injection models for faults/disturbances
- **Torque (`tor_*`)**: Mechanical torque models for generators
- **Two-port (`twop_*`)**: Two-port network models (e.g., SVC, STATCOM)
- **Discrete Control (`dctl_*`)**: Discrete control system models

## Step-by-Step Integration Process

### 1. Download and Setup
```bash
# Download the latest release from:
# https://github.com/SPS-L/URAMSES/releases/
# Extract to your desired location
```

### 2. Add Your Models
Place your generated `.f90` model files (created by CODEGEN) into the `my_models/` directory. Each model should have:
- A Fortran source file (`.f90`)
- An optional parameter file (`.txt`) defining model parameters

### 3. Enable Models in Association Files
Edit the appropriate association file in `src/` to register your models:

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

### 4. Open in Visual Studio
1. Open `URAMSES.sln` in Microsoft Visual Studio
2. Ensure the Intel Fortran compiler is properly configured

### 5. Add Model Files to Project
1. Right-click on the `dllramses` project in Solution Explorer
2. Select "Add" → "Existing Item"
3. Navigate to `my_models/` and select your `.f90` files
4. Click "Add" to include them in the project

### 6. Compile the Projects
You can compile individual projects or the entire solution:

**Option A: Compile all projects (recommended)**
1. Right-click on the `URAMSES` solution in Solution Explorer
2. Select "Build Solution"
3. This will create:
   - `ramses.dll` (from dllramses project)
   - `dynsim.exe` (from exeramses project)
   - `ramsesmdl.dll` (from MDL project)

**Option B: Compile individual projects**
- **For PyRAMSES integration**: Right-click `dllramses` → "Build"
- **For STEPSS integration**: Right-click `exeramses` → "Build"
- **For model library**: Right-click `MDL` → "Build"

All compiled files will be created in `Release_intel_w64/`

### 7. Use Your Models

**Option A: With PyRAMSES (Python)**
```python
import pyramses

# Initialize simulation with your custom DLL
ram = pyramses.sim(r'C:\path\to\your\URAMSES\Release_intel_w64')

# Your models are now available for use in simulations
```

**Option B: With STEPSS (Java)**
```java
// Use the same ramses.dll with STEPSS Java interface
// Your custom models will be available in STEPSS simulations
```

**Option C: Standalone Simulation**
```bash
# Run simulations directly using the compiled executable
cd Release_intel_w64
./dynsim.exe
```
## Troubleshooting

### Common Issues
1. **Compilation Errors**: Ensure Intel Fortran compiler is properly installed and configured
2. **Missing Models**: Verify model names match exactly in association files
3. **DLL Loading**: Check that the path to `ramses.dll` is correct in PyRAMSES
4. **Model Parameters**: Ensure parameter files are properly formatted

### Debug Tips
- Check Visual Studio output window for compilation errors
- Verify model subroutine names match exactly in association files
- Test with simple models first before complex implementations

## Examples

The `my_models/` directory contains several example models:
- `exc_ENTSOE_lim.f90`: ENTSO-E exciter model with limiters
- `exc_GENERIC3.f90`: Generic exciter model type 3
- `inj_AIR_COND1_mod.f90`: Air conditioning load model
- `tor_ENTSOE_simp.f90`: Simplified ENTSO-E torque model

## Documentation

For comprehensive PyRAMSES documentation, visit:
[https://pyramses.sps-lab.org](https://pyramses.sps-lab.org)

## License

This project is licensed under the Academic Public License. See [LICENSE.rst](LICENSE.rst) for details.

## Support

For issues and questions:
- Check the PyRAMSES documentation
- Review example models in `my_models/`
- Ensure all prerequisites are properly installed

## Contributing

When contributing models:
1. Follow the existing naming conventions
2. Include parameter files (`.txt`) for your models
3. Test thoroughly before submission
4. Document any special requirements or dependencies 