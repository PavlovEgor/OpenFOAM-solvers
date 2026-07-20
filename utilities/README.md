# OpenFOAM Chemistry Model Utilities

## chemkin2foam.py

Converts CHEMKIN mechanism files to OpenFOAM format for use with chemistryModel.

### Features

- Parses CHEMKIN thermodynamics (.CKT) and kinetics (.CKI) files
- Extracts NASA 7-coefficient polynomial data
- Generates OpenFOAM-format files:
  - `thermo.compressibleGas` — species thermodynamic properties
  - `reactions` — reaction definitions
- No external dependencies (built-in CHEMKIN parser)

### Usage

```bash
./chemkin2foam.py <thermo_file.CKT> <kinetics_file.CKI> [output_dir]
```

### Arguments

- `thermo_file.CKT`: CHEMKIN thermodynamic database file
- `kinetics_file.CKI`: CHEMKIN kinetics/reactions file
- `output_dir` (optional): Output directory (default: current directory)

### Example

Convert C1_C3_HT_NOX mechanism:

```bash
./utilities/chemkin2foam.py \
  external/Kinetic-Mechanisms/Gas-Phase/CoreMechanism_C0-C4/Soot-NOx/C1_C3_HT_NOX_159_2459/thermo.CHEMKIN.CKT \
  external/Kinetic-Mechanisms/Gas-Phase/CoreMechanism_C0-C4/Soot-NOx/C1_C3_HT_NOX_159_2459/kinetics.CHEMKIN.CKI \
  tutorials/DLBFoam-Hydrogen-Tutorials/mech/thermo/
```

### Output

- `thermo.compressibleGas`: Contains species definitions with:
  - Molecular weight
  - Elemental composition
  - NASA polynomials (high and low temperature ranges)
  - Transport properties
  
- `reactions`: Contains reaction definitions with:
  - Element list
  - Species list
  - Reaction equations and Arrhenius parameters (A, β, Tₐ)

### Supported Formats

- NASA 7-coefficient polynomial format (most common in CHEMKIN)
- Arrhenius and reversible reaction rate expressions

## pyjac4foam_options_gen.py

Generates OpenFOAM Make files (`Make/files`, `Make/filesCuda`, `Make/options`) for any pyJac-generated chemistry mechanism.

### Usage

```bash
./pyjac4foam_options_gen.py <mechanism_path> <output_dir>
```

### Arguments

- `mechanism_path`: Path to the directory containing pyJac-generated .cu files (typically `pyJacSource/<mechanism_name>/out/`)
- `output_dir`: Output directory where `Make/files`, `Make/filesCuda`, and `Make/options` will be generated (typically `src/<solver>/Make/`)

### Example

Generate Make files for the C1_C3_HT_NOX mechanism in kodesChemistryModel:

```bash
./utilities/pyjac4foam_options_gen.py \
  src/kodesChemistyModel/pyJacSource/C1_C3_HT_NOX_159_2459/out \
  src/kodesChemistyModel/Make/
```

This will:
1. Scan the mechanism directory for all `.cu` files (recursively)
2. Identify include directories (for `-I` flags)
3. Generate Make files with:
   - KODES library paths
   - Mechanism source files
   - Include paths for all generated subdirectories
   - OpenFOAM library dependencies

### Output Files

- `Make/files` — lists of CUDA source files to compile
- `Make/filesCuda` — identical to files (for consistency with OpenFOAM build system)
- `Make/options` — compiler flags, include paths, and library links

### How it Works

The script:
1. Recursively finds all `.cu` files in the mechanism directory
2. Computes relative paths from the output directory to mechanism files
3. Identifies all subdirectories that should be included via `-I` flags
4. Adds standard OpenFOAM library includes and KODES paths
5. Generates wmake-compatible Make files with continuation backslashes
