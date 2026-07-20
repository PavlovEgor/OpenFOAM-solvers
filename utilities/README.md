# OpenFOAM Chemistry Model Utilities

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
