# OpenFOAM Chemistry Model Utilities

## chemkinToFoam.py

Converts a CHEMKIN-III mechanism (kinetics + thermo + transport) to OpenFOAM format by
driving the real, upstream `chemkinToFoam` utility instead of re-parsing CHEMKIN syntax
by hand. Produces correct reaction types (reversible / third-body / Troe / Lindemann
fall-off), real per-specie NASA thermo data, and real per-specie Sutherland transport
coefficients (`As`, `Ts`) fitted from the mechanism's own Lennard-Jones transport data via
Cantera.

Two things about real-world CHEMKIN files that trip up OpenFOAM's `chemkinToFoam` are
handled automatically:

- **PLOG** (pressure-dependent Arrhenius) auxiliary lines aren't supported by OpenFOAM's
  chemkinReader at all ("unknown third-body specie PLOG"). They're stripped, so those
  reactions fall back to their own leading (high-pressure-limit) Arrhenius parameters —
  pressure-dependence is lost for the affected reactions.
- The `THERMO`/`THERMO ALL` temperature-range line must be an exact fixed-width 3×10
  column record (Fortran `3F10.0`); free-format thermo files are reformatted to match.

`chemkinToFoam` also doesn't read a CHEMKIN `.TRAN` transport file directly — its
`<CHEMKINTransport>` argument must already be an OpenFOAM dictionary of per-specie
`sutherlandTransport` coefficients, which this script computes.

### Requirements

- `cantera` (Python package) — used for mechanism parsing and transport property
  evaluation.
- The OpenFOAM `chemkinToFoam` binary on `$PATH` (source an OpenFOAM environment first).

### Usage

```bash
./chemkinToFoam.py <kinetics.CKI> <thermo.CKT> <transport.TRAN> <output_dir> \
  [--chemkinToFoam PATH] [--n-fit-points N]
```

### Example

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
./utilities/chemkinToFoam.py \
  external/Kinetic-Mechanisms/Gas-Phase/CoreMechanism_C0-C4/Soot-NOx/C1_C3_HT_NOX_159_2459/kinetics.CHEMKIN.CKI \
  external/Kinetic-Mechanisms/Gas-Phase/CoreMechanism_C0-C4/Soot-NOx/C1_C3_HT_NOX_159_2459/thermo.CHEMKIN.CKT \
  external/Kinetic-Mechanisms/Gas-Phase/CoreMechanism_C0-C4/Soot-NOx/C1_C3_HT_NOX_159_2459/TOT2003.TRAN \
  tutorials/DLBFoam-Hydrogen-Tutorials/mech/thermo/
```

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
