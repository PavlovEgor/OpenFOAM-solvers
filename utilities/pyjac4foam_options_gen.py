#!/usr/bin/env python3
"""
Generate OpenFOAM Make/files and Make/options for pyJac-generated chemistry mechanisms.

Usage:
    pyjac4foam_options_gen.py <mechanism_path> <output_dir>

Example:
    pyjac4foam_options_gen.py pyJacSource/C1_C3_HT_NOX_159_2459/out src/kodesChemistyModel/Make/
"""

import sys
import os
from pathlib import Path
from collections import defaultdict


def find_cu_files(mechanism_path):
    """Recursively find all .cu files and group by directory."""
    mechanism_path = Path(mechanism_path)
    mechanism_abs = mechanism_path.resolve()

    if not mechanism_abs.exists():
        raise FileNotFoundError(f"Mechanism path not found: {mechanism_abs}")

    files_by_dir = defaultdict(list)

    for cu_file in sorted(mechanism_abs.rglob("*.cu")):
        # Get relative path from mechanism_path
        rel_path = cu_file.relative_to(mechanism_abs)
        subdir = rel_path.parent
        files_by_dir[str(subdir)].append(str(rel_path))

    return files_by_dir


def get_include_dirs(mechanism_path):
    """Get all directories that should be included."""
    mechanism_abs = Path(mechanism_path).resolve()
    include_dirs = set()

    # Add root and all subdirectories
    include_dirs.add(".")
    for item in mechanism_abs.rglob("*"):
        if item.is_dir() and item != mechanism_abs:
            rel_path = item.relative_to(mechanism_abs)
            include_dirs.add(str(rel_path))

    return sorted(include_dirs)


def compute_relative_path(mechanism_path, output_dir):
    """Compute relative path from output_dir's parent to mechanism_path."""
    import os

    # Get absolute paths
    mechanism_abs = os.path.abspath(mechanism_path)
    output_abs = os.path.abspath(output_dir)

    # If output_abs ends with /, remove it for consistent comparison
    if output_abs.endswith(os.sep):
        output_abs = output_abs[:-1]

    # Compute from parent of output_dir (typically src/kodesChemistyModel/ if output_dir is Make/)
    output_parent = os.path.dirname(output_abs)

    # Use os.path.relpath to compute relative path
    try:
        rel_path = os.path.relpath(mechanism_abs, output_parent)
        return rel_path
    except ValueError:
        # Different drives on Windows, fall back to absolute
        return mechanism_abs


def generate_files_content(mechanism_path, output_dir, files_by_dir):
    """Generate content for Make/files"""
    rel_base = compute_relative_path(mechanism_path, output_dir)

    lines = []

    # Add KODES source files first
    kodes_files = [
        "../../external/KODES/src/basic_linalg.cu",
        "../../external/KODES/src/Resources/SeulexDeviceResources.cu",
        "../../external/KODES/src/Resources/DeviceResources.cu",
        "../../external/KODES/src/Resources/HostResources.cu",
        "../../external/KODES/src/ODESystem/pyJacSystem.cu",
    ]
    lines.extend(kodes_files)

    # Add mechanism files
    for subdir in sorted(files_by_dir.keys()):
        for filename in files_by_dir[subdir]:
            if subdir == ".":
                file_path = f"{rel_base}/{filename}"
            else:
                file_path = f"{rel_base}/{filename}"
            lines.append(file_path)

    # Add object files
    lines.extend([
        "makeKodesChemistryModelTypes.C",
        "makeKodesChemistrySolverTypes.C",
        "",
        "LIB = $(FOAM_USER_LIBBIN)/libkodesChemistryModel",
    ])

    return "\n".join(lines) + "\n"


def generate_filesCuda_content(mechanism_path, output_dir, files_by_dir):
    """Generate content for Make/filesCuda (same as Make/files for CUDA files)"""
    return generate_files_content(mechanism_path, output_dir, files_by_dir)


def generate_options_content(mechanism_path, output_dir, include_dirs):
    """Generate content for Make/options (include paths only)"""
    rel_base = compute_relative_path(mechanism_path, output_dir)

    lines = [
        "EXE_INC = \\",
        "    -I$(LIB_SRC)/finiteVolume/lnInclude \\",
        "    -I$(LIB_SRC)/meshTools/lnInclude \\",
        "    -I$(LIB_SRC)/ODE/lnInclude \\",
        "    -I$(LIB_SRC)/transportModels/compressible/lnInclude \\",
        "    -I$(LIB_SRC)/thermophysicalModels/reactionThermo/lnInclude \\",
        "    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \\",
        "    -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \\",
        "    -I$(LIB_SRC)/thermophysicalModels/chemistryModel/lnInclude \\",
        "    -I$(LIB_SRC)/thermophysicalModels/thermophysicalProperties/lnInclude \\",
        "    -I../../external/KODES/include \\",
        "    -I../../external/KODES/include/Integrators \\",
        "    -I../../external/KODES/src/Integrators \\",
        "    -I../../external/KODES/include/ODESystem \\",
        "    -I../../external/KODES/include/Resources \\",
    ]

    # Add mechanism includes
    for inc_dir in include_dirs:
        if inc_dir == ".":
            path = f"-I{rel_base}"
        else:
            path = f"-I{rel_base}/{inc_dir}"
        lines.append(f"    {path} \\")

    lines.append("")
    lines.extend([
        "LIB_LIBS = \\",
        "    -lfiniteVolume \\",
        "    -lmeshTools \\",
        "    -lODE \\",
        "    -lcompressibleTransportModels \\",
        "    -lfluidThermophysicalModels \\",
        "    -lreactionThermophysicalModels \\",
        "    -lspecie \\",
        "    -lchemistryModel \\",
        "    -lthermophysicalProperties ",
        "",
        "",
        "WMAKE_LOCAL=$(MAKE_DIR)/../../../wmake",
        "include $(WMAKE_LOCAL)/c++",
    ])

    return "\n".join(lines) + "\n"


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    mechanism_path = sys.argv[1]
    output_dir = sys.argv[2]

    try:
        # Find all .cu files
        files_by_dir = find_cu_files(mechanism_path)

        if not files_by_dir:
            print(f"Warning: No .cu files found in {mechanism_path}")

        # Get include directories
        include_dirs = get_include_dirs(mechanism_path)

        # Ensure output directory exists
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Generate and write Make/files
        files_content = generate_files_content(mechanism_path, output_dir, files_by_dir)
        files_path = output_path / "files"
        files_path.write_text(files_content)
        print(f"✓ Generated {files_path}")

        # Generate and write Make/filesCuda
        filesCuda_content = generate_filesCuda_content(mechanism_path, output_dir, files_by_dir)
        filesCuda_path = output_path / "filesCuda"
        filesCuda_path.write_text(filesCuda_content)
        print(f"✓ Generated {filesCuda_path}")

        # Generate and write Make/options
        options_content = generate_options_content(mechanism_path, output_dir, include_dirs)
        options_path = output_path / "options"
        options_path.write_text(options_content)
        print(f"✓ Generated {options_path}")

        print(f"\nSummary:")
        print(f"  Mechanism source: {mechanism_path}")
        print(f"  Output directory: {output_dir}")
        print(f"  .cu files found: {sum(len(files) for files in files_by_dir.values())}")
        print(f"  Include directories: {len(include_dirs)}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
