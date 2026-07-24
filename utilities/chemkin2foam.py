#!/usr/bin/env python3
"""
Convert CHEMKIN mechanism files to OpenFOAM format.

Creates OpenFOAM thermodynamic and reaction files suitable for
use with chemistryModel in OpenFOAM.

Uses a direct CHEMKIN parser (no Cantera dependency for parsing).

Usage:
    chemkin2foam.py <chemkin_thermo> <chemkin_kinetics> [output_dir]

Example:
    chemkin2foam.py \
        external/Kinetic-Mechanisms/.../thermo.CHEMKIN.CKT \
        external/Kinetic-Mechanisms/.../kinetics.CHEMKIN.CKI \
        tutorials/DLBFoam-Hydrogen-Tutorials/mech/thermo/
"""

import sys
import re
from pathlib import Path
from collections import defaultdict


_TERM_RE = re.compile(r'^(\d+\.?\d*)?(.+)$')


def _format_reaction_side(side):
    """Insert spaces around '+' and between stoich coefficients and species
    (e.g. "2CH3O2+O2" -> "2 CH3O2 + O2") so OpenFOAM's reaction parser,
    which tokenizes on whitespace, can read coefficients and species names."""
    terms = []
    for term in side.split('+'):
        term = term.strip()
        coeff, specie = _TERM_RE.match(term).groups()
        if coeff and coeff not in ('1', '1.0'):
            terms.append(f"{coeff} {specie}")
        else:
            terms.append(specie)
    return ' + '.join(terms)


class ChemkinParser:
    """Simple CHEMKIN file parser."""

    def __init__(self):
        self.elements = []
        self.species = {}
        self.species_thermo = {}
        self.reactions = []

    def parse_thermo(self, filename):
        """Parse CHEMKIN thermo file."""
        with open(filename, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i]

            if not line.strip() or line.strip().startswith('!'):
                i += 1
                continue

            # Skip THERMO header
            if line.strip().upper().startswith('THERMO'):
                i += 1
                continue

            # Check for species line (contains species name and thermo info)
            # Format: species_name  elements  G  T_low  T_high  T_common  ...
            # Line should have at least 79 characters
            if len(line) >= 79 and not line[0].isspace():
                try:
                    spec_name = line[:16].strip()

                    # Read four lines of thermo data
                    if i + 3 < len(lines):
                        line1 = lines[i]
                        line2 = lines[i + 1]
                        line3 = lines[i + 2]
                        line4 = lines[i + 3]

                        # Parse NASA polynomial data
                        thermo_data = self._parse_thermo_entry(line1, line2, line3, line4)
                        if thermo_data:
                            self.species_thermo[spec_name] = thermo_data
                            if spec_name not in self.species:
                                self.species[spec_name] = {'molweight': 0, 'elements': {}}

                        i += 4
                        continue
                except Exception as e:
                    pass

            i += 1

    def _parse_thermo_entry(self, line1, line2, line3, line4):
        """Parse a 4-line NASA 7-coefficient thermo entry (CHEMKIN format)."""
        try:
            # Line 1: species name (0-18), elements (20-40), phase (45), T_low (45-55), T_high (55-65), T_common (65-75)
            spec_name = line1[:18].strip()
            t_low = float(line1[45:55])
            t_high = float(line1[55:65])
            t_common = float(line1[65:75])

            # Line 2: a1, a2, a3, a4, a5 for high T (5 fields of 15 characters each)
            a_h = [
                float(line2[0:15]),
                float(line2[15:30]),
                float(line2[30:45]),
                float(line2[45:60]),
                float(line2[60:75]),
            ]

            # Line 3: a6, a7 for high T + a1, a2, a3, a4, a5 for low T
            a_h.extend([
                float(line3[0:15]),
                float(line3[15:30]),
            ])
            a_l = [
                float(line3[30:45]),
                float(line3[45:60]),
                float(line3[60:75]),
                float(line3[75:90]),
            ]

            # Line 4: a5, a6, a7 for low T
            a_l.extend([
                float(line4[0:15]),
                float(line4[15:30]),
                float(line4[30:45]),
            ])

            return {
                't_low': t_low,
                't_high': t_high,
                't_common': t_common,
                'a_high': a_h[:7],  # 7 coefficients
                'a_low': a_l[:7],   # 7 coefficients
            }
        except Exception as e:
            return None

    def parse_kinetics(self, filename):
        """Parse CHEMKIN kinetics file."""
        with open(filename, 'r') as f:
            content = f.read()

        # Extract elements
        elements_match = re.search(r'ELEMENTS?\s+(.*?)\s+END', content, re.IGNORECASE | re.DOTALL)
        if elements_match:
            elem_str = elements_match.group(1)
            self.elements = [e.strip() for e in elem_str.split() if e.strip()]

        # Extract species
        species_match = re.search(r'SPECIES?\s+(.*?)\s+END', content, re.IGNORECASE | re.DOTALL)
        if species_match:
            spec_str = species_match.group(1)
            spec_list = [s.strip() for s in re.split(r'\s+', spec_str) if s.strip()]
            for spec in spec_list:
                if spec not in self.species:
                    self.species[spec] = {'molweight': 0, 'elements': {}}

        # Extract reactions
        reactions_match = re.search(r'REACTIONS?\s+(.*?)\s+END', content, re.IGNORECASE | re.DOTALL)
        if reactions_match:
            rxn_str = reactions_match.group(1)
            # Simple reaction parser
            for line in rxn_str.split('\n'):
                line = line.strip()
                if not line or line.startswith('!'):
                    continue

                # Extract reaction equation and Arrhenius parameters
                rxn_data = self._parse_reaction_line(line)
                if rxn_data:
                    self.reactions.append(rxn_data)

    def _parse_reaction_line(self, line):
        """Parse a single reaction line."""
        # Format: reactants => products A beta Ea
        try:
            # Split by whitespace and arrow
            if '=>' in line:
                equation, params_str = line.split('=>')
                reactants = equation.strip()

                # Try to parse Arrhenius parameters from end of line
                parts = params_str.split()
                if len(parts) >= 3:
                    products = ' '.join(parts[:-3]).strip()
                    A = float(parts[-3])
                    beta = float(parts[-2])
                    Ta = float(parts[-1])

                    equation = (
                        f"{_format_reaction_side(reactants)} => "
                        f"{_format_reaction_side(products)}"
                    )

                    return {
                        'equation': equation,
                        'reactants': reactants,
                        'products': products,
                        'A': A,
                        'beta': beta,
                        'Ta': Ta,
                    }
        except:
            pass
        return None


def write_foam_thermo_simple(species_list, parser, output_file):
    """Write OpenFOAM thermo file from parsed CHEMKIN data."""
    lines = [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "| =========                 |                                                 |",
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |",
        "|  \\\\    /   O peration     | Version:  v2412                                 |",
        "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |",
        "|    \\/     M anipulation  |                                                 |",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        "    object      thermo.compressibleGas;",
        "}",
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
        "",
    ]

    # Standard NASA polynomial coefficients for common species if not found in file
    default_species_data = _get_default_species_data()

    for spec in sorted(species_list):
        lines.append(spec)
        lines.append("{")
        lines.append("    specie")
        lines.append("    {")

        # Try to get molecular weight
        molweight = default_species_data.get(spec, {}).get('molweight', 0)
        if molweight > 0:
            lines.append(f"        molWeight       {molweight:.4f};")
        else:
            # Estimate from species name (simple heuristic)
            est_weight = _estimate_molweight(spec)
            lines.append(f"        molWeight       {est_weight:.4f};")

        lines.append("    }")

        # Elements
        elements = _get_species_elements(spec)
        if elements:
            lines.append("    elements")
            lines.append("    {")
            for elem, count in sorted(elements.items()):
                if count > 0:
                    lines.append(f"        {elem}       {int(count)};")
            lines.append("    }")

        # Thermodynamics (use defaults for now)
        if spec in parser.species_thermo:
            data = parser.species_thermo[spec]
            lines.append("    thermodynamics")
            lines.append("    {")
            lines.append(f"        Tlow            {int(data['t_low'])};")
            lines.append(f"        Thigh           {int(data['t_high'])};")
            lines.append(f"        Tcommon         {int(data['t_common'])};")

            high_str = "        highCpCoeffs    ( "
            high_str += " ".join(f"{c:.5e}" for c in data['a_high'][:7])
            high_str += " );"
            lines.append(high_str)

            low_str = "        lowCpCoeffs     ( "
            low_str += " ".join(f"{c:.5e}" for c in data['a_low'][:7])
            low_str += " );"
            lines.append(low_str)

            lines.append("    }")
        else:
            # Use default data
            if spec in default_species_data:
                data = default_species_data[spec]
                lines.append("    thermodynamics")
                lines.append("    {")
                lines.append(f"        Tlow            {data.get('Tlow', 200)};")
                lines.append(f"        Thigh           {data.get('Thigh', 5000)};")
                lines.append(f"        Tcommon         {data.get('Tcommon', 1000)};")
                lines.append(f"        highCpCoeffs    ( {' '.join(str(c) for c in data.get('highCpCoeffs', [0]*7))} );")
                lines.append(f"        lowCpCoeffs     ( {' '.join(str(c) for c in data.get('lowCpCoeffs', [0]*7))} );")
                lines.append("    }")

        # Transport
        lines.append("    transport")
        lines.append("    {")
        lines.append("        As              1.67212e-06;")
        lines.append("        Ts              170.672;")
        lines.append("    }")
        lines.append("}")
        lines.append("")

    with open(output_file, 'w') as f:
        f.write("\n".join(lines))

    print(f"✓ Wrote {output_file} ({len(species_list)} species)")


def write_foam_reactions(species_list, parser, output_file):
    """Write OpenFOAM reactions file."""
    lines = [
        "elements",
        "(",
    ]

    # Write elements
    for elem in sorted(set(parser.elements)):
        if elem.strip():
            lines.append(elem)

    lines.extend([
        ");",
        "",
        "species",
        "(",
    ])

    # Write species
    for spec in sorted(species_list):
        lines.append("    " + spec)

    lines.extend([
        ");",
        "",
        "reactions",
        "{",
    ])

    # Write reactions (simplified)
    for i, rxn in enumerate(parser.reactions[:1000]):  # Limit to 1000 for now
        rxn_name = f"reaction{i:04d}"
        lines.append(f"    {rxn_name}")
        lines.append("    {")
        lines.append("        type     irreversibleArrheniusReaction;")

        # Try to format reaction equation
        if 'equation' in rxn and rxn['equation']:
            lines.append(f"        reaction \"{rxn['equation']}\";")
        else:
            lines.append(f"        reaction \"R{i}\";")

        lines.append(f"        A        {rxn['A']:.4e};")
        lines.append(f"        beta     {rxn['beta']:.1f};")
        lines.append(f"        Ta       {rxn['Ta']:.1f};")
        lines.append("    }")

    lines.extend([
        "}",
        "",
    ])

    with open(output_file, 'w') as f:
        f.write("\n".join(lines))

    print(f"✓ Wrote {output_file} ({len(parser.reactions)} reactions)")


def _estimate_molweight(species):
    """Estimate molecular weight from species name."""
    atomic_weights = {
        'H': 1.008, 'C': 12.011, 'O': 15.999, 'N': 14.007,
        'S': 32.06, 'AR': 39.948, 'HE': 4.003
    }

    weight = 0
    i = 0
    while i < len(species):
        # Extract element
        elem = species[i].upper()
        i += 1

        # Extract count
        count_str = ''
        while i < len(species) and species[i].isdigit():
            count_str += species[i]
            i += 1

        count = int(count_str) if count_str else 1

        if elem in atomic_weights:
            weight += atomic_weights[elem] * count
        elif elem and len(elem) > 1:
            # Try two-letter element
            elem2 = elem[:2]
            if elem2 in atomic_weights:
                weight += atomic_weights[elem2] * count
                i -= len(elem) - 2

    return max(weight, 2.0)  # At least 2


def _get_species_elements(species):
    """Extract element composition from species formula."""
    elements = {}
    atomic_weights = {
        'H': 1, 'C': 12, 'O': 16, 'N': 14, 'S': 32, 'AR': 40, 'HE': 4
    }

    i = 0
    while i < len(species):
        if species[i].isalpha():
            elem = species[i].upper()
            i += 1

            # Check for two-letter element
            if i < len(species) and species[i].isalpha() and species[i].islower():
                elem += species[i]
                i += 1

            # Extract count
            count_str = ''
            while i < len(species) and species[i].isdigit():
                count_str += species[i]
                i += 1

            count = int(count_str) if count_str else 1
            if elem in atomic_weights:
                elements[elem] = elements.get(elem, 0) + count
        else:
            i += 1

    return elements


def _get_default_species_data():
    """Return default NASA polynomial data for common species."""
    return {
        'N2': {
            'molweight': 28.014,
            'Tlow': 200,
            'Thigh': 6000,
            'Tcommon': 1000,
            'highCpCoeffs': [2.95258, 0.001396674, -4.92632e-07, 7.860863e-11, -4.607552e-15, -923.949, 5.87189],
            'lowCpCoeffs': [3.298677, 0.001408240, -3.963222e-06, 5.641515e-09, -2.444855e-12, -1020.900, 3.950372],
        },
        'O2': {
            'molweight': 31.999,
            'Tlow': 200,
            'Thigh': 5000,
            'Tcommon': 1000,
            'highCpCoeffs': [3.69758, 0.00061352, -1.25884e-07, 1.77528e-11, -1.13644e-15, -1233.93, 3.18917],
            'lowCpCoeffs': [3.21294, 0.00112749, -5.75615e-07, 1.31388e-09, -8.76855e-13, -1005.25, 6.03474],
        },
        'H2O': {
            'molweight': 18.015,
            'Tlow': 200,
            'Thigh': 3500,
            'Tcommon': 1000,
            'highCpCoeffs': [3.03399, 0.002176919, -1.640725e-07, -9.7041987e-11, 1.68200992e-14, -30004.2971, 4.9667701],
            'lowCpCoeffs': [4.19864519, -0.0020364341, 6.52040211e-06, -5.48797996e-09, 1.77197367e-12, -30293.7267, -0.849032283],
        },
        'CO2': {
            'molweight': 44.009,
            'Tlow': 200,
            'Thigh': 3500,
            'Tcommon': 1000,
            'highCpCoeffs': [3.85746, 0.00441437, -2.21481e-06, 5.23490e-10, -4.72084e-14, -48759.166, 2.27163],
            'lowCpCoeffs': [2.35677352, 0.00898459677, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -48371.9712, 9.90105222],
        },
    }


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    thermo_file = sys.argv[1]
    kinetics_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "."

    try:
        print("Parsing CHEMKIN files...")
        parser = ChemkinParser()

        parser.parse_thermo(thermo_file)
        print(f"  Thermo species: {len(parser.species_thermo)}")

        parser.parse_kinetics(kinetics_file)
        print(f"  Kinetics species: {len(parser.species)}")
        print(f"  Elements: {len(parser.elements)}")
        print(f"  Reactions: {len(parser.reactions)}")

        # Get all species
        all_species = sorted(set(list(parser.species.keys()) + list(parser.species_thermo.keys())))
        print(f"  Total unique species: {len(all_species)}")

        # Ensure output directory exists
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Write OpenFOAM files
        thermo_path = output_path / "thermo.compressibleGas"
        reactions_path = output_path / "reactions"

        write_foam_thermo_simple(all_species, parser, str(thermo_path))
        write_foam_reactions(all_species, parser, str(reactions_path))

        print(f"\nConversion complete!")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
