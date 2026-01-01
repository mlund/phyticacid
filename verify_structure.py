#!/usr/bin/env python3
"""
Verify phytic acid structure and calculate properties
"""

import sys

def parse_pdb(filename):
    """Parse PDB file and extract atomic information"""
    atoms = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'element': line[76:78].strip() or line[12:16].strip()[0],
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'name': line[12:16].strip()
                }
                atoms.append(atom)
    return atoms

def calculate_molecular_weight(atoms):
    """Calculate molecular weight from atom list"""
    atomic_masses = {
        'H': 1.008,
        'C': 12.011,
        'O': 15.999,
        'P': 30.974
    }

    total_mass = 0.0
    composition = {}

    for atom in atoms:
        element = atom['element']
        if element not in atomic_masses:
            # Try to infer from atom name
            element = atom['name'][0]

        if element in atomic_masses:
            total_mass += atomic_masses[element]
            composition[element] = composition.get(element, 0) + 1
        else:
            print(f"Warning: Unknown element {element}")

    return total_mass, composition

def calculate_charge(composition):
    """Estimate charge based on composition"""
    # For phytic acid at pH 7.4
    # C6H(18-x)O24P6 where x is number of deprotonated sites

    n_h = composition.get('H', 0)
    n_protons_lost = 18 - n_h  # Assuming fully protonated form has 18 H

    # Each phosphate can lose up to 2 protons
    # Typical charge at pH 7.4 is around -10

    return -n_protons_lost

def calculate_center_of_mass(atoms):
    """Calculate center of mass"""
    atomic_masses = {
        'H': 1.008,
        'C': 12.011,
        'O': 15.999,
        'P': 30.974
    }

    total_mass = 0.0
    com = [0.0, 0.0, 0.0]

    for atom in atoms:
        element = atom['element']
        if element not in atomic_masses:
            element = atom['name'][0]

        if element in atomic_masses:
            mass = atomic_masses[element]
            total_mass += mass
            com[0] += mass * atom['x']
            com[1] += mass * atom['y']
            com[2] += mass * atom['z']

    com = [c / total_mass for c in com]
    return com

def calculate_radius_of_gyration(atoms, com):
    """Calculate radius of gyration"""
    atomic_masses = {
        'H': 1.008,
        'C': 12.011,
        'O': 15.999,
        'P': 30.974
    }

    total_mass = 0.0
    rg_sq = 0.0

    for atom in atoms:
        element = atom['element']
        if element not in atomic_masses:
            element = atom['name'][0]

        if element in atomic_masses:
            mass = atomic_masses[element]
            total_mass += mass

            dx = atom['x'] - com[0]
            dy = atom['y'] - com[1]
            dz = atom['z'] - com[2]

            rg_sq += mass * (dx*dx + dy*dy + dz*dz)

    rg = (rg_sq / total_mass) ** 0.5
    return rg

def main():
    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
    else:
        pdb_file = 'phytic_acid.pdb'

    print("=" * 60)
    print("Phytic Acid Structure Verification")
    print("=" * 60)
    print()

    try:
        atoms = parse_pdb(pdb_file)
        print(f"✓ Successfully parsed {pdb_file}")
        print(f"  Total atoms: {len(atoms)}")
        print()

        # Calculate molecular weight
        mw, composition = calculate_molecular_weight(atoms)
        print("Molecular Composition:")
        for element in sorted(composition.keys()):
            print(f"  {element}: {composition[element]}")
        print()

        # Molecular formula
        formula = f"C{composition.get('C', 0)}H{composition.get('H', 0)}O{composition.get('O', 0)}P{composition.get('P', 0)}"
        print(f"Molecular Formula: {formula}")
        print(f"Molecular Weight: {mw:.3f} g/mol")
        print()

        # Reference values
        print("Reference Values:")
        print("  Fully protonated (C6H18O24P6): 660.029 g/mol")
        print(f"  This structure: {mw:.3f} g/mol")
        print(f"  Difference: {mw - 660.029:.3f} g/mol")
        print()

        # Estimate charge
        estimated_charge = calculate_charge(composition)
        print(f"Estimated net charge: {estimated_charge}")
        print("  (At pH 7.4, expected: ~-10)")
        print()

        # Calculate geometric properties
        com = calculate_center_of_mass(atoms)
        print(f"Center of mass: ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}) nm")

        rg = calculate_radius_of_gyration(atoms, com)
        print(f"Radius of gyration: {rg:.3f} nm")
        print()

        # Count specific atom types
        print("Atom type counts:")
        print(f"  Carbons (ring): {composition.get('C', 0)}")
        print(f"  Oxygens: {composition.get('O', 0)}")
        print(f"  Phosphorus: {composition.get('P', 0)}")
        print(f"  Hydrogens: {composition.get('H', 0)}")
        print()

        # Validation
        print("Validation:")
        checks = []

        if composition.get('C', 0) == 6:
            print("  ✓ Correct number of carbons (6)")
            checks.append(True)
        else:
            print(f"  ✗ Wrong number of carbons (expected 6, got {composition.get('C', 0)})")
            checks.append(False)

        if composition.get('P', 0) == 6:
            print("  ✓ Correct number of phosphorus (6)")
            checks.append(True)
        else:
            print(f"  ✗ Wrong number of phosphorus (expected 6, got {composition.get('P', 0)})")
            checks.append(False)

        if composition.get('O', 0) == 24:
            print("  ✓ Correct number of oxygens (24)")
            checks.append(True)
        else:
            print(f"  ✗ Wrong number of oxygens (expected 24, got {composition.get('O', 0)})")
            checks.append(False)

        n_protons = 18 - estimated_charge
        if 8 <= n_protons <= 12:
            print(f"  ✓ Reasonable protonation for pH 7.4")
            checks.append(True)
        else:
            print(f"  ⚠ Unusual protonation state")
            checks.append(False)

        print()
        if all(checks):
            print("✓ Structure validation PASSED")
        else:
            print("⚠ Structure validation failed some checks")

        print()
        print("=" * 60)

    except FileNotFoundError:
        print(f"Error: File '{pdb_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
