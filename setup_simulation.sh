#!/bin/bash
# Gromacs simulation setup script for phytic acid in water
# CHARMM36 force field

set -e  # Exit on error

echo "======================================================"
echo "Phytic Acid MD Simulation Setup"
echo "Force Field: CHARMM36"
echo "======================================================"
echo ""

# Check if topology files exist
if [ ! -f "phytic_acid.itp" ] || [ ! -f "topol.top" ]; then
    echo "ERROR: Topology files not found!"
    echo "Please generate topology files using CHARMM-GUI first."
    echo "See README.md for instructions."
    exit 1
fi

# Step 1: Create simulation box (4 nm cubic)
echo "Step 1: Creating simulation box (4 nm cubic)..."
gmx editconf -f phytic_acid.pdb -o boxed.gro -c -d 1.5 -bt cubic

# Step 2: Solvate with TIP3P water
echo "Step 2: Solvating with TIP3P water..."
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# Step 3: Add ions for charge neutralization
echo "Step 3: Preparing for ion addition..."
gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

echo "Adding ions (neutralize charge + 0.15 M NaCl)..."
echo "SOL" | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# Step 4: Energy Minimization
echo "Step 4: Energy minimization..."
gmx grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

echo "Checking energy minimization convergence..."
echo "Potential" | gmx energy -f em.edr -o potential.xvg

# Step 5: NVT Equilibration
echo "Step 5: NVT equilibration (100 ps at 300 K)..."
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

echo "Checking temperature equilibration..."
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg

# Step 6: NPT Equilibration
echo "Step 6: NPT equilibration (100 ps at 1 bar)..."
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

echo "Checking pressure and density equilibration..."
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
echo "Density" | gmx energy -f npt.edr -o density.xvg

# Step 7: Production MD
echo "Step 7: Preparing production MD (10 ns)..."
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

echo ""
echo "======================================================"
echo "Setup complete!"
echo "======================================================"
echo ""
echo "To run the production simulation:"
echo "  gmx mdrun -v -deffnm md"
echo ""
echo "Or with GPU acceleration:"
echo "  gmx mdrun -v -deffnm md -nb gpu"
echo ""
echo "For parallel execution (e.g., 4 threads):"
echo "  gmx mdrun -v -deffnm md -nt 4"
echo ""
