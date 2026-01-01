#!/bin/bash
# Check if all requirements are met for running the simulation

echo "======================================================"
echo "Checking Requirements for Phytic Acid MD Simulation"
echo "======================================================"
echo ""

ALL_OK=true

# Check for Gromacs
echo -n "Checking for Gromacs... "
if command -v gmx &> /dev/null; then
    VERSION=$(gmx --version 2>&1 | grep "GROMACS version" | head -n1)
    echo "✓ Found"
    echo "  $VERSION"
else
    echo "✗ NOT FOUND"
    echo "  Please install Gromacs: http://www.gromacs.org/Downloads"
    ALL_OK=false
fi
echo ""

# Check for CHARMM36 force field
echo -n "Checking for CHARMM36 force field... "
if ls charmm36*.ff/ &> /dev/null; then
    FFDIR=$(ls -d charmm36*.ff/ | head -n1)
    echo "✓ Found: $FFDIR"

    # Check if topol.top matches
    if grep -q "$FFDIR" topol.top; then
        echo "  ✓ topol.top correctly references this force field"
    else
        echo "  ⚠ WARNING: topol.top may need updating"
        echo "    Update the #include paths in topol.top to use: $FFDIR"
    fi
else
    echo "✗ NOT FOUND"
    echo "  Download from: https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs"
    echo "  Extract to this directory"
    ALL_OK=false
fi
echo ""

# Check for Python
echo -n "Checking for Python 3... "
if command -v python3 &> /dev/null; then
    PYVERSION=$(python3 --version)
    echo "✓ Found: $PYVERSION"
else
    echo "⚠ NOT FOUND"
    echo "  Python 3 is optional but recommended for analysis scripts"
fi
echo ""

# Check for Python packages (optional)
echo -n "Checking for matplotlib (optional)... "
if python3 -c "import matplotlib" &> /dev/null; then
    echo "✓ Found"
else
    echo "⚠ NOT FOUND"
    echo "  Install with: pip install matplotlib numpy"
    echo "  (Optional: needed for visualize.py)"
fi
echo ""

# Check for VMD (optional)
echo -n "Checking for VMD (optional)... "
if command -v vmd &> /dev/null; then
    echo "✓ Found"
else
    echo "⚠ NOT FOUND"
    echo "  VMD is optional but recommended for visualization"
    echo "  Download from: https://www.ks.uiuc.edu/Research/vmd/"
fi
echo ""

# Check required files
echo "Checking required files:"
REQUIRED_FILES=(
    "phytic_acid.pdb"
    "phytic_acid.itp"
    "topol.top"
    "em.mdp"
    "nvt.mdp"
    "npt.mdp"
    "md.mdp"
    "posre.itp"
    "setup_simulation.sh"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file - MISSING"
        ALL_OK=false
    fi
done
echo ""

# Summary
echo "======================================================"
if [ "$ALL_OK" = true ]; then
    echo "✓ All REQUIRED components are installed!"
    echo ""
    echo "You can proceed with the simulation:"
    echo "  ./setup_simulation.sh"
else
    echo "✗ Some REQUIRED components are missing."
    echo ""
    echo "Please install missing requirements before proceeding."
    echo "See QUICKSTART.md for installation instructions."
fi
echo "======================================================"
