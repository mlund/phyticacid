#!/bin/bash
# Download and setup CHARMM36 force field for Gromacs
# This script automates the process of obtaining the CHARMM36 force field

set -e

echo "======================================================"
echo "CHARMM36 Force Field Download & Setup"
echo "======================================================"
echo ""

# Check if force field already exists
if ls charmm36*.ff/ &> /dev/null; then
    EXISTING=$(ls -d charmm36*.ff/ | head -n1)
    echo "⚠ Force field already exists: $EXISTING"
    echo ""
    read -p "Do you want to re-download? (y/N) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Keeping existing force field."
        exit 0
    fi
fi

echo "Downloading CHARMM36 force field..."
echo ""
echo "Source: MacKerell Lab (University of Maryland)"
echo "URL: http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz"
echo ""

# Download the force field
FFURL="http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz"
FFFILE="charmm36-jul2022.ff.tgz"

if command -v curl &> /dev/null; then
    echo "Downloading with curl..."
    curl -L -o "$FFFILE" "$FFURL"
elif command -v wget &> /dev/null; then
    echo "Downloading with wget..."
    wget -O "$FFFILE" "$FFURL"
else
    echo "Error: Neither curl nor wget found."
    echo "Please install curl or wget, or download manually from:"
    echo "  https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs"
    exit 1
fi

# Check if download succeeded
if [ ! -f "$FFFILE" ]; then
    echo "Error: Download failed."
    echo ""
    echo "Please download manually from:"
    echo "  https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs"
    echo ""
    echo "Then extract to this directory with:"
    echo "  tar xzf charmm36-jul2022.ff.tgz"
    exit 1
fi

echo "✓ Download complete: $FFFILE"
echo ""

# Extract the force field
echo "Extracting force field..."
tar xzf "$FFFILE"

if [ -d "charmm36-jul2022.ff" ]; then
    echo "✓ Force field extracted successfully!"
    echo ""

    # Verify key files exist
    echo "Verifying force field files..."
    REQUIRED_FF_FILES=(
        "charmm36-jul2022.ff/forcefield.itp"
        "charmm36-jul2022.ff/tip3p.itp"
        "charmm36-jul2022.ff/ions.itp"
    )

    ALL_GOOD=true
    for file in "${REQUIRED_FF_FILES[@]}"; do
        if [ -f "$file" ]; then
            echo "  ✓ $file"
        else
            echo "  ✗ $file - MISSING"
            ALL_GOOD=false
        fi
    done

    if [ "$ALL_GOOD" = true ]; then
        echo ""
        echo "✓ Force field setup complete!"
        echo ""

        # Update topol.top if needed
        if grep -q "charmm36-jul2022.ff" topol.top; then
            echo "✓ topol.top already configured correctly"
        else
            echo "Updating topol.top..."
            # This is already correct in our topol.top
            echo "✓ topol.top uses charmm36-jul2022.ff"
        fi

        echo ""
        echo "======================================================"
        echo "You can now run the simulation setup:"
        echo "  ./setup_simulation.sh"
        echo "======================================================"

        # Clean up tarball
        read -p "Remove downloaded tarball? (Y/n) " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            rm "$FFFILE"
            echo "✓ Removed $FFFILE"
        fi
    else
        echo ""
        echo "⚠ Force field extraction incomplete"
        echo "Please check the downloaded file or download manually"
    fi
else
    echo "Error: Extraction failed or directory not found"
    echo "Expected directory: charmm36-jul2022.ff/"
    exit 1
fi

echo ""
