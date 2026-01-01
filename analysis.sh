#!/bin/bash
# Analysis script for phytic acid MD simulation

echo "======================================================"
echo "Phytic Acid MD Simulation Analysis"
echo "======================================================"
echo ""

if [ ! -f "md.tpr" ] || [ ! -f "md.xtc" ]; then
    echo "ERROR: Production MD files not found!"
    echo "Please run the simulation first."
    exit 1
fi

# Create analysis directory
mkdir -p analysis
cd analysis

# 1. RMSD analysis
echo "1. Calculating RMSD..."
echo "Backbone" | gmx rms -s ../md.tpr -f ../md.xtc -o rmsd.xvg -tu ns

# 2. Radius of gyration
echo "2. Calculating radius of gyration..."
echo "System" | gmx gyrate -s ../md.tpr -f ../md.xtc -o gyrate.xvg

# 3. Hydrogen bonds
echo "3. Analyzing hydrogen bonds..."
echo "System System" | gmx hbond -s ../md.tpr -f ../md.xtc -num hbnum.xvg

# 4. Radial distribution function (RDF) - phytic acid to water
echo "4. Calculating RDF (phytic acid - water)..."
echo "1 2" | gmx rdf -s ../md.tpr -f ../md.xtc -o rdf_phytic_water.xvg -ref "resname IP6" -sel "resname SOL"

# 5. SASA (Solvent Accessible Surface Area)
echo "5. Calculating SASA..."
echo "System" | gmx sasa -s ../md.tpr -f ../md.xtc -o sasa.xvg -tu ns

# 6. Density profile
echo "6. Calculating density profile..."
gmx density -s ../md.tpr -f ../md.xtc -o density_profile.xvg -d Z

# 7. Energy analysis
echo "7. Extracting energy components..."
echo "Potential Kinetic-En. Total-Energy Temperature Pressure Density" | gmx energy -f ../md.edr -o energy.xvg

# 8. Distance analysis (between phosphate groups)
echo "8. Calculating inter-phosphate distances..."
# Example: distance between P1 and P4
echo "a P1\na P4\n" | gmx distance -s ../md.tpr -f ../md.xtc -select -oall phosphate_distance.xvg

echo ""
echo "======================================================"
echo "Analysis complete! Results saved in analysis/ folder"
echo "======================================================"
echo ""
echo "Generated files:"
echo "  - rmsd.xvg: RMSD over time"
echo "  - gyrate.xvg: Radius of gyration"
echo "  - hbnum.xvg: Number of hydrogen bonds"
echo "  - rdf_phytic_water.xvg: Radial distribution function"
echo "  - sasa.xvg: Solvent accessible surface area"
echo "  - density_profile.xvg: Density along z-axis"
echo "  - energy.xvg: Energy components"
echo "  - phosphate_distance.xvg: Inter-phosphate distances"
echo ""
echo "To visualize XVG files, use xmgrace or Python matplotlib"
echo ""
