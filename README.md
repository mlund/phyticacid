# Phytic Acid Molecular Dynamics Simulation

This directory contains a complete setup for running a **Gromacs MD simulation** of **phytic acid** (myo-inositol hexakisphosphate, IP6) in water using the **CHARMM36 force field**.

## System Details

- **Molecule**: Phytic acid (IP6)
- **Force Field**: CHARMM36
- **Water Model**: TIP3P
- **Protonation State**: Physiological pH (~7.4), approximately 10 negative charges
- **Box Size**: 4 nm cubic box
- **Ionic Strength**: 0.15 M NaCl
- **Temperature**: 300 K
- **Pressure**: 1 bar
- **Production Run**: 10 ns

### Molecular Weight Note

The Wikipedia molecular weight of **660.029 g/mol** refers to the **fully protonated** (neutral) form: **C₆H₁₈O₂₄P₆**

At **physiological pH (~7.4)**, phytic acid is **mostly deprotonated** with the formula approximately **C₆H₁₀O₂₄P₆¹⁰⁻** (10 negative charges). The pKa values of the phosphate groups range from ~1.5 to ~10, so at pH 7.4:
- Some phosphate groups are fully deprotonated (PO₄²⁻)
- Some have one proton (HPO₄⁻)
- Total charge: approximately -10

The provided PDB structure represents this physiologically relevant, deprotonated form.

## File Descriptions

```
phytic_acid.pdb       - Initial structure (deprotonated at pH 7.4)
em.mdp                - Energy minimization parameters
nvt.mdp               - NVT equilibration parameters (100 ps)
npt.mdp               - NPT equilibration parameters (100 ps)
md.mdp                - Production MD parameters (10 ns)
posre.itp             - Position restraints for equilibration
setup_simulation.sh   - Automated setup script
analysis.sh           - Post-simulation analysis script
```

## Prerequisites

1. **Gromacs** (version 2020 or later recommended)
   ```bash
   gmx --version
   ```

2. **CHARMM36 Force Field** topology files for phytic acid

## Step 1: Generate Topology Files (CHARMM-GUI)

Since phytic acid is not a standard molecule, you need to generate topology files using **CHARMM-GUI**:

### Option A: CHARMM-GUI (Recommended)

1. Go to https://www.charmm-gui.org/
2. Select **Input Generator** → **Ligand Reader & Modeler**
3. Upload `phytic_acid.pdb`
4. Select force field: **CHARMM36**
5. Set the net charge to **-10** (physiological pH)
6. Download the generated files
7. Extract these files to your simulation directory:
   - `phytic_acid.itp` - Molecule topology
   - `phytic_acid.prm` - Force field parameters (if needed)
   - `topol.top` - Main topology file

### Option B: CGenFF (Advanced)

If you have CGenFF installed:

```bash
# Convert PDB to mol2
obabel phytic_acid.pdb -O phytic_acid.mol2

# Generate topology using CGenFF
python cgenff_charmm2gmx.py PHYTIC phytic_acid.mol2 phytic_acid.str charmm36.ff

# This creates phytic_acid.itp and phytic_acid.prm
```

### Manual Topology Template

Create `topol.top`:

```
; Include force field
#include "charmm36.ff/forcefield.itp"

; Include phytic acid topology
#include "phytic_acid.itp"

; Include water topology
#include "charmm36.ff/tip3p.itp"

#ifdef POSRES_WATER
#include "charmm36.ff/tip3p_posre.itp"
#endif

; Include position restraints
#ifdef POSRES
#include "posre.itp"
#endif

; Include ion topology
#include "charmm36.ff/ions.itp"

[ system ]
; Name
Phytic acid in water

[ molecules ]
; Compound        #mols
IP6                 1
; SOL and ions will be added by setup script
```

## Step 2: Run the Setup Script

Once you have the topology files, run the automated setup:

```bash
chmod +x setup_simulation.sh
./setup_simulation.sh
```

This script will:
1. Create a 4 nm cubic simulation box
2. Solvate with TIP3P water (~2000 water molecules)
3. Add Na⁺ and Cl⁻ ions (neutralize + 0.15 M)
4. Run energy minimization (~50,000 steps)
5. Run NVT equilibration (100 ps at 300 K)
6. Run NPT equilibration (100 ps at 1 bar)
7. Prepare production MD (10 ns)

## Step 3: Run Production MD

After setup completes successfully:

```bash
# Standard CPU run
gmx mdrun -v -deffnm md

# With GPU acceleration
gmx mdrun -v -deffnm md -nb gpu

# Parallel execution (4 threads)
gmx mdrun -v -deffnm md -nt 4

# With GPU and multiple CPUs
gmx mdrun -v -deffnm md -nb gpu -nt 8 -ntmpi 1 -ntomp 8
```

**Estimated runtime**:
- CPU only: ~2-4 hours
- With GPU: ~30-60 minutes

## Step 4: Analysis

After the simulation completes, run the analysis script:

```bash
chmod +x analysis.sh
./analysis.sh
```

This generates:
- **RMSD**: Structural stability over time
- **Radius of gyration**: Molecular compactness
- **Hydrogen bonds**: H-bond network analysis
- **RDF**: Phytic acid-water distribution
- **SASA**: Solvent accessible surface area
- **Energy**: Potential, kinetic, total energy
- **Phosphate distances**: Inter-phosphate group distances

Results are saved in the `analysis/` directory as `.xvg` files.

## Visualization

### VMD (Recommended)
```bash
vmd em.gro md.xtc
```

### PyMOL
```bash
pymol em.gro
# Load trajectory: File → Load Trajectory
```

### Gromacs native
```bash
gmx view -f md.xtc -s md.tpr
```

## Troubleshooting

### "Fatal error: No such moleculetype IP6"
- Make sure `phytic_acid.itp` is in the same directory
- Check that `topol.top` correctly includes the `.itp` file

### Energy minimization not converging
- Increase `emtol` in `em.mdp` (e.g., to 5000)
- Check for clashes in initial structure
- Verify topology charge matches actual charge

### "LINCS WARNING"
- Reduce time step in `.mdp` files (dt = 0.001)
- Check temperature coupling parameters
- Ensure proper equilibration before production

### High pressure fluctuations
- Normal during NPT equilibration
- Should stabilize around 1 bar (±200 bar fluctuations are normal)

## Expected Results

At physiological pH, phytic acid should:
- Remain highly solvated (strong H-bonds with water)
- Show relatively compact structure (multiple intramolecular interactions)
- Coordinate strongly with Na⁺ ions
- Maintain stable phosphate geometry

## References

1. **Phytic acid structure**:
   - PubChem CID: 890
   - Wikipedia: https://en.wikipedia.org/wiki/Phytic_acid

2. **CHARMM36 Force Field**:
   - Best et al., J. Chem. Theory Comput. 2012, 8, 3257-3273

3. **Gromacs**:
   - Abraham et al., SoftwareX 2015, 1-2, 19-25
   - Manual: https://manual.gromacs.org/

4. **CHARMM-GUI**:
   - Jo et al., J. Comput. Chem. 2008, 29, 1859-1865
   - Website: https://www.charmm-gui.org/

## Extending the Simulation

To run longer simulations, edit `md.mdp`:

```bash
# For 50 ns (change nsteps)
nsteps = 25000000    ; 50 ns

# For 100 ns
nsteps = 50000000    ; 100 ns
```

Then regenerate the `.tpr` file:
```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_long.tpr
gmx mdrun -v -deffnm md_long
```

## Contact & Support

For Gromacs support: http://www.gromacs.org/Documentation
For CHARMM-GUI support: https://www.charmm-gui.org/?doc=contactus

---

**Generated for CHARMM36 force field with TIP3P water model**
