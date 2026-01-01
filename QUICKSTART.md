  # Quick Start Guide - Phytic Acid MD Simulation

## Prerequisites Check

```bash
# Check if Gromacs is installed
gmx --version

# Check if you have CHARMM36 force field
ls -d charmm36*.ff/
```

**If Gromacs is not installed:**
- macOS: `brew install gromacs`
- Linux: `sudo apt install gromacs` or compile from source
- Check: http://www.gromacs.org/Downloads

**If CHARMM36 is missing:**
Download from: https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
Or the simulation will fail with "File not found" errors.

## Option 1: Quick Test (No Gromacs Required)

Just want to verify the structure?

```bash
# Verify phytic acid structure and charge
python3 verify_structure.py

# Expected output:
# - Molecular Formula: C6H8O24P6
# - Molecular Weight: 649.950 g/mol
# - Estimated net charge: -10
```

## Option 2: Full Simulation Setup

### Step 1: Verify CHARMM36 Force Field Path

The topology file expects the force field in `charmm36-jul2022.ff/`. If your force field directory has a different name, update `topol.top`:

```bash
# Find your CHARMM force field
ls -d charmm36*.ff/

# Edit topol.top to match your force field directory name
# Change: #include "charmm36-jul2022.ff/forcefield.itp"
# To:     #include "YOUR_DIRECTORY_NAME/forcefield.itp"
```

### Step 2: Run Setup

```bash
# Make scripts executable
chmod +x setup_simulation.sh analysis.sh

# Run the full setup (5-10 minutes)
./setup_simulation.sh
```

This will:
1. ✓ Create 4 nm simulation box
2. ✓ Add ~2000 water molecules
3. ✓ Add Na⁺/Cl⁻ ions (neutralize + 0.15 M)
4. ✓ Energy minimization (~2-5 min)
5. ✓ NVT equilibration (~2 min)
6. ✓ NPT equilibration (~2 min)
7. ✓ Prepare 10 ns production run

**Common Issues:**

❌ **"No such file or directory: charmm36-jul2022.ff"**
- Download CHARMM36 from https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
- Extract to simulation directory
- Update `topol.top` with correct path

❌ **"Fatal error: No such moleculetype IP6"**
- Make sure `phytic_acid.itp` is in the same directory
- Check that `topol.top` includes it correctly

### Step 3: Run Production MD

After setup completes successfully:

```bash
# Standard run (CPU, ~2-4 hours for 10 ns)
gmx mdrun -v -deffnm md

# With GPU (much faster, ~30-60 min)
gmx mdrun -v -deffnm md -nb gpu

# Check progress (in another terminal)
tail -f md.log
```

**To run in background:**
```bash
nohup gmx mdrun -v -deffnm md > md_output.log 2>&1 &

# Check progress
tail -f md_output.log
```

### Step 4: Analysis

After simulation completes:

```bash
# Run all analyses
./analysis.sh

# Results will be in analysis/ directory
ls analysis/

# Visualize results (requires matplotlib)
pip install matplotlib numpy
python3 visualize.py summary
```

## Option 3: Step-by-Step Manual Setup

If the automated script fails, run manually:

```bash
# 1. Create box
gmx editconf -f phytic_acid.pdb -o boxed.gro -c -d 1.5 -bt cubic

# 2. Solvate
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# 3. Add ions
gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# 4. Energy minimization
gmx grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# 5. NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# 6. NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# 7. Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```

## Visualization

### View trajectory with VMD
```bash
vmd em.gro md.xtc
```

### View with PyMOL
```bash
pymol em.gro
# Then: File → Load Trajectory → select md.xtc
```

### View with Gromacs
```bash
gmx view -f md.xtc -s md.tpr
```

## Expected Timeline

| Step | CPU Time | GPU Time |
|------|----------|----------|
| Setup (EM + equilibration) | 5-10 min | 3-5 min |
| Production MD (10 ns) | 2-4 hours | 30-60 min |
| Analysis | 2-5 min | 2-5 min |
| **Total** | **2-4.5 hours** | **35-70 min** |

## Key Files Reference

| File | Description |
|------|-------------|
| `phytic_acid.pdb` | Initial structure (44 atoms, -10 charge) |
| `phytic_acid.itp` | Molecule topology for CHARMM36 |
| `topol.top` | Main topology file |
| `em.mdp` | Energy minimization parameters |
| `nvt.mdp` | NVT equilibration (300 K) |
| `npt.mdp` | NPT equilibration (1 bar) |
| `md.mdp` | Production MD (10 ns) |
| `posre.itp` | Position restraints |
| `setup_simulation.sh` | Automated setup script |
| `analysis.sh` | Post-simulation analysis |
| `verify_structure.py` | Structure verification |
| `visualize.py` | Plot analysis results |

## Troubleshooting

### Simulation crashes with LINCS warnings
```bash
# Reduce time step in md.mdp
dt = 0.001   # instead of 0.002
```

### Energy minimization doesn't converge
```bash
# Increase tolerance in em.mdp
emtol = 5000.0   # instead of 1000.0
```

### "Fatal error: atomtype XXX not found"
- Your CHARMM36 version might be different
- Try using CHARMM-GUI to generate the topology instead
- See README.md for CHARMM-GUI instructions

## Next Steps

After completing the 10 ns simulation:

1. **Check stability**: Plot RMSD - should plateau
2. **Extend simulation**: Edit `md.mdp` nsteps for longer runs
3. **Analyze interactions**: Use `gmx hbond`, `gmx rdf`, `gmx mindist`
4. **Study ion coordination**: How many Na⁺ ions coordinate with phosphates?
5. **Calculate pKa**: Use constant-pH MD (advanced)

## Getting Help

- Gromacs manual: https://manual.gromacs.org/
- Gromacs forums: https://gromacs.bioexcel.eu/
- CHARMM-GUI: https://www.charmm-gui.org/
- Tutorial: http://www.mdtutorials.com/gmx/

## Citations

If you use this setup, please cite:

1. **Gromacs**: Abraham et al., SoftwareX 2015, 1-2, 19-25
2. **CHARMM36**: Best et al., J. Chem. Theory Comput. 2012, 8, 3257-3273
3. **Phytic acid**: Torres et al., J. Agric. Food Chem. 2005, 53, 2618-2623

---

*Generated for CHARMM36 force field - Phytic acid at pH 7.4*
