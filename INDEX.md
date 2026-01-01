# Phytic Acid MD Simulation - File Index

## Quick Navigation

### 🚀 Getting Started
1. **First time?** → Read [QUICKSTART.md](QUICKSTART.md)
2. **Check requirements** → Run `./check_requirements.sh`
3. **Download force field** → Run `./download_charmm36.sh`
4. **Run simulation** → Run `./setup_simulation.sh`

### 📚 Documentation Files

| File | Purpose | Read When |
|------|---------|-----------|
| [QUICKSTART.md](QUICKSTART.md) | Quick start guide | Starting the project |
| [README.md](README.md) | Complete documentation | Need detailed info |
| [SIMULATION_NOTES.md](SIMULATION_NOTES.md) | Technical details | Understanding the science |
| [INDEX.md](INDEX.md) | This file | Finding specific info |

---

## File Categories

### 🧬 Core Structure Files

**phytic_acid.pdb** (3.6 KB)
- Initial 3D structure of phytic acid
- 44 atoms total (6C, 24O, 6P, 8H)
- Protonation state: pH 7.4 (-10 charge)
- Use: Starting structure for simulation

**phytic_acid.itp** (6.3 KB)
- CHARMM36 topology for phytic acid
- Contains: atoms, bonds, angles, dihedrals
- Defines force field parameters
- Use: Referenced by topol.top

**topol.top** (721 B)
- Main topology file
- Includes force field and molecule definitions
- Lists system composition
- Use: Required for all Gromacs commands

**posre.itp** (1.6 KB)
- Position restraints for equilibration
- Applied during NVT and NPT phases
- Restrains heavy atoms with 1000 kJ/mol/nm²
- Use: Automatically included when -DPOSRES is defined

---

### ⚙️ Simulation Parameter Files

**em.mdp** (1.2 KB)
- Energy minimization parameters
- Steepest descent, 50,000 max steps
- Tolerance: 1000 kJ/mol/nm
- Use: First step to remove bad contacts

**nvt.mdp** (2.0 KB)
- NVT equilibration (constant volume)
- Duration: 100 ps at 300 K
- V-rescale thermostat
- Use: Temperature equilibration with restraints

**npt.mdp** (1.6 KB)
- NPT equilibration (constant pressure)
- Duration: 100 ps at 1 bar, 300 K
- Parrinello-Rahman barostat
- Use: Pressure/density equilibration

**md.mdp** (1.4 KB)
- Production MD parameters
- Duration: 10 ns (5,000,000 steps)
- No restraints, full dynamics
- Use: Main production simulation

---

### 🤖 Automation Scripts

**setup_simulation.sh** (2.7 KB)
- Automated setup pipeline
- Runs: box creation → solvation → ions → EM → NVT → NPT → prep MD
- Checks for topology files first
- Use: `./setup_simulation.sh`

**analysis.sh** (2.4 KB)
- Post-simulation analysis
- Generates: RMSD, Rg, H-bonds, RDF, SASA, energy
- Creates analysis/ directory
- Use: `./analysis.sh` (after MD completes)

**check_requirements.sh** (3.1 KB)
- Verify all dependencies installed
- Checks: Gromacs, Python, force field, files
- Shows what's missing
- Use: `./check_requirements.sh`

**download_charmm36.sh** (NEW)
- Automated CHARMM36 download
- Downloads and extracts force field
- Verifies installation
- Use: `./download_charmm36.sh`

---

### 🔬 Analysis Tools

**verify_structure.py** (6.6 KB)
- Validate phytic acid structure
- Calculates: MW, composition, charge, Rg
- Checks atom counts
- Use: `python3 verify_structure.py phytic_acid.pdb`

**visualize.py** (6.2 KB)
- Plot XVG files from Gromacs
- Create summary figures
- Statistical analysis
- Use:
  - `python3 visualize.py analysis/rmsd.xvg`
  - `python3 visualize.py summary`
  - `python3 visualize.py all`

---

### 📖 Documentation

**README.md** (6.8 KB)
- Complete project documentation
- Includes: setup, usage, troubleshooting
- Force field generation instructions
- References and citations

**QUICKSTART.md** (6.0 KB)
- Fast-track guide for experienced users
- Three setup options (test, full, manual)
- Timeline estimates
- Common issues and solutions

**SIMULATION_NOTES.md** (7.2 KB)
- Detailed scientific notes
- System specifications and expectations
- Analysis checklist
- Extension suggestions
- Known limitations

**.gitignore** (721 B)
- Git ignore rules
- Excludes: simulation outputs, backups, force field
- Use: Automatic when using git

---

## Workflow Diagram

```
1. START
   ↓
2. check_requirements.sh ← Check if ready
   ↓
3. download_charmm36.sh  ← Get force field
   ↓
4. verify_structure.py   ← (Optional) Verify structure
   ↓
5. setup_simulation.sh   ← Run full setup
   ↓                        (creates: em, nvt, npt, md.tpr)
6. gmx mdrun -v -deffnm md ← Run production MD
   ↓                        (creates: md.xtc, md.edr, md.log)
7. analysis.sh           ← Analyze results
   ↓                        (creates: analysis/*.xvg)
8. visualize.py summary  ← Create plots
   ↓
9. DONE - Review results in analysis/
```

---

## Quick Reference Commands

### Setup
```bash
./check_requirements.sh           # Check dependencies
./download_charmm36.sh            # Get CHARMM36 force field
python3 verify_structure.py       # Verify structure
./setup_simulation.sh             # Full automated setup
```

### Run Simulation
```bash
gmx mdrun -v -deffnm md           # CPU only
gmx mdrun -v -deffnm md -nb gpu   # With GPU
```

### Analysis
```bash
./analysis.sh                     # Run all analyses
python3 visualize.py summary      # Create summary plot
python3 visualize.py all          # Plot all XVG files
```

### Visualization
```bash
vmd em.gro md.xtc                 # VMD
pymol em.gro                      # PyMOL
gmx view -f md.xtc -s md.tpr      # Gromacs viewer
```

---

## File Size Summary

| Category | Files | Total Size |
|----------|-------|------------|
| Structure/Topology | 4 | ~12 KB |
| Parameters (.mdp) | 4 | ~6 KB |
| Scripts | 4 | ~11 KB |
| Python Tools | 2 | ~13 KB |
| Documentation | 4 | ~27 KB |
| **Total** | **18** | **~69 KB** |

Output files (after simulation) will be much larger:
- Trajectory (.xtc): ~50-200 MB
- Energy (.edr): ~5-10 MB
- Log files: ~1-5 MB

---

## Common Questions

**Q: Which file should I read first?**
A: Start with [QUICKSTART.md](QUICKSTART.md)

**Q: The simulation failed, where do I look?**
A: Check the last lines of the .log file. See README.md troubleshooting section.

**Q: How do I extend the simulation to 50 ns?**
A: Edit md.mdp, change `nsteps = 25000000`, then regenerate with grompp. See SIMULATION_NOTES.md.

**Q: Can I change the pH or protonation state?**
A: Yes, but requires regenerating phytic_acid.pdb and .itp. See SIMULATION_NOTES.md for details.

**Q: What force field should I use?**
A: CHARMM36 is recommended for phosphate groups. See README.md for alternatives.

**Q: How do I cite this work?**
A: See README.md and SIMULATION_NOTES.md for references.

---

## Updates and Versions

**Version**: 1.0
**Created**: 2026-01-01
**Force Field**: CHARMM36 (July 2022)
**Gromacs**: 2020+ compatible
**pH**: 7.4 (physiological)
**Charge**: -10

---

## Getting Help

1. **Documentation**: Read README.md thoroughly
2. **Gromacs Manual**: https://manual.gromacs.org/
3. **CHARMM-GUI**: https://www.charmm-gui.org/
4. **Tutorials**: http://www.mdtutorials.com/gmx/
5. **Forums**: https://gromacs.bioexcel.eu/

---

*For a complete understanding of the system, read all three main docs:*
*QUICKSTART.md → README.md → SIMULATION_NOTES.md*
