# Simulation Notes - Phytic Acid at pH 7.4

## System Specifications

### Molecule Details
- **Name**: Phytic acid (myo-inositol hexakisphosphate, IP6)
- **Formula**: C₆H₈O₂₄P₆ (at pH 7.4)
- **Molecular Weight**: 649.950 g/mol
- **Net Charge**: -10 (at pH 7.4)
- **Total Atoms**: 44 (6C + 24O + 6P + 8H)

### Protonation State
At physiological pH (~7.4):
- **10 deprotonated phosphate oxygens** → -10 charge
- **2 protonated phosphate groups** (monoprotonated)
- This represents a realistic biological state

**Reference**: Fully protonated form (C₆H₁₈O₂₄P₆) = 660.029 g/mol
**Mass difference**: -10.079 g/mol (10 fewer H atoms)

### Force Field Parameters
- **Force Field**: CHARMM36
- **Water Model**: TIP3P (CHARMM-compatible)
- **Long-range Electrostatics**: PME (Particle Mesh Ewald)
- **Cutoffs**: 1.2 nm (Coulomb and vdW)

### System Composition (After Setup)
- **Phytic acid molecules**: 1
- **Water molecules**: ~2000 (estimated for 4 nm box)
- **Na⁺ ions**: ~10 (for charge neutralization)
- **Additional NaCl**: 0.15 M ionic strength
- **Box size**: 4 nm cubic
- **Total atoms**: ~6,500-7,000

## Simulation Protocol

### 1. Energy Minimization
- **Algorithm**: Steepest descent
- **Max steps**: 50,000
- **Tolerance**: 1000 kJ/mol/nm
- **Purpose**: Remove bad contacts and clashes

### 2. NVT Equilibration
- **Duration**: 100 ps
- **Temperature**: 300 K
- **Thermostat**: V-rescale (τ = 0.1 ps)
- **Restraints**: Position restraints on phytic acid heavy atoms
- **Purpose**: Equilibrate temperature

### 3. NPT Equilibration
- **Duration**: 100 ps
- **Temperature**: 300 K
- **Pressure**: 1 bar
- **Barostat**: Parrinello-Rahman (τ = 2.0 ps)
- **Restraints**: Position restraints on phytic acid heavy atoms
- **Purpose**: Equilibrate pressure and density

### 4. Production MD
- **Duration**: 10 ns (5,000,000 steps)
- **Timestep**: 2 fs
- **Temperature**: 300 K
- **Pressure**: 1 bar
- **Restraints**: None
- **Output frequency**: Every 10 ps

## Expected Behavior

### Structural Properties
1. **Compactness**: Phytic acid should remain relatively compact due to:
   - Intramolecular hydrogen bonds
   - Electrostatic repulsion between phosphate groups

2. **Hydration**: Expect extensive hydration shell:
   - Each phosphate group can accept 4-6 H-bonds from water
   - Total ~30-40 water molecules in first hydration shell

3. **Ion Coordination**: Strong Na⁺ coordination:
   - Each phosphate can coordinate 1-2 Na⁺ ions
   - Competitive binding between different phosphate groups
   - May observe transient chelation complexes

4. **Dynamics**:
   - Phosphate groups: High flexibility, large fluctuations
   - Inositol ring: Rigid, chair conformation maintained
   - Overall RMSD: Expected 0.2-0.5 nm

### Physical Chemistry Considerations

**pKa Values of Phytic Acid:**
- pKa₁, pKa₂ ≈ 1.5-2.0 (always deprotonated at pH 7.4)
- pKa₃-pKa₆ ≈ 2.5-5.5 (always deprotonated at pH 7.4)
- pKa₇-pKa₁₀ ≈ 6.5-9.5 (partially protonated at pH 7.4)
- pKa₁₁, pKa₁₂ ≈ 10-12 (mostly protonated at pH 7.4)

**Current Model**: 2 protonated sites → -10 charge (reasonable for pH 7.4)

## Analysis Checklist

After simulation, check:

### Stability
- [ ] RMSD plateaus after equilibration
- [ ] Temperature stable at 300 ± 5 K
- [ ] Pressure fluctuates around 1 ± 200 bar
- [ ] Density stable at ~1000 kg/m³
- [ ] Total energy conserved (no drift)

### Structure
- [ ] Inositol ring maintains chair conformation
- [ ] No phosphate group protonation changes (fixed in classical MD)
- [ ] Radius of gyration stable
- [ ] SASA relatively constant

### Solvation
- [ ] RDF shows well-defined water peaks around phosphates
- [ ] First hydration shell at ~0.25-0.35 nm
- [ ] Hydrogen bond count stable (30-50 total)
- [ ] Na⁺ ions strongly coordinated

### Specific Analyses to Run

```bash
# 1. Coordination number of Na+ around phosphates
gmx rdf -s md.tpr -f md.xtc -ref "resname IP6 and name P*" -sel "resname NA"

# 2. Hydrogen bonds between phytic acid and water
gmx hbond -s md.tpr -f md.xtc -num hbnum.xvg

# 3. Inositol ring planarity (distance matrix)
gmx distance -s md.tpr -f md.xtc -select "resname IP6 and name C*"

# 4. Phosphate-phosphate distances
gmx pairdist -s md.tpr -f md.xtc -ref "resname IP6 and name P1" -sel "resname IP6 and name P4"

# 5. Water residence time around phosphates
gmx trajectory -s md.tpr -f md.xtc -select "resname SOL within 0.35 of resname IP6"
```

## Known Limitations

### Force Field Limitations
1. **Fixed protonation**: Cannot model pH-dependent protonation changes
2. **Polarization**: Non-polarizable force field (CHARMM36)
   - High charge density may require polarizable models for accuracy
3. **Metal coordination**: Classical force field may not capture all Na⁺-phosphate interactions

### Simulation Length
- **10 ns**: Good for initial stability assessment
- **50-100 ns**: Better for conformational sampling
- **100+ ns**: Recommended for publication-quality results

### System Size
- **4 nm box**: Adequate for single molecule
- Minimum distance from phytic acid to box edge: ~1.5 nm
- Box is large enough to avoid self-interactions

## Extending the Simulation

### For Better Sampling
1. Increase production run to 50-100 ns
2. Run multiple replicas with different initial velocities
3. Perform umbrella sampling for conformational free energies

### For pH Effects
1. Use constant-pH MD (requires specific Gromacs compilation)
2. Simulate multiple protonation states separately
3. Calculate relative free energies

### For Ion Binding
1. Increase concentration of Ca²⁺ or Mg²⁺ (strong phytate chelators)
2. Calculate binding free energies using alchemical methods
3. Analyze coordination geometries and residence times

### For Biological Context
1. Add protein (e.g., phytase enzyme)
2. Simulate membrane binding (phytic acid can interact with lipids)
3. Include additional metabolites or cofactors

## Troubleshooting Tips

### If energy minimization fails:
- Structure may have clashes → check PDB visually
- Increase emtol to 5000
- Check topology charges sum correctly

### If NVT temperature unstable:
- Reduce tau_t (more aggressive coupling)
- Check that position restraints are applied
- Verify gen_vel = yes

### If NPT pressure oscillates wildly:
- Normal for small systems (large fluctuations expected)
- Extend equilibration time
- Check compressibility value

### If production crashes:
- LINCS warnings → reduce timestep to 1 fs
- Check for exploding atoms (coordinates > 999 nm)
- Verify checkpoint files exist

## References & Further Reading

### Phytic Acid Structure
- Graf, E. (1983) "Applications of phytic acid" J. Am. Oil Chem. Soc. 60, 1861-1867
- Cowieson et al. (2011) "Phytate and microbial phytase" J. Nutr. 141, 1331-1336

### MD Simulation
- GROMACS tutorial: http://www.mdtutorials.com/gmx/
- CHARMM36 paper: Best et al. JCTC 2012, 8, 3257-3273
- Phosphate parameters: Klauda et al. J. Phys. Chem. B 2010, 114, 7830-7843

### Phytic Acid Chemistry
- Torres et al. (2005) J. Agric. Food Chem. 53, 2618-2623
- Bohn et al. (2008) "Phytate: Impact on environment and human nutrition" Crit. Rev. Food Sci. Nutr. 48, 747-784

---

**Last Updated**: 2026-01-01
**Force Field**: CHARMM36 (July 2022 or later)
**Gromacs Version**: 2020+ recommended
