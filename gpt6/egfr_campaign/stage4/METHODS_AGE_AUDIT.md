# Stage 4 historical-method audit

## Scope

Scientific method classes used to generate and evaluate the current candidate
lineage are from 2010 or earlier. Modern implementations, Python, numerical
libraries and modern datasets are permitted by the user's clarification. New
source files were written during this campaign; this is not a claim that those
files, the chat model, operating system or every dependency existed in 2010.

| Operation | Historical basis | Use in this stage |
| --- | --- | --- |
| Generic nanobody framework | Vincke et al., 2008; PDB3EAK | Framework inherited. All original binding-loop sequences and coordinates had been removed in candidate ancestry. Native loops used only in control trajectories. |
| New loop geometry | CCD protein-loop construction,2003 | Inherited from the regenerated MT19937 lineage; not a new generator here. |
| Sequence changes | Geometry-compatible Gly-to-Ala/X-to-Pro entropy strategy, Matthews et al.,1987 | Deterministic enumeration; no theoretical entropy bonus treated as measured stability. |
| Additional disulfides | Geometric cystine design and antibody stabilization before2011, including Saerens2008 | Finite rotamers, bond/angle/torsion geometry, Powell1964 optimization. Correct oxidation is not established by a modeled bond. |
| Atom completion and packing | Empirical rotamers and classical nonbonded energies, pre2011 | Historical methods applied to permitted modern structural data. |
| All-atom scoring | Amber99SB2006; OBC generalized Born2004; Coulomb/Debye and van der Waals physics | OpenMM8.4 is a permitted modern implementation. No learned potential. |
| Local refinement | Gradient-based molecular-mechanics minimization, pre2011 | Matched additional refinement of parent controls; classical chirality restraints. |
| Protonation dependence | Binding polynomials/Wyman1964, continuum microstate methods and Henderson-Hasselbalch weighting | Explicit neutral histidine tautomers and protonated histidine; carboxylate/terminal-cap sensitivities; conditional model only. |
| Additional His probe | Classical donor geometry and multistate negative design, pre2011 | F107H proposed and tested, not presumed acid-on. Matched extra-refinement controls added. |
| Thermal stress tests | Velocity-Verlet1967; Andersen1980; RATTLE1983; Maxwell statistics | Synchronized velocity-Verlet with projected constraint-group collisions. MT19937 velocities/collisions. No position restraints or target in production; chirality restraints retained. |
| Pseudorandom sampling | MT19937,1998 | Explicit NumPy MT19937 and compatible standard seeding. |
| Coordinate comparison | Kabsch least-squares superposition,1976 | Framework motion and loop-internal deformation reported separately. |
| Similarity screening | Gapped BLAST1997; BLOSUM621992; Levenshtein1965 | Modern BLAST2.17 and RapidFuzz implementations. Composition adjustment disabled. No modern learned antibody numbering tool run. |
| Proton-linkage bound | Classical positive binding polynomials | Independent numerical property tests, not experimental validation. |

## Exclusions and limitations

No modern protein generator, sequence language model, structure predictor,
learned potential or post-2010 binder-design algorithm was run in stage4. The
supplied modern AlphaFold mouse coordinates are an input dataset, not inference
performed here. PLAbDab/Thera-SAbDab/Swiss-Prot/PDB are screening data, not sources
of copied binding loops.

The original exploratory campaign had a PCG64 sampler. That pilot was explicitly
excluded and the candidate lineage regenerated with MT19937 in stage2. Preserve
that disclosure: a blanket certificate covering literally every exploratory
operation ever attempted would be false.

Three early stage4 trajectories used half-step Verlet velocities with collisions
and exhibited a temperature bias. They are preserved but excluded. The corrected
synchronized integrator passed independent harmonic/constraint checks. Very short
trajectories, artificial chirality restraints and incomplete physical models do
not establish folding, equilibrium stability or binding.

Canonical-loop class constraints (Chothia/Lesk1987) were considered as a next
strategy but not implemented in this stage. MCCE, PROPKA, constant-pH dynamics,
Rosetta and calorimetry were not run. Do not attribute their published validation
accuracy to this custom workflow.

Primary-source identifiers and source-use distinctions are in
`reference/sources.json`; inherited methods also have the stage3 audit.
