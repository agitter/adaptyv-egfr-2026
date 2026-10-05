# Stage 3 historical-method audit

Scope: scientific design, physical scoring, and sequence-screening operations in
this stage and the current candidate lineage. The user explicitly permits modern
programming languages, implementations of historical tools, and modern datasets.
This is not a claim that the operating system, compiler, chat model, or every
software dependency existed in 2010.

| Operation | Historical basis | Implementation here |
| --- | --- | --- |
| New loop backbones | Cyclic coordinate descent for protein loops, 2003 | Own implementation from stage 2; all three original framework CDRs removed |
| Coordinate superposition / rigid docking | Classical rigid geometry and least-squares superposition, pre-2011 | NumPy/SciPy coordinate calculations; no learned docking model |
| Sequence and rotamer search | Empirical rotamers, Metropolis sampling (1953), simulated annealing (1983), fixed-backbone interface design | Own one-/two-body energy search; current PDB residues used only as rotamer data |
| Positive / negative multistate design | Havranek and Harbury, 2003, DOI 10.1038/nsb877 | Human acidic, human neutral, mouse, intrabinder and glycan terms |
| Pseudorandom sampling | MT19937, 1998 | Explicit NumPy MT19937 generator, seeds retained |
| All-atom molecular mechanics | Amber99SB, 2006, DOI 10.1002/prot.21123 | OpenMM 8.4, permitted modern implementation |
| Implicit solvent and desolvation | OBC generalized Born, 2004; classical Coulomb/Debye, Lennard-Jones, surface-area models | Exact matched-partner MM/GBSA proxy; independent OpenMM arithmetic checks |
| Protonation / tautomer enumeration | Binding polynomials and Wyman linkage, 1964; continuum-electrostatic microstate methods before 2011 | Own finite-state partition sums; HID, HIE, HIP and two neutral carboxylate states |
| Local minimization | Classical gradient-based molecular-mechanics minimization, pre-2011 | OpenMM local minimizer; classical chirality / geometry restraints |
| Short unbound trajectory | Verlet, 1967; constrained molecular dynamics, pre-2011 | Deterministic Verlet, explicit MT19937 velocities; no stochastic learned potential |
| Solvent-accessible area | Bondi radii, 1964; Shrake-Rupley, 1973 | Own sphere-point calculation with 192/512-point sensitivity check |
| Novelty sequence search | BLAST, 1990; gapped BLAST, 1997; Levenshtein edit distance, 1965 | BLAST+ 2.17 under user's explicit modern-implementation allowance; exact loop edit-distance comparisons |
| Controls | FcRn/Fc structure 1I1A and mechanism, 2001 | Calibration/rejection control only, not a design seed |
| Framework | Universal humanized nanobody framework publication, 2008; 3EAK deposited/released 2008 | Framework only; original CDR sequences and coordinates not retained as designed loops |

## Important limits on certification

No modern learned protein generator, sequence model, structure predictor, learned
potential, or post-2010 binder-design algorithm was run in stage 3. Current
UniProt, PDB, PLAbDab and therapeutic-antibody datasets, the supplied AlphaFold
mouse coordinates, and the supplied modern human coordinates are used as data,
not as new model inference.

The original exploratory campaign included a PCG64-based sampler. That was
identified and explicitly excluded; stage 2 regenerated the loop/docking lineage
with MT19937. Stage 3 derives from those regenerated backbones, not from the
excluded pilot. Preserve this disclosure. An unqualified statement that literally
*every operation ever executed in this conversation* met the historical cutoff
would be false, even though the current design lineage uses the eligible methods.

The new code and combinations of historical calculations were written in 2026;
the historical cutoff applies to scientific method classes, with implementation
versions permitted by the user. Software execution does not validate the model's
biological accuracy.

## Sources

Primary-source identifiers and current rule checks are retained in
`reference/plan.json` and `reference/live_rules_and_literature_check.json`.
Relevant historical identifiers include:

- Canutescu and Dunbrack (2003), cyclic coordinate descent, DOI 10.1110/ps.0242703.
- Havranek and Harbury (2003), DOI 10.1038/nsb877.
- Hornak et al. (2006), DOI 10.1002/prot.21123.
- Wyman (1964), DOI 10.1016/S0065-3233(08)60190-4.
- Georgescu, Alexov and Gunner (2002), DOI 10.1016/S0006-3495(02)73940-4.
- Mongan, Case and McCammon (2004), DOI 10.1002/jcc.20139.
- Martin et al. (2001), FcRn/Fc, DOI 10.1016/S1097-2765(01)00230-1.
- Vincke et al. (2008), humanized VHH scaffold, DOI 10.1074/jbc.M806889200.

MCCE and constant-pH molecular dynamics are historical precedents, not programs
run here. Their published validation accuracy must not be attributed to this
custom fixed-geometry enumeration.
