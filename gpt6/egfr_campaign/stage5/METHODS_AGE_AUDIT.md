# Stage 5 historical-method audit

## Certification scope

The scientific design and evaluation method classes used in **stage 5** are from 2010 or earlier. The newly written source files implement those historical approaches. Modern Python, numerical libraries, compilation, OpenMM/BLAST implementations, and current input datasets are within the user's explicit allowance. This is not a claim that the source files, chat model, operating system, or every software release existed before 2011.

No modern protein generator, learned structure predictor, sequence language model, learned potential, or post-2010 antibody-numbering/classification method was run. The supplied mouse prediction is an input dataset, not inference performed here.

| Actual operation | Historical basis and scope |
| --- | --- |
| Independent rebuilt loop backbones | Internal-coordinate geometry, non-antibody fragments, CCD2003, Ramachandran1963 and classical steric filters. MT19937 seeds recorded. |
| Generic VHH framework | Vincke2008/3EAK. Original binding-loop coordinates and sequences were removed in the current ancestry; native loops only calibrate controls. |
| Sequence/side-chain design | Finite empirical rotamers, classical sterics/Coulomb/contact heuristics and pre-2011 packing/search ideas; no learned design model. |
| Continuous glutamate proposals | Powell1964 torsion optimization and classical donor geometry. Proposal restraints are absent from physical refinement/scoring. |
| All-atom evaluation | Amber99SB2006/OBC2004 with the user-authorized current OpenMM implementation. Fixed-coordinate energy differences are not calibrated binding energies. |
| Protonation-state enumeration | Henderson-Hasselbalch/Wyman1964 linkage and classical continuum electrostatics. Histidine and carboxylate tautomers explicitly enumerated. Six-site guard expansion changes scope, not algorithm. |
| Independent prior sensitivity | Products of binding-polynomial weights, exhaustive enumeration and log-sum-exp arithmetic, all classical. pKas and free-tautomer populations remain assumptions. |
| Thermal tests | Velocity-Verlet, Andersen1980 collisions, RATTLE1983 constraints, Maxwell velocities and MT19937. No modern thermostat or learned potential. |
| Coordinate comparison | Kabsch1976, signed volumes, bond lengths, peptide and disulfide torsions. |
| Glycan flexibility stress | Rigid rotations about measured bond axes, finite torsion grid, graph-based covalent exclusions and steric distances. Not GLYCAM sampling, learned glycan inference or equilibrium weighting. |
| Novelty | Gapped BLAST1997/BLOSUM621992, exact Levenshtein1965/1966 edit distance, current permitted datasets. No MMseqs2/Foldseek/ANARCI execution. |

## Corrections, excluded history, and unimplemented methods

The earlier shorthand "velocity-Verlet1967" is made more precise here: Verlet's original scheme is1967, while a directly relevant velocity-form precedent is Swope et al.1982 and Andersen's RATTLE1983. Both precede the cutoff. The actual current dynamics algorithm and results were not changed during this bibliographic clarification.

The original exploratory campaign contained a PCG64 pilot. It was excluded and candidate ancestry regenerated with MT19937 in stage2; that provenance remains intact. A blanket certificate that **every historical exploratory operation** complied would be false. Stage4's biased half-step-velocity pilots remain excluded. Stage5 uses the corrected synchronized implementation; independent seeds were declared before trajectories.

Canonical-loop classification was considered but not implemented. MCCE, PROPKA, Rosetta, constant-pH dynamics, learned numbering, folding prediction and experimental binding/stability measurements were not run. Their published accuracy must not be attributed to this workflow.

The library includes failed designs and unvalidated hypotheses. Historical compliance does not certify folding, affinity, species cross-reactivity, neutral-pH nondetection, correct oxidation, expression or official novelty acceptance. Final delivered-sequence ancestry must be checked again when the final100 are selected.

Primary identifiers and source-use distinctions are in reference/source_register.json. Input and output checksums, seeds, code and failures are preserved. No font files or newly downloaded executables are distributed in the checkpoint.
