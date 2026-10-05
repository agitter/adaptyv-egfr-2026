# EGFR design campaign - stage 2 checkpoint

**Date:** October 1, 2026. **Status:** sequence-generation and falsification checkpoint, not a final submission. **Phone-a-Friend:** two used, one unused.

## Result

There are now **366 actual, independently initialized sequence designs**, each with binding-loop sequences, atom coordinates, selected empirical rotamers, random seeds, and scores. They are not merely docked alanine placeholders. The library contains 346 unique CDR3 sequences across 41 docking poses; sequence lengths are 119-127 residues.

The desired final deliverable remains 100 ranked sequences, with a top-20 competition subset. It has not been released because the new tests exposed meaningful false positives in the initial pH ranking. No sequence in this checkpoint has measured affinity, measured pH selectivity, or demonstrated mouse binding. Positive computational scores do not establish the requirement of no detectable human binding at pH 7.4.

## Independent checks completed

| Check | Result |
|---|---:|
| Fresh classical loop construction | 14,400 attempts; 152 accepted individual loops |
| Fresh rigid-body docking | 108,000 orientations; 474 retained poses |
| Sequence/side-chain optimization | 366 distinct sequences on 41 poses |
| Human/mouse/glycan geometry prefilter | 211 pass; 155 rejected |
| Acid-on direction in all 12 coupled-protonation scenarios | 185 of the 211 geometry survivors |
| Geometry plus robust coarse direction plus no N-X-S/T sequon | 127 |
| Local expanded CDR3 novelty test | All 366 below 70% edit identity; maximum 61.11% |
| Classical scoring/geometry unit checks | 264 pass |
| Coupled binding-polynomial checks | 132 pass |
| Separated-partner MM/GBSA subtraction control | Pass; residual -0.0191 kcal/mol at 1,000 nm separation |

These gates test different things. In particular, the 185 positive-direction results are **not** 185 successful pH-switch binders, and the 127 additional sequence-liability survivors are **not** a validated final candidate pool.

The 155 geometry rejections include overlapping counts of 128 glycan clashes, 27 human-protein clashes, and 18 mouse-protein clashes. A covalent-link audit separates EGFR sugars from antibody sugars in the human reference: 150 target-glycan atoms are retained; 14 antibody-linked atoms are excluded. The initial sequence search had not screened side-chain/glycan collisions directly, so these failures are preserved as a corrective finding.

An additional receptor-conformation test found that all 211 present geometry survivors collide with domain I in ligand-bound 1IVO. Of those, 154 remain clash-free against the unliganded 1NQL conformation. They are conformation-specific hypotheses; active-state compatibility must not be claimed.

## The most important result: the coarse pH model was too encouraging

Candidate P00024_07, the strongest worst-case coarse pH candidate, was examined more deeply. Two monomer controls had already passed geometry checks. Joint minimization of this candidate and human domain III, with fixed EGFR glycan obstacles, then passed the independent audit: 260 alpha stereocenters, 43 beta stereocenters, no bond-length outliers, and no severe heavy-atom clashes.

Despite that geometrically clean result, the matched-component Amber99SB/OBC diagnostic did not support an affinity claim. Its neutral-reference interaction proxy was **+21.77 kcal/mol**, not a calibrated binding free energy. Protonation effects also differed from the initial coarse-model impression:

| Protonated histidine | Change in binding-energy proxy, kcal/mol | Diagnostic direction |
|---|---:|---|
| Human receptor H370 | -1.17 | Favors acid-on |
| Human receptor H433 | -6.22 | Favors acid-on |
| Binder H108 | +1.91 | Opposes acid-on |
| Binder H109 | +3.23 | Opposes acid-on |

Each number compares complex-minus-partner energies at matched coordinates and protonation. The separated-partner control tests this subtraction, not predictive accuracy. The earlier, less reliable monomer-only probe is also retained; it had a much larger interface-strain artifact. No raw energies with different proton counts were treated as pH-dependent affinities.

**Decision:** do not promote this candidate merely because it leads the coarse pH ranking. Test and redesign the responsible contacts. A clean structure is necessary but insufficient.

## Novelty coverage

The expanded local collection contains **446,966 unique antibody sequences and 203,818 CDR3 segments**. Validated numbering supplied by the database was used as data; no modern antibody-numbering inference was run. Every generated CDR3 was compared to the collection using normalized Levenshtein distance. The maximum nearest-reference identity was 61.11%.

Broad BLAST searches against Swiss-Prot (575,748 sequences) and PDB (1,165,667 sequences) completed for all 366 candidates. Among returned hits covering at least 80% of a query, maximum identities were 75.00% and 74.194%, respectively. Shared antibody frameworks are expected, so these are screening statistics, not a binary organizer-eligibility decision. The broad search against the older, smaller antibody index was stopped at the checkpoint; it should be replaced with a search against the augmented collection after redesign. Its incomplete status is explicit in `reference/blast_status.json`.

Local annotation, reference coverage, and thresholds are not an exact reproduction of the organizers' novelty pipeline. Upload acceptance is not certified.

## Revised scientific direction

The evidence favors testing interfaces in which **acidic binder groups contact protonatable receptor histidines**, rather than relying on a cluster of introduced binder histidines. New contact patterns must retain favorable packing and hydrogen bonds; counting charges is not sufficient.

The target alignment identifies conserved receptor H358, H370, and H433. Human H383 is arginine in mouse. That substitution might support acid-on recognition of human EGFR and constitutive recognition of mouse EGFR, provided the larger arginine is sterically accommodated. This is an untested design hypothesis, not a demonstrated binding mechanism. Mouse pH switching is not required, so future optimization should not penalize stronger mouse binding with a symmetric species-difference term.

The next stage is to place target glycans directly in rotamer selection, benchmark and extend the all-atom diagnostic, optimize receptor-histidine-gated interfaces, evaluate multiple target conformations, and rerun sequence/novelty filters before final ranking.

## Historical methods and checkpoint integrity

The new candidate lineage uses MT19937, CCD loop closure, empirical rotamers, rigid docking, simulated annealing, classical multistate design, binding polynomials, Coulomb/Born/SASA approximations, Amber99SB/OBC, Levenshtein distance, and BLAST. Those method families predate 2011. Modern implementations and supplied modern datasets are used under the user's explicit allowance; no contemporary binder generator or new learned folding inference was run.

An inherited PCG64 exploratory pilot is explicitly excluded from the new sequence lineage. Because that earlier pilot did run, an unqualified claim about every operation in the entire conversation would be inaccurate; the final method audit must disclose the correction.

The previously linked large stage-1 archive was absent from the reset runtime. This stage reconstructed required functionality from the retained earlier checkpoint and original user inputs rather than claiming to have recovered missing files. The compact checkpoint now includes all current source, small intermediates, real sequence/structure designs, failed diagnostics, test results, logs, and file manifests. Large immutable inputs and rebuildable runtime/database/index files are represented by hashes and restoration instructions, not duplicated. The archive-only 366-sequence FASTA is clearly marked **NOT FOR SUBMISSION**.

All campaign workers have completed or been stopped; nothing is left running after this checkpoint. Detailed continuation instructions are in `stage2/STATE.md`.
