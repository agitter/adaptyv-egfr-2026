# EGFR campaign: stage 2 checkpoint, 2026-10-01

## Read this first

This is a working checkpoint, NOT a final binder submission. There are 366 actual sequence/structure hypotheses, but none is experimentally validated and the best coarse-model candidate did not pass a convincing all-atom affinity/pH interpretation. The final ranked 100-sequence FASTA, public codenames, competition CSV, top-20 subset, and final method certification have NOT been produced. Do not submit the internal 366-record FASTA.

Phone-a-Friend requests used: 2 of 3. No third request was used in stage 2.

## Runtime restoration and provenance

The previously linked `egfr_stage1_checkpoint.zip` was NOT present in this reset runtime. Only its report was retained. We restored the earlier `egfr_campaign_checkpoint.zip`, `egfr_request2_audit.zip`, and both user-supplied input archives. Stage 2 therefore independently reconstructs required functionality rather than claiming a bit-identical stage-1 continuation. All current new work is under `stage2/`; inherited pilot work remains in the parent directory.

The compact stage-2 checkpoint includes source, small source inputs, new loops/poses/sequences/structures, controls, failures, evaluations, logs, and manifests. It excludes recoverable software wheels/binaries, compressed databases, expanded reference databases, BLAST indexes, and caches. Those are covered by the excluded-file and original-archive hash manifests and regeneration scripts. Do not delete the original input ZIPs.

Paths used by current scripts:

```
R=/mnt/data/egfr_campaign
S=$R/stage2
export PYTHONPATH=$S/code:$R/code:$S/runtime/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
```

To reconstruct missing runtime/reference files, first extract the original transfer archives under `inputs/transfer1` and `inputs/transfer2`, preserving their internal folder names; then run:

```
python "$S/code/recover_runtime.py"
python "$S/code/normalize_antibodies.py"
python "$S/code/augment_novelty.py"
sh "$S/code/build_blast_indexes.sh"
```

The first two BLAST indexes are built by the shell script. For the expanded antibody index, explicitly use `stage2/intermediate/novelty/antibodies_augmented.fasta` with `makeblastdb -dbtype prot`; the earlier completed/partial broad BLAST searches used the INITIAL antibody index, not the augmented index. Read `reference/blast_status.json` for final search completion status.

## Completed new scientific branch

- Generic humanized VHH framework from 3EAK chain A. The three complete original binding-loop regions were removed. No donor CDR sequence or coordinates initialized new loops.
- Fresh MT19937-based torsion sampling and cyclic coordinate descent: 14,400 attempts, 152 accepted individual loops.
- Fresh docking: 450 assembled backbones from 5,713 assembly trials, 108,000 orientations, 474 retained poses.
- Classical pairwise atom-contact, hydrogen-bond, electrostatic, packing and multistate simulated-annealing sequence search: 41 poses, 366 unique sequences, 346 unique CDR3s, lengths 119-127.
- Alpha chirality, peptide bonds, human/mouse target contacts and crystallographic glycan checks: 211 survive. 155 fail, with overlapping counts of 128 glycan, 27 human protein, and 18 mouse protein clashes.
- Full coupled histidine binding polynomials use the complex partition function divided by both unbound-partner partition functions. They include intrapartner terms, approximate charge burial, and histidine coupling. Twelve parameter scenarios: 185 of 211 have the acid-on DIRECTION in every scenario; 26 do not. Direction alone is not the required experimental switch magnitude or neutral-pH nonbinding.
- Of the 185, 127 also lack N-X-S/T sequons. Other liabilities are recorded, not silently ignored.
- 264 classical tests and 132 binding-polynomial tests passed. The FcRn-IgG 1I1A positive control has the expected acid-on direction in all twelve coarse scenarios, but this does not calibrate EGFR affinity or validate switch magnitudes.
- All 211 human/mouse geometry survivors clash with domain I of the ligand-bound 1IVO arrangement; 154 are clash-free in the unliganded 1NQL arrangement. Current poses are conformation-specific hypotheses, not established binders across receptor states.
- The glycan audit follows actual covalent links: 150 EGFR glycan atoms retained, 14 antibody-linked sugar atoms excluded. The initial sequence optimizer did NOT explicitly exclude glycan side-chain clashes; this was caught by independent screening and is a required next-stage correction.

## Novelty checks

Expanded references contain 446,966 unique antibody sequences and 203,818 CDR3 segments. The expansion uses provided, validated database numbering as DATA; no modern antibody-numbering model was run. 50,466 aligned heavy-chain rows failed the annotation validation rule and were not falsely accepted as annotated CDR3s; raw sequences remain available through other screening routes.

All 366 generated CDR3s are below 70% normalized Levenshtein identity to this local CDR3 collection. The highest nearest-reference value is 61.1111%. This is not an exact reproduction of IMGT or the organizer's database snapshots, and is not eligibility certification. Swiss-Prot and PDB broad searches supplement, not replace, the antibody CDR screen. Read the final BLAST status file rather than assuming a running job completed.

## All-atom findings: do not bury these

Amber99SB/OBC and classical restrained minimization were used. Monomer controls P00001_01 and P00024_07 passed independent heavy-atom geometry checks. The jointly relaxed P00024_07/human domain-III complex passed 260 alpha and 43 beta stereocenters, bond bounds, and severe-clash checks. The target crop is 334-505; 150 fixed target-glycan obstacles were included during joint repair.

The monomer-only docking diagnostic had a large steric artifact, so a joint-complex refinement was necessary. The joint MM/GBSA diagnostic is retained alongside the failed/less useful earlier diagnostic. It uses identical protonation and coordinates in complex and isolated partners; it does NOT compare raw energies of molecules with different proton counts.

For the jointly relaxed candidate:

- Neutral-reference MM/GBSA proxy: +21.77 kcal/mol. This is not a measured binding free energy or a calibrated KD; it does not support a strong-affinity claim.
- Protonating target H370 changes the binding proxy by -1.17 kcal/mol.
- Protonating target H433 changes it by -6.22 kcal/mol.
- Protonating binder H108 changes it by +1.91 kcal/mol.
- Protonating binder H109 changes it by +3.23 kcal/mol.

Thus two engineered binder histidines oppose acid-on binding in this diagnostic, contrary to the encouraging initial coarse ranking. The candidate must NOT be promoted to a final lead solely on its coarse pH score. Neither geometry nor this single-site MM/GBSA probe establishes the final experimental criteria. Missing ingredients include reliable bound/unbound ensembles, complete protonation coupling/tautomer sensitivity, entropy, salt and glycan uncertainty, and experimental validation.

## Required next stage

1. Audit the MM/GBSA calculation using matched-coordinate separated-partner and historical FcRn controls, then test more than one candidate. Preserve failures.
2. Put EGFR-attached glycans and all fixed framework heavy atoms directly in rotamer selection and docking exclusion, instead of relying on late rejection.
3. Create a receptor-histidine-gated branch: use acidic binder groups to recognize protonated receptor histidines; disfavor added binder histidines unless all-atom diagnostics support the desired sign. Repack aromatic/hydrophobic and polar contacts rather than optimizing charge count alone.
4. Explore patches spanning H358/H370/H433 and H370/H383/H433. H358, H370 and H433 are conserved as histidine in mouse. Human H383 is mouse arginine: this may offer useful acidic-human/constitutive-mouse recognition if sterically compatible. It is a hypothesis, not a demonstrated mechanism. Mouse pH switching is not a competition requirement, so do not penalize stronger mouse binding by a symmetric species-difference penalty.
5. Compare alternative target conformations explicitly. Binding the unliganded state may be relevant, but do not claim active-state compatibility for the present poses.
6. Advance only candidates with repeatable geometry, appropriate pH direction across models, plausible contact energetics, mouse compatibility, and novelty/sequence checks. Complete/refine the broad novelty searches against the augmented antibody collection before final selection.
7. Only then select 100 diverse, ranked candidates, assign cryptic public codenames, emit FASTA/CSV and top-20 files, and write a qualified method audit. Do not invent KD values, pH selectivity measurements, or success probabilities.

## Historical-method boundary

Scientific method families in the new lineage are pre-2011: classical covalent geometry and molecular mechanics; torsion/rotamer sampling; CCD loop closure (2003); rigid docking; simulated annealing; multistate positive/negative design; Wyman binding polynomials; Coulomb/Debye/Born/SASA models; Amber99SB (2006); OBC (2004); Levenshtein (1965); BLAST (1990); MT19937 (1998). Modern language/library/compiler releases are allowed implementations, not modern learned design methods. Supplied modern structures/databases are permitted data; no new AlphaFold or modern binder-generator inference was run.

Inherited exploratory PCG64 pilot outputs are explicitly EXCLUDED from current sequence ancestry. Because that older pilot did execute, do not issue an unqualified claim that literally every historical operation in the entire conversation used a pre-2011 algorithm. The current branch was regenerated to fix this. The final audit must disclose it.

`reference/framework_definition.json` inherited planned fields `framework_FERA_option` and `terminal_append`; these are not records of applied mutations. Current sequences use the unmodified retained framework and its existing VTVS terminus. There is no extra appended serine.
