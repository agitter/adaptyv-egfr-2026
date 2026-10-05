# Continuation state: EGFR stage 3

Read REPORT.md, METHODS_AGE_AUDIT.md, and output/stage3_summary.json first.
This is an interim checkpoint, NOT the final ranked 100-sequence delivery.
Phone-a-Friend usage: 2 of 3. No third request has been made.

## Current scientific state

1,114 sequence/structure records, 1,097 distinct sequences, 1,044 CDR3s. Do not
count the 17 duplicate structural hypotheses as distinct final sequences.
The G00117_C04 family, especially C00003 and C00009, has the best expanded
sensitivity evidence. C00017/C00024 also have promising default and mouse
results, but do not yet have the same acid-state, capped-target and alternative-
geometry coverage. C00003/C00009 share a flexible-loop family: parent G00117_C04
moves substantially in a short unbound trajectory. Do not claim folding or
preorganization from the local-minimization passes.

Detailed human calculations: 49 native geometries, 27 acid-on in all 54 model
assumptions. Mouse: 7. Geometric checks pass for all 49 refined complexes;
G00012_R00 has one peptide-omega flag. A00008 has a glycosylation sequon and must
not be promoted. A-series pH improvements generally lose favorable human acidic
interaction energy. G00053_C05 is geometry/conformation sensitive.

C00003: own B00013/G00117 family, K36Q/K105Q. C00009 adds N52E. Acid-stress tests
include both neutral carboxylate forms, 81 and 243 microstates respectively.
All free-pKa values and state priors remain assumptions. Use output JSONs for
exact final values, not recollection of preliminary logs.

## Runtime and limits

Container CPU quota is 8 and hard memory cap approximately 4 GiB. Unrelated
artifact_rpc processes exist and must not be killed. Use at most 4-6 single-
threaded OpenMM processes, and check memory before large concurrent jobs.
No GPU, external download DNS unavailable. Browser research works. All needed
original inputs for the current workflow have been supplied.

```sh
R=/mnt/data/egfr_campaign
S=$R/stage3
export PYTHONPATH=$S/code:$R/stage2/code:$R/code:$R/stage2/runtime/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
```

The archive contains everything small and campaign-generated needed for this
state. Large uploaded inputs, runtime binaries/wheels, BLAST indexes and expanded
reference databases are excluded with hashes. Restore the original archives:

- resources.zip -> inputs/provided (consult existing source metadata; do not
  blindly duplicate internal nesting).
- egfr_inputs.zip -> inputs/transfer1, retaining legacy_fetch/.
- egfr_refinement_inputs.zip -> inputs/transfer2, retaining egfr_refinement_fetch/.

Then, if missing, rebuild using stage2/code/recover_runtime.py,
normalize_antibodies.py, augment_novelty.py and build_blast_indexes.sh. Construct
an augmented antibody index with:

```sh
$R/stage2/runtime/bin/makeblastdb \
  -in $R/stage2/intermediate/novelty/antibodies_augmented.fasta \
  -dbtype prot -out $R/stage2/intermediate/blast/antibodies_augmented
```

Use source and archive hashes in the checkpoint manifests. Original stage-1
large checkpoint was absent in this runtime; stage 2 explicitly reconstructed
what was needed. The earlier PCG64 pilot is excluded from candidate ancestry.

## Key files

- intermediate/designs/{id}.json: sequence, CDRs, atom coordinates, origins,
  parent edits, seeds, score-validity labels.
- intermediate/evaluation/{id}.json: initial geometry/sequence/coarse gates.
- intermediate/refined/{id}_{human6ARU|mouseAF}_*: actual all-atom PDBs, audits,
  explicit proton microstate energies and conditional results.
- intermediate/acid_ensemble/: carboxylate plus histidine stress tests.
- intermediate/capped/: target termini capped with ACE/NME.
- intermediate/common_sites/: cutoff sensitivity; lead four-His audit included.
- intermediate/controls/fcrn/: historical independent positive control.
- intermediate/unbound/: local free-binder relaxation, parent short trajectory.
- intermediate/surface/: whole-target buried area and alternative-receptor clashes.
- intermediate/novelty/: full-loop screen, complete selected BLAST output,
  explicitly partial all-library antibody search retained as PARTIAL.
- output/: audit-derived comparisons and unranked NOT_FOR_SUBMISSION FASTA.

The final all-library CDR3 search has 203,818 reference segments from 446,966
unique antibody sequences. All records below70%, max64.706%. Whole-chain BLAST
completed for 61 records/60 sequences against each of Swiss-Prot, PDB and
augmented antibodies. Do not imply whole-chain screening of all1,114 completed.
The organizer's exact annotation/database pipeline is not reproduced.

## Commands and pitfalls

Rebuild summaries without rerunning science:

```sh
python "$S/code/stage3_audit.py"
```

Evaluate one existing candidate against human or mouse:

```sh
python "$S/code/refine_new.py" C00003 human6ARU
python "$S/code/refine_new.py" C00003 mouseAF
```

refine_new.py may reuse cached complete outputs. Do not overwrite scientific
comparisons casually. The batch runner accepts JSON entries that are strings
or dictionaries {"candidate_id":"...","species":"mouseAF"}; NOT two-item
arrays. An earlier invalid array launch is retained in logs and reference.

```sh
python "$S/code/acid_histidine_ensemble.py" C00003 54,103
python "$S/code/acid_histidine_ensemble.py" C00009 52,54,103
python "$S/code/refine_state.py" C00003
python "$S/code/cap_target.py" C00003
python "$S/code/free_binder.py" C00003
python "$S/code/interface_surface.py" C00003 192
```

Review each script's main signature before invoking it; do not infer argument
contracts from a similar script. free_binder.py optional steps must currently
be a multiple of 1000 and writes its final audit only after the entire trajectory.
Its no-thermostat velocity initialization is not a 300 K equilibration protocol.
A future dynamics study needs staged equilibration and matched controls before
interpreting flexibility quantitatively. The existing 10-ps trajectory remains
valid as a disclosed exploratory warning, not a stability certificate.

G/B/C variants with inherited or placeholder coarse scores are explicitly
excluded by the new coarse_rank_eligible guard. Default evaluate_new.py
selection is an EXPLORATORY queue, not the final rank. All manual selections,
failed inputs, and pre-patch source versions are preserved.

## Immediate next stage

Keep the two-His recognition hypothesis but improve loop preorganization,
particularly the mobile second CDR. Test only geometrically compatible Pro or
Gly replacements, intrabinder packing/H bonds, or viable extra disulfide
constraints using methods predating2011. Preserve all contact residues until
alternative contacts are explicitly evaluated. Do not assign an entropy bonus
and call it experimental folding evidence. Compare free and bound state
geometries, avoid changing neutral binding favorability unnoticed, and retain
alternative-pose diversity rather than filling the final20 with one motif.

Historical references for this next stage (NOT methods already executed):
Matthews/Nicholson/Becktel1987, DOI10.1073/pnas.84.19.6663; Saerens et al.2008,
"Disulfide Bond Introduction for General Stabilization of Immunoglobulin
Heavy-Chain Variable Domains" (verify precise bibliographic details before
final citation); interloop disulfide in antibody mimics,2007,
DOI10.1016/j.jmb.2007.02.029. These support considering constrained loops, not
assuming their effect transfers to these sequences.

Expand distinct-sequence detailed coverage, perform the same novelty audit,
and produce the final100 ranked FASTA plus CSV and top20 subset. Final header
names should be opaque biochemical puns; keep the design-ID/public-name key
private in provenance, not in public FASTA descriptions. Never submit100 to a
track allowing20. Do not promise that an untested sequence satisfies biological
criteria. State the evidence and uncertainty in the final delivery.

## Completion and provenance

The postflight JSON records no running campaign jobs and the memory OOM counters.
The archive checksum and independent per-file verification accompany delivery.
Large inputs are referenced, not duplicated; retain the original three uploads.
Scientific method classes satisfy the historical constraint in the current
lineage under the user's modern-implementation allowance. Preserve the excluded
PCG64 pilot disclosure rather than issuing a false campaign-wide blanket claim.
