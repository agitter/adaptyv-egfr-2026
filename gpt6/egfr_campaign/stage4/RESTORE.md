# Restore and reproduce stage4

## Restore data rather than rerun it accidentally

The checkpoint retains stage2/stage3 scientific files plus stage4 code, new
sequences and coordinates, controls, failed pilots, trajectories, analyses,
reference metadata and manifests. The original uploaded inputs and regenerable
software/indexes are excluded. Keep the original archives as well:
`resources.zip`, `egfr_inputs.zip`, `egfr_refinement_inputs.zip`.
Their SHA-256 hashes are in `reference/recovery.json` and the packaging manifests.

Restore the checkpoint into `/mnt/data`; it contains `egfr_campaign/...`.
Extract resources.zip to `egfr_campaign/inputs/provided/`, egfr_inputs.zip to
`egfr_campaign/inputs/transfer1/`, and egfr_refinement_inputs.zip to
`egfr_campaign/inputs/transfer2/`. Inspect member paths and sizes before extraction;
do not overwrite unrelated work. The two transfer zips already contain their
`legacy_fetch/` and `egfr_refinement_fetch/` subdirectories.

The OpenMM8.4 CPython3.13 Linux wheel is in the second transfer. In the current
campaign it was unpacked (not fetched online) to `stage2/runtime/python/`.
The original environment also has NumPy,SciPy,Biopython,RapidFuzz and Numba. The
NumPy recovery wheel is in transfer2. Do not install a different forcefield.

Use these environment paths in a fresh isolated restoration:

```bash
R=/mnt/data/egfr_campaign
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
export PYTHONPATH=$R/stage4/code:$R/stage3/code:$R/stage2/code:$R/code:$R/stage2/runtime/python
```

The BLAST binaries and indexes can be rebuilt by `stage4/code/blast_campaign.py`
from the archived transfer1 BLASTtarball and original compressed sequences.
This also rebuilds the augmented antibody FASTA/source map from transfer2.
`novelty_screen.py` rebuilds the CDR3 reference normalization. Hashes of all
excluded regenerated files are retained in the excluded-file manifest.

## Scientific execution dependency order

Do not run this whole sequence over the preserved results: several scripts cache
or refuse existing results and others are intentionally one-time provenance
patches. Reproduce in a separate clean campaign copy/container, preserving
stage3 inputs and using an empty stage4 intermediate/output directory.

1. `inspect_loops.py`, then `stabilize_design.py`: generate the first61 variants.
   `evaluate_designs.py` audits format, paired cysteines and starting geometry.
2. `refine_candidate.py <id> [human6ARU|mouseAF]`: joint450-iteration refinement and
   explicit microstate calculations. Initial IDs are in human_selection.json.
   `launch_tasks.py <config> <max_parallel>` can execute preserved selections.
3. `matched_control.py C00003`: same additional refinement with no mutation.
   This is a comparator, not an additional sequence.
4. `dynamics.py <id> 198002 4 2`: corrected short thermal trajectories. The
   source preparation and native-loop mask are in the current code. The three
   seed198001 pilot outputs are retained but scientifically excluded; do not
   regenerate them as accepted results. Native, oldC03/C09, S17/S22/S28 and
   MATCHED_C00003 controls are explicitly distinguished.
5. Focused mouse/acid-state refinement selections are in
   focused_refinement_selection.json. `sensitivity.py` uses the indicated roots
   for carboxylate and terminal-cap checks. See sensitivity_selection.json.
6. `dynamics_extended.py <id>` performs2pswarmup+20psproduction, seed198003.
   Extended native, matchedparent and S28 were run first. S22 was subsequently
   predeclared in final_proline_predeclaration.json and run by
   final_proline_comparison.py. No seeds were selected based on better outcomes.
7. `probe_histidine_opportunities.py`, `add_histidine_candidates.py`: create
   S62/S63 F107H probes, then refine them. `histidine_matched_controls.py
   S00062|S00063` gives each an unmutated comparator with equal refinement rounds.
8. Run `audit_refined.py`, `audit_final_context.py`, `surface_comparison.py`,
   `analyze_trajectories.py` and `trajectory_contacts.py`. The minimized-geometry
   omega/cystine thresholds are not automatic exclusions on thermal snapshots.
9. Run CDRnovelty and BLASTmain/followup scripts. The saved original run searched
   the initial61 plus a two-sequence followup with one exact-reference control
   in each database. `compile_novelty.py` verifies all63 query IDs in all3sets.
10. `test_stage4.py`, `test_integrator.py`, `proton_linkage_bound.py`,
    `summarize_partial.py`, `compile_results.py`, and `write_report.py` generate
    checks and the final stage summary. Checkpoint packaging is only allowed
    after all workers exit; independent SHA-256verification follows packaging.

## Reproducibility limits and provenance

Mutation enumeration is deterministic. Stochastic components explicitly use
MT19937 with saved seeds. OpenMM CPU floating-point reductions and native
library/hardware differences can prevent bitwise-identical physical trajectories.
The archived results and exact hashes are the definitive records of this run.

One-time patch scripts and superseded snapshots preserve the actual development
history; they are not extra scientific steps to rerun on corrected source. The
report lists caught errors and excludes pilots without deleting them. Every
candidate records its parent sequence file, parent refined geometry, hashes,
mutation and rotamer provenance. Scores invalidated by mutation are not reused
for ranking. Controls are kept in separate directories, not added to FASTA counts.
