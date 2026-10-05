# Restore and inspect the final campaign checkpoint

Extract `egfr_campaign_final_checkpoint.zip` into `/mnt/data` on Linux. The public FASTA/CSV and private evidence are already present under `egfr_campaign/stage6/output`; inspection does not require running the scientific code.

Verify the archived payloads with Python's standard library:

```bash
python /mnt/data/egfr_campaign/stage6/code/restore_final.py
```

For numerical work, retain the original user-supplied `resources.zip`, `egfr_inputs.zip`, and `egfr_refinement_inputs.zip`. Their required SHA-256 hashes are recorded in `reference/final_checkpoint_inputs.json`. Restore their payloads and the supplied OpenMM wheel using:

```bash
python /mnt/data/egfr_campaign/stage6/code/restore_final.py --restore-inputs
```

The numerical environment was Linux x86_64, Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, Biopython 1.86, Numba 0.65.1, RapidFuzz 3.14.3, and the supplied OpenMM 8.4.0.post2 wheel. Modern implementations were permitted by the user. Do not substitute a modern learned potential or protein generator. Scientific scripts use `/mnt/data/egfr_campaign` paths. Set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, and OPENMM_CPU_THREADS=1. The observed run had eight CPU cores of quota and 4 GiB RAM.

BLAST indexes and normalized input FASTAs are regenerable with `stage6/code/build_databases6.py`; they and unpacked executables are not duplicated here. Final search queries/results, controls, hashes, and source database counts are retained. The supplied mouse structure is a modern input dataset; no structure prediction was run locally.

For reproduction use a separate clean working copy. The exact job plans, commands, seeds, source ancestry, and completed statuses are in `stage6/reference/`. Saved-result caches intentionally prevent silent overwriting. The observed outputs are the provenance record; a rerun can change floating-point minimization results and is not promised to be bit-for-bit identical.

The computational order for this final stage is: restoration and prior-hash checks; inventory and uniform old-model audit; matched J contact combinations and human/mouse calculations; K alternative-pose checks and L mouse completions; expanded cluster/protonation analysis and its independent tests; unbound local diagnostics with the documented postprocessing correction; local novelty query reconciliation; complete glycan screening; evidence-stratified selection; FASTA/CSV export; independent final postflight.

Read the earlier stage reports for ancestry before stage6. Failed/obsolete script variants are historical records, not an instruction to execute them against the final outputs. `restore_stage6.py` is the original pre-stage6 restoration script; use `restore_final.py` for this final snapshot because the root state file has intentionally changed.

The checkpoint archive has a complete included-payload manifest and a separate list of omitted regenerable stage6 files. The original input archives are external dependencies, not missing computed results. The independent archive verification is delivered alongside the ZIP. The excluded modern-RNG pilot is preserved and explicitly disclaimed in METHODS_CERTIFICATE.md.
