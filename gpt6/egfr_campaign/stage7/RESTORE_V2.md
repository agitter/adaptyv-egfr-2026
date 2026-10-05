# Restore, inspect, and reproduce the EGFR v2 campaign

Extract `egfr_campaign_v2_checkpoint.zip` into `/mnt/data`. Public files are in
`egfr_campaign/stage7/output`. The stage6 release remains unchanged historical
material. `CURRENT_RELEASE.md` points to v2; the old `STATE.md` was not rewritten.

Verify every manifested payload with only Python's standard library:

```bash
python /mnt/data/egfr_campaign/stage7/code/restore_v2.py
```

Original large archives remain external dependencies. Keep `resources.zip`,
`egfr_inputs.zip`, and `egfr_refinement_inputs.zip`; hashes and destinations are
recorded in stage6/reference/final_checkpoint_inputs.json and
stage7/reference/stage7_original_inputs.json. The v2 ZIP contains the earlier
archived scientific work, so the old checkpoint ZIP is not required to inspect
or restore it. To restore original data and the supplied OpenMM implementation:

```bash
python /mnt/data/egfr_campaign/stage7/code/restore_v2.py --restore-inputs
```

Scientific calculations used Linux x86_64, Python3.13, OpenMM8.4.0.post2 and the
recorded modern NumPy/SciPy/Biopython/Numba/RapidFuzz stack. Modern implementations
were expressly permitted, but scientific algorithms must remain pre-2011.
Set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and OPENMM_CPU_THREADS=1; the wrappers
also explicitly set the CPU platform thread default. The workspace quota was
eight CPU cores and 4GiB RAM. No GPU or learned protein model was used.

Regenerable reference FASTAs, BLAST indexes and unpacked binaries are omitted,
with hashes in checkpoint_v2_omitted_manifest.json. Stage7/code/build_databases.py
reconstructs the BLAST environment from the original supplied archive.
Stage7/code/novelty_screen.py regenerates the CDRH3 normalization from retained
references. Actual search queries, result tables, exact-reference controls,
positive-control logs and candidate novelty outcomes are included.

For numerical reproduction, use a clean working copy rather than overwriting
archived observations. Scientific scripts deliberately use absolute campaign
paths; caches and no-overwrite checks are part of the original workflow. The
step order is: lead inspection; glycine-replacement enumeration; matched human
refinement; mouse follow-up; expanded protonation; local unbound minimization;
thermal comparisons without added chirality forces; independent CA-trace
analysis; phase-matched endpoint quenches; novelty reconciliation; evidence
compilation; panel selection; release postflight; reporting and packaging.

The exact plans, seeds, sequence ancestry, mutation records, source hashes,
commands and completed task statuses are in stage7/reference and logs. Primary
thermal seeds are198005 and198006. The extra-restraint control trajectories and
all raw endpoint geometry failures remain preserved. The correction to apply
minimized-geometry thresholds after a matched quench was recorded before any
unrestrained trajectory completed.

Canonical current wrappers are in stage7/code. Draft, failed, or superseded
versions under stage7/reference are history, not instructions to run them.
`postflight_v2.py` checks file/sequence/coordinate consistency, not biology.
`make_submission_v2.py` is only a rank-preserving formatter for replacing names
rejected by Adaptyv. It never performs or bypasses organizer novelty screening.

Floating-point and native-library differences can change minimization and MD
results. Reproducing the algorithm is not promised to reproduce every coordinate
bit-for-bit. The saved outputs and SHA-256 manifests document the actual run.
All biological objectives and official novelty/classification remain unverified.
