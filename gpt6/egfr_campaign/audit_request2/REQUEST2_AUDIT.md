# Request 2 dependency and capability audit

Date: 2026-09-29. This revises Phone-a-Friend 2; it does not initiate request 3.

## Decision

Keep request 2, but supersede the original downloader with version 2.0.
The revised transfer, together with the first transfer and the supplied files,
contains all currently identified external dependencies of the selected CPU-only
classical workflow. This is not an end-to-end execution guarantee, a binding
claim, or an organizer eligibility decision. One request remains unspent.

## Verified locally

- Python 3.13.5, glibc 2.41, NumPy 2.3.5, SciPy 1.17.0, Biopython 1.86,
  Numba 0.65.1, llvmlite 0.47.0, pip 25.1.1 and packaging 25.0 are installed.
- libstdc++, libgomp, libgcc_s and libOpenCL are discoverable.
- OpenMM is not installed; no OpenMM energy calculation was performed.
- DNS resolution of PyPI, files.pythonhosted.org and OPIG still fails.
- The original checkpoint was safely restored. It records prior BLAST database
  construction. Large databases/indexes are not present in the checkpoint and
  must be re-extracted/rebuilt from the already uploaded first transfer when
  needed. No additional user transfer is needed for that restoration.
- The revised downloader passes 29 offline tests with fabricated HTTP responses.
  The tests do not establish live download availability or binary compatibility.

## Version 2 additions and hardening

1. Thera-SAbDab sequence-only CSV, in addition to PLAbDab paired, unpaired and
   nanobody sequence datasets. The organizer's novelty description explicitly
   includes a therapeutic-antibody database.
2. Pinned NumPy 2.3.5 Linux CPython 3.13 wheel for offline dependency recovery.
   Installed NumPy will not be replaced without need; installation will be
   campaign-local rather than modifying shared packages.
3. Standalone copies of amber99sb.xml, amber99_obc.xml, hydrogens.xml and
   residues.xml from OpenMM's 8.4.0 source tag. These support parameter auditing
   and a transparent fallback implementation, but a fallback is not yet written
   or validated. Loading raw XML here does not execute XML script elements.
4. Wheel checks for version, CPU plugin, core library, force-field/topology
   contents and uncovered base dependency names. glibc and CPython tags are
   checked for the remote runtime, not the downloading computer. Dependency
   version constraints and native linking still require checking on arrival.
5. Reuses validated cached files. Records URLs, SHA-256 hashes, byte counts,
   success/failure, and remaining limitations. Does not install or execute
   downloaded software on the user's computer.

## Workflow and dependency boundary

- Backbone/side-chain construction, sequence search, simulated annealing,
  glycan steric checks, sequence-liability filters and diversity selection:
  custom classical implementations on the installed numerical stack. Some
  modules remain unwritten; no external design service is assumed.
- Molecular mechanics: requested OpenMM wheel, explicit Amber99SB/OBC settings.
  Test local installation, native linking, CPU or Reference contexts, finite
  energies/forces, and minimization of a historical control before bulk use.
  Do not invoke a generic built-in test that silently selects a later force
  field. No GPU or paid/license-gated design package is assumed.
- pH selectivity: explicit protonation-state ensembles and thermodynamic
  binding polynomials, including unbound references. OpenMM's automatic
  hydrogen assignment is not a pH-dependent binding calculation. This custom
  model, control checks and parameter sensitivity analysis remain pending.
- Protonation variants: explicitly account for both neutral histidine
  tautomers and the charged variant. Reconstruct missing atoms and terminal
  caps with internal-coordinate methods; do not silently introduce PDBFixer.
- Glycans: retain conservative steric exclusion. Do not claim full glycan
  energetics from a protein-only force field or assume unobserved glycans are
  absent. No new carbohydrate simulation package is required by this workflow.
- Novelty: local BLAST plus conservative CDR boundary extraction/alignment and
  edit-distance checks. Existing/new annotations can be used as input data;
  ANARCI, MMseqs2, Foldseek and modern structure predictors are not executed.
- Delivery: ranked FASTA of 100 unique candidate sequences, companion ranked
  CSV and a 20-sequence submission subset; record all ranking uncertainty.

## Explicit residual limits

- Imported binary execution has not been tested and CPU throughput is unknown.
  Native CPU failure may allow OpenMM's Reference platform; this is not
  guaranteed before actual testing. Preserve request 3 for a demonstrated
  blocker, not an already planned routine transfer.
- The sequence databases are not a complete public patent dump and do not
  reproduce the organizer's exact database snapshots or private decisions.
  Local novelty screening cannot guarantee upload acceptance.
- New loops, docking scores, minimized energies and protonation models cannot
  establish expression, folding, affinity, cross-species binding, or absence of
  binding at neutral pH. No scoring result is converted to absolute KD without
  an independently justified calibration.
- Amber99SB/OBC and the engine do not complete the sequence-design algorithm;
  development and testing remain substantive work, not a mere installation.
- Scientific algorithm dating must be audited separately from modern runtime
  permissions. The checkpoint records PCG64 in an earlier pilot; subsequent
  final stochastic production should use MT19937, and any affected production
  must be regenerated or excluded before an unqualified method certification.

## Primary sources checked in this audit

- https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
- https://www.adaptyvbio.com/blog/novelty
- https://opig.stats.ox.ac.uk/webapps/plabdab
- https://opig.stats.ox.ac.uk/webapps/plabdab-nano
- https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/
- https://pypi.org/project/OpenMM/8.4.0.post2/
- https://pypi.org/project/numpy/2.3.5/
- https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/setup.py
- https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/openmm/app/data/amber99sb.xml
- https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/openmm/app/data/amber99_obc.xml
- https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/openmm/app/data/hydrogens.xml
- https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/openmm/app/data/residues.xml
- https://docs.openmm.org/latest/userguide/application/01_getting_started.html
- https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html

Latest documentation was consulted for infrastructure compatibility, not as
permission to use every scientific method included in a modern software release.
