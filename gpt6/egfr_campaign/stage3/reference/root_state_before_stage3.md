# EGFR campaign checkpoint - 2026-09-29

## Status

No final candidate sequences have been produced. No experimental binding, pH
selectivity, folding stability, or organizer novelty decision is claimed.
The current phase has produced a target ensemble, regenerated antibody-loop
backbones, and geometrically filtered docking poses. Sequence/rotamer design,
all-atom refinement, novelty screening of actual designs, and final ranking
are pending. Phone-a-Friend request 1 is complete. Request 2 is the accompanying
fetch_egfr_refinement_inputs.py transfer. One request remains after request 2.

## User contract

Deliver 100 ranked sequences in FASTA with cryptic codenames. Also prepare the
competition's CSV columns name,sequence,molecule_class, and a top-20 subset
because Track 3 permits at most 20 actual submissions. The 100 are a reserve
pool, not permission to submit 100. Scientific methods must date to 2010 or
earlier. The user expressly permits modern Python and modern implementations
of old tools, and modern input datasets. No modern protein-generative or
folding-model inference is permitted. Keep all code, inputs, intermediates,
failures, provenance, scoring assumptions, and output mapping.

## Competition priorities

Human EGFR at pH 6.5, no detectable human binding at pH 7.4, and mouse binding.
Selection prioritizes pH selectivity, then mouse cross-reactivity, then human
affinity. Sequence length 10-250. Common nanobody frameworks are expressly
permitted; modifying a known binder is not. Eligibility cannot be guaranteed
by a local partial database search. In particular, do not substitute a count
of histidines for a calculation of proton-linked binding.

## Strategy comparison and decision

1. Known EGF/cetuximab/EGFR-affibody/G532 variants: rejected as design seeds.
   Existing binders can be controls but cannot provide the new binding loops.
2. Unconstrained short peptides: poor fold/presentation confidence; not selected
   as the main approach.
3. Binary-patterned helical bundles and de novo alpha/beta designs: historically
   feasible but substantial sequence/fold novelty-classification risk under
   Proteinbase's structural-similarity criteria. Retain as a conceptual option,
   not an executed branch.
4. Generic VHH framework with wholly rebuilt CDR1, CDR2 and CDR3: selected.
   The framework is from the published 2008/2009 generic humanized scaffold
   study, PDB 3EAK chain A. ALL original loop sequences and loop coordinates
   at source residues 26-38,53-67,100-117 were deleted. Framework coordinates
   are retained. Loops were sampled from non-antibody backbone-torsion
   fragments and rebuilt by CCD. No parent antibody CDR seed was used.
5. Binder histidines paired to target acids: selected for explicit testing,
   not assumed sufficient. Neutral histidines can also hydrogen bond, and
   burial/desolvation may oppose protonation.
6. Binder acidic groups paired to conserved receptor histidines H370/H433:
   selected for explicit testing. Receptor-based switching may avoid putting
   all protonation liability in the binder. H383 is not conserved (mouse R).
7. Combined receptor and binder titration sites: selected as a multistate
   hypothesis, with electrostatic/solvation uncertainty and neutral-pH
   negative design to be included. Avoid pH-dependent target-conformation
   confounding: compare the same structural ensemble at both pH values.

## Work actually completed

- Safe extraction and SHA-256 input manifest.
- BLAST+ 2.17.0 executables extracted/tested (blastp/makeblastdb/blastdbcmd).
- Swiss-Prot 2026_03 BLAST database: 575748 sequences.
- PDB seqres BLAST database: 1165667 sequences.
- Full human/mouse global affine-gap alignment and target coordinate mapping.
- Domain-III superposition of human structures and supplied mouse model.
- SASA, sequon and resolved-glycan analyses.
- Generic framework defined and all original CDRs removed.
- 11200 random loop starts; 129 pass closure, torsion and steric prefilters.
- 500 accepted assembled backbones; 120000 rigid-body orientations tested.
- 389 poses retained, involving 207 distinct assembled backbone trials.
- Empirical sidechain conformers extracted from non-antibody controls/targets.
- Motif-based provisional novelty set: 5094 unique putative heavy CDR3s.
  This is NOT full IMGT annotation or complete patent/antibody coverage.
- Offline downloader tests: 11 passed. No live download was verified here.

## Important target facts calculated from supplied/transferred data

Numbering below is full-length human UniProt P00533, not mature/PDB numbering.
Human H370/H433 are conserved in mouse. Their sidechain charge-center distance
in human 6ARU is about 9.8 A. Human domain III and aligned mouse have CA RMSD
about 0.499 A over the chosen alignment; this does not establish equivalent
binding energetics. Human N361 sequon is absent in mouse, a local glycan risk.
Target experimental glycans and potential sequons are included in geometric
exclusions. The recommended domain-III MM crop 334-505 retains disulfides
337-362 and 470-499; cap chain ends and also check complete ECR sterics.

## Planned next steps (not yet executed)

- Acquire OpenMM Linux CPython-3.13 wheel and PLAbDab paired/unpaired/nano
  sequence data through request 2. Do not download antibody model archives.
- Install OpenMM locally, benchmark CPU minimization. Use only Amber99SB
  (2006) and OBC GB (2004), not modern ML potentials or later force fields.
- Sequence/rotamer optimization by classical pairwise physical scoring and
  simulated annealing, including self-packing and positive/negative states.
- Histidine protonation microstates/binding polynomials; solvent and pKa
  sensitivity. Do not turn uncalibrated energies into absolute KD values.
- Refine/rerank with actual force-field calculations and multiple target
  conformations, glycan checks, sequence liabilities and diversity controls.
- BLAST and CDR edit-distance novelty searches including transferred patent
  and antibody sequence datasets. Annotated input data are allowed; do not
  execute ANARCI (2015), Foldseek, MMseqs2, or modern protein-design models.
- Emit validated 100-candidate FASTA and CSV plus top20 CSV. Certify only
  the historical-method constraint; biological success is unvalidated.

## Code order and execution notes

00_ingest.py; 01_prepare.py; 02_targets.py; 03_loops.py;
04_rotamers_and_cdr_db.py; 05_dock.py. legacy_geometry.py supplies loop geometry.
Failed versions/logs are retained. 02_targets_v1_failed.py and
05_dock_v1_failed.py are superseded, not successful production versions.
05_dock.py is importable without rerunning its main search, for assemble().
Loop-generation pilot used NumPy default_rng/PCG64 as numerical infrastructure;
docking used MT19937. Prefer MT19937 for further stochastic scientific runs.
No inference model, trained modern protein generator, or modern scoring model
has been run. NumPy/SciPy/BioPython/Numba are modern implementations/tooling.

The container cannot resolve outside download/package hosts; browser research
is separate. CPU only; Python3.13.5, Linux x86-64, glibc2.41. About28GiB free
at this checkpoint. Use container tools for execution. Interactive streaming
sessions are unavailable; batch processes were run with nohup and inspected
within the same working turn. No task is promised to run after a response.

## Checkpoint scope

The checkpoint ZIP includes generated code, small intermediates, raw structure
inputs, supplied documentation, logs, manifests and this state record. It
excludes large FASTA copies, BLAST indexes/binaries, downloaded installation
archives and caches because these are reproducible from the user's retained
resources.zip and egfr_inputs.zip. Exclusions and hashes are enumerated in the
checkpoint manifest. No font files are distributed.

## Request-2 audit update (2026-09-29)

Original request-2 downloader is superseded by fetch_egfr_refinement_inputs_v2.py
(version 2.0). It adds a NumPy 2.3.5 wheel, Thera-SAbDab CSV, standalone
Amber99SB/OBC/topology XML files, and stricter validation. 29 offline tests
passed; no live downloads or OpenMM execution verified. Detailed remaining
work and limitations: audit_request2/REQUEST2_AUDIT.md. Request 2 has not yet
been fulfilled. Request 3 remains unspent. The checkpoint has been restored
in this current runtime, but large BLAST indexes/binaries remain to be restored
from already uploaded egfr_inputs.zip; previous successful logs are historical.
