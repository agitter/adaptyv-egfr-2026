# Stage4 continuation state

Checkpoint complete; not the final submission. Two Phone-a-Friend requests used.

Primary machine-readable evidence: `output/final_evidence_inventory.json`.
Decision ledger: `output/decision_ledger.json`. Internal63-sequence FASTA is not
for submission. Read REPORT.md and METHODS_AGE_AUDIT.md before continuing.

63 new unique125aa sequences S00001-S00063. Stage3+4:1177records/1160unique;
all historical stage2-4 including rejected:1543records/1526unique. Do not inflate
current successful-candidate counts with archived failures or parent controls.
14 new detailed human models,9strict geometry pass;3mouse models,2pass.
All63 local CDR3and3-database BLAST searches complete, with exact controls.

C00003's matched additional refinement is in
`intermediate/matched_controls/intermediate/refined/`. This is the correct
baseline for initial stage4 mutations; not the old stage3 minimized geometry.
S22=N55P+N59P on C03;S28=Q57C+Y60C. Protected D54,Y100,E103 remain.
S28 shows no extended thermal advantage over matched parent.
The proline variant's 20-ps mean loop displacement is +0.184 A relative to the matched parent; the disulfide variant differs by +0.034 A. These are descriptive differences from one trajectory per sequence at that duration, not statistical significance or a stability free energy.

F107H probes S62/S63 have second-round matched unmutated controls under
`intermediate/histidine_matched_controls/S00062|S00063/`. Both reduce the
conservative contrast. Do not promote them as extra-proton successes.

11 corrected thermal runs:2pswarmup+4or20psproduction; seeds198002/198003.
No production target/position restraints; chirality restraints remain.
Three seed198001 half-step-velocity pilots are excluded and preserved.
New independent minimized-structure gate excludes S35,S36,S53,S59 (omega)
and S60 (cystine torsion). Thermal endpoint omega excursions are not the same
as minimized-structure failures, and occur in the native control too.

Actual neutral no-detection,affinity,mouse binding,folding/oxidation and official
novelty acceptance remain unverified. No final100or top20 submission artifact yet.
Continue with more preorganized loop backbones and better proton-only contacts,
not score inflation or endless histidine/cysteine substitutions. Current-source
patch scripts preserve development history and are not to be reapplied blindly.

Restore/runtime notes are in RESTORE.md. Large immutable inputs,OpenMM/BLAST
runtime and rebuildable indexes are excluded from the checkpoint but hashed and
recoverable from the original user archives. No worker is left running after
successful packaging. All code execution has been in this container.
