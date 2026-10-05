"""Record method scope and restoration instructions; do not rewrite old history."""
from pathlib import Path
import json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sources=[
 {'method':'CCD loop closure','year':2003,'doi':'10.1110/ps.0242703','use':'executed independent de novo loops; no antibody donor loop'},
 {'method':'Generic VHH framework','year':2008,'doi':'10.1074/jbc.M806889200','structure':'3EAK','use':'framework only in ancestry; native loops exclusively controls'},
 {'method':'Ramachandran geometry','year':1963,'use':'classical dihedral/steric heuristics, not a trained canonical classifier'},
 {'method':'Canonical antibody loop structures','year':1987,'doi':'10.1016/0022-2836(87)90412-8','use':'considered rationale only; canonical classes not assigned or used as template labels'},
 {'method':'Powell derivative-free optimization','year':1964,'doi':'10.1093/comjnl/7.2.155','use':'continuous carboxylate proposal torsions; objective removed for physical tests'},
 {'method':'Amber99SB','year':2006,'doi':'10.1002/prot.21123','use':'all-atom physical potential'},
 {'method':'OBC implicit solvation','year':2004,'doi':'10.1002/prot.20033','use':'fixed-coordinate interaction subtraction and minimization'},
 {'method':'Wyman proton linkage','year':1964,'doi':'10.1016/S0065-3233(08)60190-4','use':'binding-polynomial arithmetic; free pKas are assumptions, not inferred values'},
 {'method':'Atomic-detail continuum electrostatics precedent','year':1990,'doi':'10.1021/bi00496a010','use':'historical precedent only; its pKa-prediction procedure or accuracy not claimed'},
 {'method':'Velocity-Verlet precedent','year':1982,'doi':'10.1063/1.442716','use':'synchronized positions/velocities; more precise than earlier shorthand Verlet1967'},
 {'method':'RATTLE','year':1983,'doi':'10.1016/0021-9991(83)90014-1','use':'position and velocity constraints; primary article https://www.sciencedirect.com/science/article/pii/0021999183900141'},
 {'method':'Andersen thermostat','year':1980,'doi':'10.1063/1.439486','use':'stochastic collisions; custom synchronized implementation'},
 {'method':'MT19937','year':1998,'doi':'10.1145/272991.272995','use':'all new stochastic scientific steps'},
 {'method':'Kabsch alignment','year':1976,'doi':'10.1107/S0567739476001873','use':'independently checked rigid framework fit'},
 {'method':'Gapped BLAST','year':1997,'doi':'10.1093/nar/25.17.3389','use':'modern2.17 user-permitted implementation; composition adjustment disabled'},
 {'method':'BLOSUM62','year':1992,'doi':'10.1073/pnas.89.22.10915','use':'sequence-search substitution matrix'},
 {'method':'Levenshtein distance','year':1965,'english_translation_year':1966,'use':'exact edit identity; modern RapidFuzz implementation'},
 {'method':'Rigid glycan-link torsion grid','year_upper_bound':2010,'use':'classical bond-axis rotations and geometry checks only; no modern glycan predictor or free-energy claim'}]
(S/'reference/source_register.json').write_text(json.dumps({'historical_sources':sources,'current_policy_urls':json.loads((S/'reference/rules_audit.json').read_text())['sources'],'web_access_notes':'Publisher DOI redirects for Swope1982 and Bashford1990 were attempted but unavailable; bibliographic metadata and historical primary references retained. No inaccessible article is claimed to have been fully read in this stage. Existing provided publications were inspected in earlier stages.'},indent=2))
audit='''# Stage 5 historical-method audit

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
'''
(S/'METHODS_AGE_AUDIT.md').write_text(audit);Path('/mnt/data/egfr_stage5_methods_audit.md').write_text(audit)
state='''# Stage5 continuation state

Completed finite checkpoint. Two Phone-a-Friend requests used, one unused.
Final100 FASTA/opaque codenames/officialCSV/top20 remain pending. Internal276 FASTA
is explicitly NOT_FOR_SUBMISSION and includes failures; do not upload as ranked.

Read REPORT.md, METHODS_AGE_AUDIT.md, RESTORE.md, output/final_evidence_inventory.json
and output/decision_ledger.json. All finished jobs are recorded in postflight.

276 new unique119-125aa designs,139 uniqueCDR3s;13439backbone attempts,94accepted.
231initialpass,275localCDRpass,230both. N3081002 is excluded(71.43%CDRidentity).
All276 completed all3localBLASTdatasets with nine exact positive controls.
Failed first DNA/PDBcontrol and interrupted antibody-run exits remain preserved;
repair_pdb_control/resume_antibody_search provide final complete recovery proofs.

B00000 is unchanged C00003-derived matched control, not a new sequence.
24new neutral human models,19combinedgeometrypass;4acid human models inclcontrol;
7mouse inclcontrol,6pass;30of36totalstrictcontextpasses.
H00011(Y104E) and H00021(D54E,E103Q,Y104E) retain checked human+mousegeometry.
H00020 has humanpH gains but fails mousepeptide100-101 omega(-153.697deg): keepgate.
T3070400/T3080401 clash with fixed1IVO; do not promote despite pHscores.

Independent engineered-acid+His minima (kcal/mol conditional contrast):
B00000 .600896615;H00011 .697139445;H00020 .709198805;H00021 .665495891.
H00011:729microstates/39366independentprior scenarios;others243/13122.
No predicted/measuredpKas, no calibratedKD and no neutral-pH nondetection claim.
H21 mouseacid-off scenarios are not a challenge disqualification:mousepHswitchnotrequired.
4matched unbound localminimizations complete;noHvariant thermal runs yet.

8corrected dynamics runs:2pswarmup+20psproduction,2seeds percontrol/variant.
Two-seed mean loop displacement:B1.8723,N7 2.0515,N9 1.8666,native1.2830A.
No demonstratedstabilization;N9 additionally has12/18N352glycangridclashes vs6/18parent.
Glycangrid is unweighted geometry,notprobabilities. 150resolvedheavyatoms only.

Next scientific priorities:compare priorC00009 under same expanded scopes;
separate observed modelgains from unverifiedfunction;repairgeometry before promotion;
consider H11/H21 thermal/conformational sensitivity and broader distinct-sequence
validation before finalranking. Do not endlessly pad a family with trivialvariants.
No knownbinder seeds, no modern generator/predictor. PriorPCGpilot excluded,disclosed.
'''
(S/'STATE.md').write_text(state)
root_before=S/'reference/root_state_before_stage5.md'
if not root_before.exists():root_before.write_bytes((R/'STATE.md').read_bytes())
(R/'STATE.md').write_text('# EGFR campaign current state\n\nCurrent continuation point: stage5/STATE.md.\nRead stage5/REPORT.md, METHODS_AGE_AUDIT.md, RESTORE.md and output/final_evidence_inventory.json.\nThis is a completed checkpoint, not the final100 sequence submission.\nTwo Phone-a-Friend requests used, one remains.\nEarlier root state is preserved in stage5/reference/root_state_before_stage5.md.\n')
restore='''# Restore and reproduce stage5

## Restore existing data without regenerating it

Extract this checkpoint to /mnt/data; it contains egfr_campaign/.... The supplied
archives resources.zip,egfr_inputs.zip,egfr_refinement_inputs.zip are still needed
for immutable target/dataset/software inputs. Their exact hashes are retained in
stage5/reference/restoration.json and checkpoint_archive_inputs.json.
The archive does not duplicate input payloads, installed runtimes,BLAST indexes,
reference normalization or caches. All scientific queries/results remain.

stage5/code/restore.py performs path-checked extraction to inputs/provided,
inputs/transfer1 and inputs/transfer2 and unpacks the supplied OpenMM wheel to
stage2/runtime/python. Do not install an unrecorded new forcefield or generator.
Use OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1. Existing scripts
add their campaign-local module paths. The user allowed current implementations
of historical methods. Database rebuilding is in build_databases.py;CDR references
are regenerated by novelty_screen.py,not replaced with a modern classifier.

## Scientific dependency order (reproduction must use a separate clean copy)

1. Inspect parent and declare strategies/rules. build_return_loops.py and
   design_returns.py/design_returns_neutral.py generate R/N branches. All seeds,
   attempts, failed backbones and selected sequence/rotamer metadata are retained.
2. audit_select.py and novelty_screen.py establish the predeclaredhuman subset.
   refine.py gives450neutraliterations and an independent strictaudit;B00000 is
   unchanged control with the same extra refinement,not additionalcandidate.
3. dynamics_run.py uses the corrected inherited synchronizedintegrator. The exact
   eight jobs/seeds are in dynamics_batch.json and predeclaration. Both20ps runs
   for eachsequence/nativecontrol completed. No stage4biasedpilots are accepted.
4. contact_geometry.py identifiedE106;build_preserved_tail.py/tail_extension.py
   and design_preserved_tail.py preserve thatcontact. build_beta_returns.py is a
   failed finite search,not acceptedcanonical templates. tail_audit/followup retain
   strictgeometry/contextfailures rather than movingthresholds.
5. dual_nitrogen_design.py and continuous_carboxylate.py construct21Hvariants.
   Eight were selected before neutralMM. AllHIPrefine_acid.py and mouse followups
   compareB00000/H11/H20/H21. H20mousepeptide gate fails and remains flagged.
6. sensitivity_three_acids.py corrects the old2acid scope. chemistry_acid_scope.py
   includes4acids forH11(729states),3forH20/H21(243). Samehistoricalphysics;
   acid_histidine_six_sites.py onlyraises the explicit exhaustiveguard.
   independent_mixed_priors.py writes allscenarioarrays with per-site independent
   assumedpKas/tautomerfractions and recovers oldtiedpriors exactly.
7. sensitivity.py cap,unbound_minimize.py,independent_priors.py,context_audit.py,
   glycan_flexibility.py and target_reaudit.py provide distinct sensitivity checks.
   They must not be combined into a claim of oneexhaustivevalidation calculation.
8. BLASTmain/tail/chemistry scripts searched all276. OriginalDNAcontrol and
   interruptedantibodyrun are failures, preserved;repair_pdb_control.py and
   resume_antibody_search.py document successfulrecovery. compile_novelty.py
   verifies querysets andcontrols;do not infer itfrom oldnonzeroexitlogs alone.
9. analyze_dynamics.py,test_stage5.py,catalog_campaign.py,compile_evidence.py and
   write_report.py create finaldata/report. Reporting-onlycorrection scripts are
   history,not steps to blindlyrepeat on correctedsource. Same for failedversions.
10. postflight.py must pass with no campaignworkers before make_checkpoint.py.
    verify_checkpoint.py independently checks every archivedpayloadSHA256/size.

The decisive completed results are the archived files and hashes. Hardware and
floating-point changes can alter energyminimization/MD, despite recordedseeds.
Do not overwrite the checkpoint's observedresults with a silent rerun. Allnew
models remain computational hypotheses;final100selection and upload are pending.
'''
(S/'RESTORE.md').write_text(restore)
print('Method audit, source register and continuation/restoration state written')
if __name__=='__main__':pass
