"""Generate the checkpoint report from archived computed evidence."""
from pathlib import Path
import json
S=Path('/mnt/data/egfr_campaign/stage4');R=S.parent

def table(headers,rows):return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,row))+' |\n' for row in rows)
def main():
 f=json.loads((S/'output/final_evidence_inventory.json').read_text());dyn=f['dynamics'];h=f['human'];ctrl=f['controls']
 mdrows=[]
 for cid,label in [('NATIVE_3EAK','Native experimental framework control'),('MATCHED_C00003','Equal-refinement unmutated parent'),('S00028','Extra-disulfide variant'),('S00022','Two-proline variant')]:
  a=dyn[cid+'_v2_s198002'];b=dyn[cid+'_v2_s198003'];mdrows.append([label,f"{a['production_mean_loop_RMSD_A']:.3f}",f"{b['production_mean_loop_RMSD_A']:.3f}",f"{b['production_mean_temperature_K']:.2f}"])
 phrows=[]
 for label,cid in [('Equal-refinement parent','C00003'),('Two-proline variant','S00022'),('Extra-disulfide variant','S00028')]:
  m=ctrl['matched_C00003'] if cid=='C00003' else h[cid]
  phrows.append([label,f"{m['minimum_contrast_kcal']:.3f}",f"{m['representative_acid_proxy_kcal']:.2f}",f"{m['representative_neutral_proxy_kcal']:.2f}"])
 hirows=[]
 for probe,label in [('S00062','Uncrosslinked background'),('S00063','Disulfide background')]:
  a=ctrl['histidine_probe_parent_'+probe]['minimum_contrast_kcal'];b=h[probe]['minimum_contrast_kcal'];hirows.append([label,f'{a:.3f}',f'{b:.3f}',f'{b-a:+.3f}'])
 allrows=[]
 for cid,m in sorted(h.items()):
  d=json.loads((S/'intermediate/designs'/(cid+'.json')).read_text());edits='F107H' if cid in ['S00062','S00063'] else ','.join(e['from']+str(e['position'])+e['to'] for e in d['stabilization']['edits']);allrows.append([cid,edits,f"{m['minimum_contrast_kcal']:.3f}",'pass' if m['independent_geometry']['pass'] else 'exclude current geometry'])
 sensrows=[]
 for label,d in f['sensitivity'].items():sensrows.append([label,*[f"{d[k]['minimum_contrast_kcal']:.3f}" for k in ['native54','capped54','histidine_carboxylate162']]])
 parent=dyn['MATCHED_C00003_v2_s198003']['production_mean_loop_RMSD_A'];pro=dyn['S00022_v2_s198003']['production_mean_loop_RMSD_A'];ss=dyn['S00028_v2_s198003']['production_mean_loop_RMSD_A']
 delta=f'The proline variant\'s 20-ps mean loop displacement is {pro-parent:+.3f} A relative to the matched parent; the disulfide variant differs by {ss-parent:+.3f} A. These are descriptive differences from one trajectory per sequence at that duration, not statistical significance or a stability free energy.'
 text='''# EGFR design campaign: stage4

**Date:** October1,2026. **Status:** computational design and falsification
checkpoint, not the final submission package. **Phone-a-Friend:** two used, one
unused. All jobs must finish before archive packaging.

## Executive result

This stage generated63 new, distinct125-residue variants and subjected14 to
explicit human-interface molecular mechanics and protonation calculations.
Nine of those14 pass the stricter minimized-structure peptide/cystine gate.
Three have detailed mouse calculations; two pass the corresponding geometry
checks. No design has experimentally established human affinity, mouse binding,
or absence of neutral-pH binding.

The main scientific result is a more trustworthy comparison, not an assertion
that more mutations necessarily produced better binders. The apparent short-test
advantage of an extra disulfide did not persist in the longer matched-parent
comparison. An added histidine worsened the conservative pH contrast on both
backgrounds even after controlling additional refinement. A simple proline
variant was carried through the same extended test rather than being promoted
from its short trajectory alone.

The requested final ranked100 FASTA, opaque public codenames, competition CSV,
and top20 subset remain pending. The included63-sequence FASTA is deliberately
marked NOT_FOR_SUBMISSION and contains failed designs for provenance.

## 1. Criteria and scope

The live challenge/FAQ was rechecked. The primary biological objective is human
binding at pH6.5 with no detectable binding at pH7.4; mouse cross-reactivity and
human affinity also matter. No quantitative detection threshold is supplied on
the challenge page. A positive computed pH contrast is not that experimental
criterion. Track3 permits at most20 submissions; the100 requested sequences are
a ranked reserve pool, not authorization to upload100. Length10-250 and common
nanobody frameworks are allowed; existing binder seeds are prohibited. Our
variants derive from our own untested de novo loops, not an experimental binder.

Source: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr

The current stage3-plus-stage4 generations contain1177 sequence/structure records
and1160 distinct sequences. Counting the archived, mostly rejected stage2
exploratory generation as well gives1543 records and1526 distinct sequences.
Those broader historical counts are not counts of eligible or successful binders.
The63 sequences generated here have9 distinct CDR3s; they are a focused repair
branch, not a diverse final100 by themselves.

## 2. Historical design choices and actual generation

Starting from the refined C00003/C00009 own-design parents, I enumerated
geometry-compatible Gly-to-Ala and X-to-Pro substitutions, plus proposed extra
cystines. Finite rotamers, local sterics, mouse compatibility and fixed target
glycans were checked before accepting a proposal. Glycines with incompatible
positive phi were not casually replaced, and proline ring closure was checked.
The protected D54,Y100,E103 contact identities were retained. All63 output
sequences match their atom-residue records, contain only standard amino acids,
and lack N-X-S/T glycosylation sequons.

The stabilization rationale is pre2011: Matthews etal1987 demonstrated selected
entropy-reducing substitutions in lysozyme, while Saerens etal2008 showed that
extra nanobody disulfides are site-dependent rather than universally beneficial.
Neither paper establishes the effect of these particular modifications.
DOIs:10.1073/pnas.84.19.6663 and10.1016/j.jmb.2008.01.022.

Sixty-one initial variants were generated. Eight of the final63 sequences carry
an intended additional disulfide, including the later histidine probe. Modeled
bonds do not establish correct oxidative folding or expression yield.

Two further variants test F107H. A classical donor-direction/rotamer search found
only one plausible new histidine donor contact, to receptor E344. It found no
geometry supporting a two-nitrogen/two-distinct-carboxylate clamp. The possibility
that a neutral tautomer already satisfies the single contact was declared before
testing rather than ignored.

## 3. Stricter geometry catches false positives

The new independent gate checks signed alpha chirality, peptide omega, actual
modeled disulfides, sulfur bond lengths, CB-S-S angles and cystine chi3.
Four of the14 deeply evaluated variants fail a peptide omega threshold; S00060
fails the cystine torsion threshold despite an acceptable sulfur-sulfur distance.
Its chi3 is about155degrees, far from the intended near90degree geometry.

The original C00003/C00009 structures also fail the stricter25degree omega
threshold at the35-36 peptide bond. The old stage3 audit used a looser30degree
warning criterion, so its absence of that flag was not proof of ideal geometry.
An unmutated parent given the same extra450-iteration refinement as the new
variants passes the stricter gate. This control prevents attributing a benefit
of extra minimization to a sequence change.

'''
 text+=table(['New candidate','Changes on own-design background','Minimum pH contrast, kcal/mol','Minimized geometry'],allrows)
 text+='''
A separate final-coordinate audit tested23 refined complexes/control structures
against full target context, both fixed alternate human conformations1IVO/1NQL,
and the resolved glycan obstacles. All23 pass the2A heavy-atom overlap check in
those comparisons. This does not rescind peptide/cystine failures, prove binding
to every receptor state, or model unobserved glycan conformations.

The geometrically calculated interface areas for matched parent,S22,S28 are
approximately486,493,461square Angstroms, respectively, using the same
Shrake-Rupley sphere sampling. Surface area is not an affinity estimate.

## 4. Matched thermal stress tests, including a falsified apparent gain

The corrected protocol uses Amber99SB2006/OBC2004, explicit MT19937,
synchronized velocity-Verlet/RATTLE, and Andersen-style constraint-group
velocity collisions at300K. There are2ps of positional-restraint warmup, then
4or20ps production without the target or positional restraints. Chirality
restraints remain, which is an important artificial constraint.

Three earlier pilot trajectories refreshed half-step leapfrog velocities without
proper synchronization and produced a309-312K bias. They are retained but
excluded. The corrected integrator passed independent harmonic-oscillator and
position/velocity-constraint checks. The native control uses native loops only
for calibration, never as a source of designed binding loops.

The table compares average loop C-alpha displacement after framework alignment.
The4ps and20ps columns use different predeclared seeds; each duration is compared
across matched sequences. They are not independent-frame statistical samples.

'''
 text+=table(['Structure','4ps production mean, A','20ps production mean, A','20ps mean temperature,K'],mdrows)+'\n'+delta+'\n'
 text+='''
The extra-disulfide design appeared better in the4ps comparison but was not better
than the equally refined parent over20ps. It must not be called stabilized on
that evidence. Loop-internal deformation is separately fitted to the initial
recorded unbound frame for every structure, avoiding the earlier mismatch between
native and designed reference choices. Native loops remain the stronger control
for internal geometric preservation.

There are11 corrected trajectories,108ps total production and130ps aggregate
simulation including warmup. This is not one130ps trajectory and certainly not
an equilibrium folding study. No folding free energy, melting temperature,
kinetic rate or long-timescale stability is inferred. Thermal endpoints may
cross the minimized-structure omega/cystine thresholds even for the native
control; those instantaneous excursions are diagnostics, not automatic sequence
exclusions. All corrected endpoints pass the separate bond/severe-clash/chirality
checks used for trajectory integrity.

Framework-aligned endpoint contact audits are also preserved. Receptor geometry
is held fixed in those overlays; displaced free loops can have new overlaps.
No initial carboxylate-His contact was below the prechosen4A distance threshold
in that audit, so the code reports zero eligible contacts rather than moving the
threshold to manufacture a favorable retention percentage. Continuous distances
are retained. These overlays are not induced-fit or binding simulations.

## 5. pH contrast survives, but neutral rejection remains unproved

All14 new detailed human models retain an acid-on direction across54 salt,
dielectric,pKa and neutral-tautomer assumptions. This is a conditional model
sensitivity exercise, not54 independent validations. Representative absolute
interaction proxies remain favorable at neutral pH for the main family.

'''
 text+=table(['Matched model','Worst contrast across54','Representative acid proxy','Representative neutral proxy'],phrows)
 text+='''
All energy entries are kcal/mol. The absolute interaction proxies omit key
binding entropy and ensemble terms and are not calibrated binding free energies.
They must not be exponentiated into claimed experimental dissociation constants.
Positive contrast alone does not establish no detectable neutral binding.

Explicit all-protonated interface refinement produces a second structural
hypothesis for S17,S22,S28. The corresponding minimum contrasts are about
1.116,1.129,1.131kcal/mol; the more favorable structural hypothesis is not
cherry-picked as the final score. The lower contrast over hypotheses remains
relevant. Mouse calculations are favorable in interaction-energy terms for all
three, but S17's mouse geometry fails the stricter peptide gate. S22 and S28
pass it. Mouse pH switching is not required and was not imposed as a penalty.

Additional carboxylate protonation and artificial target-terminus capping checks
were completed for matched parent and S28:

'''
 text+=table(['Background','Native54 minimum','Capped54 minimum','His+carboxylate162 minimum'],sensrows)
 text+='''
The F107H probes were compared not just with their old parents but with parents
given exactly the additional refinement round used for each probe:

'''
 text+=table(['Background','Matched parent minimum','F107H minimum','Difference'],hirows)
 text+='''
Neither is promoted as a pH improvement. This is stronger evidence than comparing
against an unrefined or differently refined parent, but still conditional on
this physical model and sampled geometries.

A separate classical proton-linkage calculation clarifies the design challenge.
If binding changes proton uptake by at most two across pH6.5to7.4, the maximum
association-constant ratio is10^(2*0.9), about63.1, and the maximum contrast at
300K is2.471kcal/mol. This is a conditional thermodynamic bound, not a universal
limit for proteins and not a predicted affinity ratio for these designs.
Additional ionizable groups or coupled conformations can change the assumptions.
There is no published detection threshold here from which to deduce that63-fold
is either necessary or sufficient. The1,200 algebraic property checks validate
the code's bound, not a biological outcome. Classical linkage references are
Wyman1964 and Baker/Murphy1996, DOIs10.1016/S0065-3233(08)60190-4 and
10.1016/S0006-3495(96)79403-1.

## 6. Novelty screening is complete locally for the63 new sequences

Every query was searched against rebuilt Swiss-Prot575748records,
PDB1165667records, and augmented antibody446966records. These total2188381
indexed records, not distinct proteins across databases; the supplied PDB
sequence archive includes polymer records beyond proteins. All three exact
reference positive controls passed. Query-ID parsing verifies63/63 searches
in every database, rather than inferring success from a log message.

Among reported hits covering at least90percent of a query, maximum identities
are61.6percent in Swiss-Prot,72.0percent in PDB,70.4percent in the antibody set.
These are top20 BLAST-hit statistics, not exhaustive global alignments. Shared
antibody framework similarity is expected and is evaluated separately from CDR3.

Exact edit-distance screening against203818 reconstructed/provided CDR3 segments
finds all63 below70percent, with maximum50percent nearest-reference identity.
Nine distinct CDR3queries account for the63 focused variants. The reference
normalization matches the previous-stage count and its hash is recorded.

The organizers' modern classifier, numbering and exact database snapshots were
not replicated. Generic-framework novelty and antibody-CDR novelty must not be
confused. Local screening is not an upload-acceptance certificate.
Source:https://www.adaptyvbio.com/blog/novelty

## 7. Tests, provenance and methodological limits

The207 sequence/geometry/reproducibility unit checks pass. Three dedicated
integrator checks, the rigid-alignment check,1,200 binding-polynomial property
checks and three exact BLASTpositive controls also pass. These implementation
checks do not establish model accuracy or experimental function.

All new candidate records retain parent sequence and refined-coordinate hashes,
mutation/rotamer provenance and explicit invalidation of inherited scores.
Unmutated controls are separate records outside the candidate FASTA. Failed
proposals, original temperature-biased pilots, stricter geometry rejections,
one-time code corrections and native controls remain in the archive. Earlier
stage3 checkpoint payload hashes were independently verified on restoration.

The algorithm-age audit is in METHODS_AGE_AUDIT.md. No modern learned protein
generator,predictor or potential was run. Current datasets and implementations
are within the user's allowance. The excluded early PCG64 campaign pilot remains
disclosed; present candidates trace to the MT19937-regenerated lineage. A blanket
claim that literally every historical exploratory operation complied would be
false. This checkpoint is not the final FASTA certification.

## 8. Design decision and continuation

Do not promote any variant because it contains prolines, an extra disulfide or
an extra histidine. The current evidence ledger separates geometry-clean
hypotheses, explicit failures, untested properties and nonpromotion decisions.
The repaired parent and simple proline/disulfide variants remain hypotheses,
not measured binders; neutral-pH rejection is still the central unsolved design
requirement in this campaign.

The next useful search should improve loop backbone preorganization and
proton-dependent contact geometry rather than produce another large collection
of near-identical substitutions. Pre2011 canonical-loop constraints and shorter
CCD-generated alternatives can be tested while preserving the intended receptor
contacts. Chothia/Lesk1987 is a historical basis to investigate, not a method
already implemented here. Additional proton coupling should require contacts
that neither neutral histidine tautomer can satisfy, then be evaluated in the
same full multistate model. Broader epitope/loop diversity is necessary before
selecting100 reserves and20 actual submissions.

See STATE.md,RESTORE.md and the machine-readable evidence inventory for precise
continuation paths. All current workers are checked for completion before
packaging. No third Phone-a-Friend request was consumed.
'''
 (S/'REPORT.md').write_text(text);Path('/mnt/data/egfr_stage4_report.md').write_text(text)
 state=f'''# Stage4 continuation state

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
{delta}

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
'''
 (S/'STATE.md').write_text(state);Path('/mnt/data/egfr_stage4_methods_audit.md').write_text((S/'METHODS_AGE_AUDIT.md').read_text());print('Wrote report',len(text.split()),'words')
if __name__=='__main__':main()
