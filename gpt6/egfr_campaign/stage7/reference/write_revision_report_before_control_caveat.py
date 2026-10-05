"""Write revision documents directly from completed computational evidence."""
from pathlib import Path
import json,datetime,shutil,collections
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';O=S/'output'
def fmt(x,n=3):return f'{x:.{n}f}'
def main():
 summary=json.loads((O/'revision_summary.json').read_text());panel=json.loads((O/'final_panel_private_v2.json').read_text())
 evidence=json.loads((S/'reference/evidence_compilation.json').read_text());by={r['candidate_id']:r for r in evidence['records']}
 thermal=json.loads((S/'reference/thermal_analysis.json').read_text());paired=json.loads((S/'reference/paired_protonation_comparison.json').read_text());post=json.loads((S/'reference/release_validation_v2.json').read_text())
 quench=[json.loads(p.read_text()) for p in sorted((S/'intermediate/endpoint_quench').glob('*_audit.json'))];assert len(quench)==12
 raw=json.loads((S/'reference/thermal_final_geometry.json').read_text());assert raw['completed_endpoints']==12
 now=datetime.datetime.now(datetime.timezone.utc).isoformat();mapping={r['source_internal_key'].split(':')[1]:r['public_name'] for r in panel if r['new_this_revision']}
 def label(cid):return 'Matched unchanged parent (M00000)' if cid=='M00000' else mapping.get(cid,cid)
 promoted=summary['promoted_by_complete_comparison']
 decision=(f"{len(promoted)} new sequence(s) met the complete predeclared descriptive promotion rule and were placed ahead of the retained panel. This is a prioritization decision, not proof of improved biological function." if promoted else "No new sequence met the complete descriptive promotion rule. The original top-20 ordering is therefore retained; the revised reserve pool contains additional computationally characterized candidates rather than an invented new lead.")
 lines=[f'''# EGFR campaign revision 2: evidence-led sequence refinement

Completed UTC: {now}. Track 3. Two Phone-a-Friend requests used, one unused.

## 1. Deliverables and decision

The revised file contains **100 unique, ranked candidate sequences**, {summary['length_min']}-{summary['length_max']} residues long. **{summary['new_sequences_in_pool']} sequences in the pool are new in this revision; {summary['new_sequences_in_top20']} are in the top 20.** The remaining retained sequences keep their original public identifiers and exact amino-acid sequences. File order defines rank; any numeric prefix in a retained identifier is an immutable name from the earlier release, not a revised rank.

{decision}

The files are `egfr_ranked_100_v2.fasta`, the corresponding `egfr_track3_top20_v2.csv`, a top-20 FASTA, and a clearly labeled 100-entry reserve CSV. Use at most 20 active Track-3 entries, not an additional 20 on top of an earlier submission. No submission was performed by this campaign. All experimental outcomes and official novelty/classification remain unverified.

## 2. Scientific change and ancestry

This round generated **22 new sequences plus one unchanged-sequence control** from our own untested J00008 design, all 125 residues. Its parent already traces to the independently generated MT19937/CCD binding loops, not a known experimental EGFR binder. The generic nanobody framework was preserved. The native 3EAK structure is a separate fold control; its original binding loops were not copied into a design.

The new hypothesis was that selected glycine replacements could reduce excessive loop flexibility. Geometry-compatible glycine-to-alanine substitutions have a pre-2011 experimental precedent [3]; that precedent supplies a rationale, not a transferred stabilization result. Serine alternatives test compatible small polar side chains. We examined six negative-phi positions, 35, 51, 99, 102, 111 and 113, across saved human-bound, mouse-bound and unbound coordinates. Positive-phi glycines were retained. Finite empirical rotamer enumeration evaluated the worst packing across these conformations. No modern protein generator, learned potential or structure predictor was run.

All proposed pH-contact positions (52, 53, 54, 103, 104 and 106) and framework amino acids remained unchanged. There were no new cysteines, prolines or N-X-S/T sequons. Preserving contact identities does not guarantee preservation of their spatial presentation.

## 3. Structural checks and cross-species scope

All **23 human models** and **10 mouse models** completed all-atom refinement. All **33** pass the independently checked minimized-structure gate and the checked receptor-context exclusions. These checks include stereochemistry, peptide geometry, disulfide strain, clashes, full-reference receptor context, the tested alternative human receptor conformations, and resolved glycan obstacles. The same numerical thresholds were retained. Glycan conformer overlap fractions are unweighted geometric stress-test summaries, not steric probabilities or accessibility measurements.

These are model-quality checks, not demonstrations that the nanobody folds or binds either species. Mouse calculations use the supplied predicted mouse structure as an allowed dataset; no mouse structure prediction was run. The force-field binding calculations retain the inherited target crop, residues 334-505. Full receptor/glycan context is assessed geometrically, not by a fully sampled full-receptor free-energy calculation.

## 4. Expanded protonation comparison

The unchanged control and the two principal variants received the same expanded eight-site analysis: two receptor histidines and six interface carboxylates, with explicit histidine and carboxylate protonation alternatives. There are 6,561 microstates and 354,294 combinations of assumed independent parameters per design. These combinations are sensitivity analysis, not independent experiments, confidence intervals or a probability distribution over success. Other ionization sites and conformational entropy are not exhaustively sampled.

| Record | Minimum conditional pH contrast, kcal/mol | Maximum held-out energy error, kcal/mol |
|---|---:|---:|''']
 for cid in ['M00000','M00018','M00021']:
  x=by[cid]['expanded_protonation'];lines.append(f"| {label(cid)} [{cid}] | {fmt(x['proton_model']['minimum_contrast_kcal'],6)} | {fmt(x['max_heldout_energy_error_kcal'],6)} |")
 lines.append('''
Positive contrast means an acid-on direction in this model. It is **not** an experimental binding free energy, dissociation constant, binding ratio or demonstration of neutral-pH nondetection. The matched unchanged control is the proper comparison; small differences from the prior J00008 score are additional-refinement effects, not mutation gains.

Both variants have slightly lower conservative contrast than the matched parent. Their interaction proxies also become more favorable at neutral pH, so stronger modeled attraction cannot be presented as a specificity improvement. The exact paired-parameter comparisons are:

| Record | Change in minimum contrast | Range of paired contrast changes, kcal/mol | Parameter combinations with larger contrast |
|---|---:|---:|---:|''')
 for r in paired['comparisons']:
  lines.append(f"| {label(r['candidate_id'])} | {fmt(r['difference_of_minima_kcal'],6)} | {fmt(r['paired_parameter_contrast_change_min_kcal'],6)} to {fmt(r['paired_parameter_contrast_change_max_kcal'],6)} | {r['number_of_parameter_settings_with_larger_contrast']:,} / {r['scenario_count']:,} |")
 lines.append('''
The intended gain, if supported, is improved loop preorganization while retaining an acid-on model direction, **not** a stronger pH switch. No candidate is marked as satisfying the biological no-detectable-binding requirement.

## 5. Thermal tests and the artificial-restraint confound

The primary comparison consists of two separately seeded 20-picosecond production trajectories for each of the unchanged parent, M00018, M00021 and native fold control. Every run has a 2-picosecond warmup, 300 K target temperature, a 2-femtosecond timestep, Amber99SB/OBC, synchronized velocity-Verlet/RATTLE integration and Andersen-style constraint-group collisions generated by MT19937. Warmup position restraints are removed in production.

The inherited protocol adds artificial chirality restraints. Glycine replacements create additional stereocenters and therefore additional restraints, which could falsely appear as improved rigidity. Before any mutant thermal runs, the protocol was amended so the **primary comparison has no added chirality restraints**. Standard force-field improper terms remain. The four earlier added-restraint parent/native runs are retained as controls and are not substituted for the primary matched parent.

All saved CA-trace RMSDs were independently reconstructed by framework alignment and checked against the simulation records. The complete, unrestrained-production results are:

| Record | Seed | Mean loop RMSD, A | Late-half loop RMSD, A | Mean framework RMSD, A | Mean temperature, K |
|---|---:|---:|---:|---:|---:|''')
 primary=[r for r in thermal['replicas'] if not r['added_chirality_restraints']]
 for r in sorted(primary,key=lambda r:(r['candidate_id'],r['seed'])):
  lines.append(f"| {label(r['candidate_id'])} | {r['seed']} | {fmt(r['loop_mean_A'])} | {fmt(r['loop_late_mean_A'])} | {fmt(r['framework_mean_A'])} | {fmt(r['temperature_mean_K'],2)} |")
 lines.append('''
The predeclared descriptive rule requires at least 0.15 A lower mean loop displacement in **each** seed, no worsening in either late half, temperatures within 10 K of the matched parent, retained expanded pH contrast within 0.10 kcal/mol, and geometry/novelty safeguards. This rule is not a statistical significance test. Equal seed labels do not create identical atomic noise in different proteins.

| Variant | Mean loop reduction, seed 198005 / 198006, A | Late-half reduction, A | Descriptive thermal rule |
|---|---:|---:|---|''')
 for c in thermal['comparisons_without_extra_restraints']:
  ps=sorted(c['paired_seed_descriptions'],key=lambda z:z['seed']);lines.append(f"| {label(c['candidate_id'])} | {' / '.join(fmt(p['loop_mean_reduction_A']) for p in ps)} | {' / '.join(fmt(p['loop_late_reduction_A']) for p in ps)} | {'Pass' if c['descriptive_thermal_rule_pass'] else 'Not met'} |")
 lines.append('''
An additional descriptive check examines backbone movement at the six intended pH-contact positions. It was added during execution, before the complete comparison, and is not misrepresented as the original primary endpoint. Whole-loop RMSD can conceal worse movement at functional contacts.

| Variant | Contact-backbone reduction, seed 198005 / 198006, A |
|---|---:|''')
 for c in thermal['comparisons_without_extra_restraints']:
  ps=sorted(c['paired_seed_descriptions'],key=lambda z:z['seed']);lines.append(f"| {label(c['candidate_id'])} | {' / '.join(fmt(p['hotspot_backbone_mean_reduction_A']) for p in ps)} |")
 lines.append('''
These are exceptionally short, nonequilibrated stress tests, not folding simulations, melting temperatures, thermodynamic stability estimates or evidence of expression. Changed amino-acid masses and the limited trajectory duration also limit interpretation. CA motion is not a direct measurement of side-chain hydrogen-bond occupancy or preservation of the binding mode.

### Phase-matched endpoint geometry

Applying the minimized-structure 25-degree peptide-angle cutoff directly to single thermal endpoints rejected all four initial control endpoints, including both native fold controls. This exposed a phase mismatch in the proposed audit, not a reason to call the native protein unfolded. The raw outcomes were preserved. Before any complete unrestrained trajectory was available, we recorded a correction: quench **all twelve** endpoints identically with 500 minimizer iterations and no added positional/chirality restraints, then apply the unchanged strict numerical geometry thresholds to those minimized endpoints. No loop-motion threshold or pH criterion was relaxed.
''')
 lines.append(f"Raw endpoint strict passes: **{raw['strict_pass_count']} / {raw['completed_endpoints']}**. After phase-matched quenching, basic geometry passes: **{sum(r['basic_geometry']['pass'] for r in quench)} / 12**; independent strict passes: **{sum(r['strict_geometry']['pass'] for r in quench)} / 12**. All failed endpoints remain in the archive. These checks concern endpoints only; only CA traces, not full all-atom intermediate trajectories, were retained.")
 lines.append('''
### The unbound-minimization warning remains

Before these thermal tests, local unbound minimization rearranged the two variants more than the parent:

| Record | Loop CA displacement on unbound minimization, A |
|---|---:|''')
 for cid in ['M00000','M00018','M00021']:
  u=by[cid]['unbound_minimization'];val=u.get('unbound_metrics',u)
  # The exact saved key is validated rather than inferred from an unrelated score.
  def locate(obj):
   if isinstance(obj,dict):
    if 'loop_CA_RMSD_A' in obj:return obj['loop_CA_RMSD_A']
    for k in ['displacement','metrics','unbound_displacement','alignment','relaxation_metrics']:
     if k in obj:
      x=locate(obj[k])
      if x is not None:return x
   return None
  displacement=locate(u)
  if displacement is None:
   # Saved free-binder audit stores its alignment under this explicit field.
   for value in u.values():
    if isinstance(value,dict) and 'loop_CA_RMSD_A' in value:displacement=value['loop_CA_RMSD_A'];break
  assert displacement is not None,(cid,list(u))
  lines.append(f"| {label(cid)} | {fmt(displacement,6)} |")
 lines.append('''
Those local-reorganization results are retained even when another diagnostic is favorable. They do not by themselves measure folding free energy, but they prevent describing every assessment as an improvement.

## 6. Novelty and competition rules

All 23 distinct query sequences, including the unchanged control, completed searches against local Swiss-Prot, PDB sequence and augmented antibody databases with exact-reference positive controls and reconciled query identifiers. All 22 new sequences are below the local 70% CDRH3 edit-identity safeguard; the maximum is 50%. The local CDRH3 reference set contains 203,818 segments. Search inputs, results, positive controls, source mappings and hashes are retained.

The current Adaptyv policy has a separate antibody novelty scale: a conserved framework is not automatically disqualifying, and CDRH3 similarity is central once a sequence is classified as antibody-like [2]. We did not run the modern organizer classifier, IMGT numbering implementation, learned predictor, MMseqs2 or Foldseek. Local CDR extraction/edit-distance and BLAST are safeguards, not a reproduction of Adaptyv's exact pipeline. Official classification, novelty level and acceptance remain pending.

Track 3 permits **at most 20 designs**. Upload format is a ranked CSV with `name,sequence,molecule_class`; all delivered entries are intended nanobodies. The original de novo binding-loop lineage is preserved, but organizer eligibility is not self-certified. No public identifier contains an instruction, internal design label or methodological claim. Participant registration, rights, licensing and other personal attestations remain the submitter's responsibility.

The live challenge deadline has been extended to **October 6, 2026 at 23:59 AoE (UTC-12), which is October 7 at 06:59 in America/Chicago** [1]. Historical reports retain their earlier deadline observations; this revision corrects the current instruction. Human binding at pH6.5 with no detectable binding at pH7.4 remains the central criterion, with mouse cross-reactivity and human affinity also required. Mouse pH switching is not imposed by this campaign as an extra success criterion.

## 7. Selection, testing and provenance
''')
 lines.append(f"The v2 pool contains **{summary['unique_CDR3']} distinct designed CDRH3 sequences**, with at most **{summary['max_CDR3_repetition']}** entries sharing a CDRH3. This is a count, not a claim that all structural hypotheses are independent. Evidence tiers in the 100 are: {summary['evidence_counts']}. Tier 0 is expanded protonation plus both-species geometry, tier 1 narrower protonation plus both-species geometry, tier 2 detailed human evidence with more limited mouse checks, and tier 3 exploratory reserves. No tier certifies biological function.")
 lines.append(f"\nIndependent release validation passed **{post['passed']} checks**, including **{post['coordinate_sequence_checks']} coordinate/sequence identity checks**. Earlier in this revision, 369 generation/lineage/sequence-search checks, 104 finite-polynomial tests and 13 formatter/error-path tests passed. These counts describe implementation and consistency checks, not independent biological validation.")
 lines.append('''
The first parent-refinement attempt exceeded its command timeout before CPU thread limits were set explicitly; the inputs, source version and timeout record are preserved; no completed audit or protonation result was produced by that attempt, and a completed rerun is separately recorded. Source-review corrections (including a draft tautological framework assertion, fixed before execution) and a NumPy-boolean JSON-serialization failure are also retained. No failed computation has been silently relabeled a pass.

All 7,737 historical stage-6 checkpoint payloads were verified before continuation and reverified unchanged. New source files, design records, models, trajectories, energy tensors, controls, plans, exclusions, logs and hashes are isolated in stage7. The new checkpoint includes the previous archived scientific record plus this revision; regenerable runtimes, large reference FASTAs and indexes are omitted with paths, sizes, hashes and restoration instructions. Original uploaded archives remain external reproducibility inputs.

The optional `make_submission_v2.py` only formats the next eligible-by-name reserves after the user supplies rejection identifiers. It neither runs nor bypasses official novelty assessment. Keep the revised set to 20 total active entries. Retained identifiers always map to the same sequence.

## 8. Historical-method certification and remaining uncertainty

For this revision and the accepted lineage of the delivered v2 sequences, the executed scientific method classes predate 2011, using the expressly permitted modern implementations and datasets. No modern protein-binder generator or learned folding/interaction predictor was run. The early, excluded PCG64 pilot in the inherited campaign remains disclosed; it is not part of the accepted sequence ancestry. Therefore the certificate concerns the delivered lineage, not a false assertion that every exploratory action in the entire historical campaign met the cutoff.

No experiment has established any candidate's acidic binding, neutral-pH nondetection, mouse binding, affinity, folding, oligomeric state, expression yield, therapeutic safety or efficacy. The conditional models still allow undesired neutral binding. This revision supplies better-characterized experimental hypotheses and an explicit evidence-based priority order, not a certificate that all competition biological criteria are satisfied.

## Sources

[1] Live challenge and FAQ, checked October4,2026: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr

[2] Adaptyv novelty definitions, checked October4,2026: https://www.adaptyvbio.com/blog/novelty

[3] Matthews BW, Nicholson H, Becktel WJ. Enhanced protein thermostability from site-directed mutations that decrease the entropy of unfolding. PNAS84:6663-6667 (1987). https://doi.org/10.1073/pnas.84.19.6663

[4] Shifman JM, Mayo SL. Exploring the origins of binding specificity through the computational redesign of calmodulin. PNAS100:13274-13279 (2003). https://doi.org/10.1073/pnas.2234277100

Prior force-field, sampling, loop-construction, proton-linkage and cluster-expansion sources are preserved in the stage6 certificate and source registers. Numerical tables above are actual campaign computations, not values taken from these publications.
''')
 (O/'REVISION_REPORT.md').write_text('\n'.join(lines))
 certificate='''# Historical-method certificate: EGFR v2 panel

## Accepted-lineage certification

For the 100 delivered v2 sequences and their accepted design/evaluation lineage, the scientific method classes run or implemented are established in 2010 or earlier. Modern programming languages, implementations and datasets were used under the user's express allowance. No modern protein-binder generator, learned sequence model, learned structure predictor or learned interaction potential was run.

The revision uses geometry-compatible Gly-to-Ala/Ser substitution hypotheses, finite empirical rotamers, classical steric and multi-conformation packing, Amber99SB(2006), OBC generalized Born(2004), local energy minimization, classical explicit proton-state binding polynomials, a pre-2011 finite cluster-expansion strategy, Kabsch alignment(1976), velocity-Verlet/RATTLE(1967/1983), Andersen-style collisions(1980), MT19937(1998), BLAST/BLOSUM62(1997/1992) and Levenshtein edit distance(1960s). The sequence ancestry retains the prior de novo CCD loops and generic framework; original native binding loops were not design seeds.

The new glycine-replacement rationale is supported by Matthews, Nicholson and Becktel, PNAS(1987), DOI10.1073/pnas.84.19.6663. This is not evidence that our substitutions improve experimental stability. The primary new thermal comparisons omit the inherited additional chirality restraints to avoid artificial-stiffening bias; unchanged classical force-field improper terms remain. Matched endpoint quenching uses the same pre-2011 physical framework and does not establish folding free energies.

## Explicit historical exception

An early, discarded campaign pilot used modern PCG64 random sampling. That exception remains preserved and was superseded by regenerated MT19937 construction. Delivered sequences trace to the accepted regenerated lineage. It would be false to certify that every exploratory operation in the entire history complied; the affirmative statement is scoped to the delivered design/evaluation ancestry.

## Infrastructure and scope

Current ChatGPT assisted reasoning, code authoring and reporting. Protein sequences were emitted by the recorded classical construction and mutation-search scripts, not a modern protein generator. Current Python, NumPy, SciPy, Biopython, OpenMM, RapidFuzz, BLAST and other numerical/file utilities are implementations, not claims that those exact releases existed in 2010. The supplied mouse predicted structure and modern sequence databases were used as allowed data; no prediction was run to create them.

This certificate does not certify experimental binding, pH switching, folding, expression, safety, official novelty/classification, selection for synthesis or competition eligibility. Adaptyv's external submission assessment is separate from the restricted local workflow.
'''
 (O/'METHODS_CERTIFICATE_V2.md').write_text(certificate)
 readme=f'''# EGFR revised candidate files (v2)

100 unique ranked hypotheses, {summary['length_min']}-{summary['length_max']} residues. {summary['new_sequences_in_pool']} new sequences in this pool; {summary['new_sequences_in_top20']} new sequences in its top20. None is experimentally validated.

- `egfr_ranked_100_v2.fasta`: complete reserve pool, ranked by file order.
- `egfr_track3_top20_v2.csv`: the at-most20 Track3 submission set, exactly name,sequence,molecule_class.
- `egfr_track3_top20_v2.fasta`: same20 for inspection.
- `egfr_100_v2_RESERVES_NOT_SINGLE_UPLOAD.csv`: reserves, NOT an instruction to upload100 to Track3.
- `make_submission_v2.py`: optional offline replacement formatter.

Use v2 as a replacement set, not an additional20 on top of an existing submission. Retained public names always mean the same sequence. Numeric prefixes in old identifiers are historical names; current rank is file order. Official novelty/classification and all binding/folding/expression outcomes remain unverified.

Deadline: October6,2026,23:59AoE(UTC-12), corresponding to October7 at06:59Chicago. Check the portal for later organizer changes.

After Adaptyv rejects identifiers, place one exact identifier per line in `rejected.txt`, then run:

```bash
python make_submission_v2.py --fasta egfr_ranked_100_v2.fasta --reject-file rejected.txt --output egfr_track3_v2_revised_top20.csv
```

This utility does not assess or bypass novelty. Unknown identifiers and insufficient reserves raise an error. Its synthetic formatting tests are not protein-design candidates.

The private report, evidence codebook and complete scientific checkpoint are separate from the public submission-file bundle. Sharing methodological metadata is optional and may disclose the design approach; no metadata is submitted automatically.
'''
 (O/'README.md').write_text(readme);shutil.copy2(S/'code/make_submission_v2.py',O/'make_submission_v2.py')
 (R/'CURRENT_RELEASE.md').write_text('# Current EGFR release: v2\n\nRead stage7/output/README.md and REVISION_REPORT.md. The earlier STATE.md and stage6 release remain unchanged historical records. Current public files are stage7/output/egfr_ranked_100_v2.fasta and egfr_track3_top20_v2.csv. All biological outcomes and official novelty/classification are unverified. Two Phone-a-Friend requests used, one unused.\n')
 print('Wrote v2 report, certificate, README and current-release pointer.')
if __name__=='__main__':main()
