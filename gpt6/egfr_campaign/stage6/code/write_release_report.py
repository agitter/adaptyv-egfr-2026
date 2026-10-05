"""Write the release report and explicitly scoped historical-method certificate."""
from pathlib import Path
import json,collections,hashlib,shutil,datetime
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';O=S/'output'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 panel=read(O/'final_panel_private.json');v=read(O/'release_validation.json');models=read(S/'reference/new_model_summary.json');comb=read(S/'reference/combination_generation.json');newseq=sum(r.get('is_new_sequence',False) for r in comb['records']);gen=read(S/'reference/final_selection.json');allhist={};records=0
 for stage in ['stage2','stage3','stage4','stage5','stage6']:
  for p in (R/stage/'intermediate/designs').glob('*.json'):
   d=read(p)
   if isinstance(d,dict) and 'sequence' in d:
    allhist.setdefault(d['sequence'],[]).append(stage+':'+d.get('candidate_id',p.stem));records+=1
 clusters={p.stem.split('_cluster')[0]:read(p) for p in (S/'intermediate/cluster_acids').glob('*_cluster.json')};assert len(clusters)==11
 controls={'J00000':'Unchanged C00003 parent','J00001':'Unchanged C00009 sequence','J00002':'Unchanged C00017 sequence','J00003':'Unchanged H00011 sequence','J00004':'Unchanged H00021 sequence'}
 table=[]
 for cid in ['J00000','J00001','J00002','J00003','J00004','J00005','J00006','J00008','J00009','J00010','J00014']:
  d=clusters[cid];t=d['proton_model'];table.append(f'| {cid} | {controls.get(cid,"New contact combination")} | {t["minimum_contrast_kcal"]:.4f} | {len(d["sites"])} | {d["held_out_max_absolute_energy_error_kcal"]:.5f} |')
 top=[]
 for r in panel[:20]:top.append(f'| {r["rank"]} | {r["public_name"]} | {r["candidate_key"]} | {r["tier"]} | {r["ranking_contrast_kcal_proxy"]:.4f} |')
 cm=[m for m in models if m['key'].startswith('stage6:K') and m['species']=='human'];altpass=[m for m in cm if m['combined_geometry_pass'] and m['independent_priors']['all_acid_on']];human=[m for m in models if m['species']=='human'];mouse=[m for m in models if m['species']=='mouse'];unbound=[]
 for cid in ['J00000','J00001','J00002','J00003','J00004','J00005','J00008','J00009','J00010','J00014']:
  d=read(S/'intermediate/unbound'/f'{cid}_free_audit.json');a=read(S/'intermediate/unbound'/f'{cid}_free_independent.json');unbound.append({'candidate_id':cid,'loop_displacement_A':d['displacements']['loop_CA_RMSD_A'],'potential_relaxation_kcal':d['potential_relaxation_kcal'],'native_geometry_pass':d['geometry_audit']['pass'],'independent_geometry_pass':a['pass']})
 summary={'release_utc':v['utc'],'historical_sequence_records_including_duplicates_and_failures':records,'historical_unique_sequences_including_failures':len(allhist),'new_unique_sequences_this_stage':newseq,'stage6_sequence_structure_records':len(list((S/'intermediate/designs').glob('[JKL]*.json'))),'stage6_detailed_models':len(models),'stage6_human_models':len(human),'stage6_mouse_models':len(mouse),'stage6_uniform_geometry_pass':sum(m['combined_geometry_pass'] for m in models),'stage6_human_uniform_geometry_pass':sum(m['combined_geometry_pass'] for m in human),'stage6_mouse_uniform_geometry_pass':sum(m['combined_geometry_pass'] for m in mouse),'alternative_pose_human_pass_geometry_and_His_direction':[m['key'] for m in altpass],'cluster_models':len(clusters),'cluster_numerical_validation_pass':sum(d['numerical_validation_pass'] for d in clusters.values()),'cluster_max_held_out_energy_error_kcal':max(d['held_out_max_absolute_energy_error_kcal'] for d in clusters.values()),'cluster_total_exact_energy_evaluations':sum(d['exact_energy_evaluations'] for d in clusters.values()),'unbound_local_checks':unbound,'selected_panel':v,'candidates_in_otherwise_eligible_pool':gen['eligible_count'],'phone_a_friend_used':2,'phone_a_friend_unused':1}
 (O/'campaign_release_summary.json').write_text(json.dumps(summary,indent=2))
 tests=read(S/'reference/stage6_tests.json');ft=read(S/'reference/formatter_tests.json');lead=panel[0]
 report=f'''# EGFR classical-design campaign: ranked 100-sequence release

**Release date:** October 2, 2026. **Track:** 3. **Status:** computational designs for experimental testing; no experimental binding or expression data. **Phone-a-Friend:** two requests used, one unused.

## 1. Deliverables and decision

The release contains **100 distinct, ranked amino-acid sequences**, with public codenames and separate private evidence. Their lengths are **{v['length_range'][0]}–{v['length_range'][1]} residues**. All are intended single-domain nanobodies on a common VHH framework. The public upload files contain no internal design identifiers, instructions, or embedded prompts.

`egfr_ranked_100.fasta` is the requested complete reserve pool. `egfr_track3_top20.csv` is the ordered, 20-entry submission file. Its columns are exactly `name,sequence,molecule_class`, with `nanobody` as the class. `egfr_100_RESERVE_POOL_NOT_SINGLE_UPLOAD.csv` is a convenience export, **not a valid recommendation to upload all 100 in Track 3**. The top-20 FASTA is also supplied.

The ranking is an experimental-priority ordering, **not a prediction that any candidate wins or meets all three biological criteria**. Every candidate retains explicit UNVERIFIED fields for acidic human binding, neutral-pH nondetection, mouse binding, affinity, folding, and expression. Neutral-pH interaction-energy proxies remain favorable; no numerical score is relabeled as evidence of neutral-pH nondetection.

The first-ranked candidate is **{lead['public_name']}** (private source `{lead['candidate_key']}`). Its conservative ranking contrast is **{lead['ranking_contrast_kcal_proxy']:.4f} kcal/mol within the stated conditional model**, not a measured binding free energy or a calibrated dissociation constant.

## 2. Current rules and novelty policy

The live challenge and updated novelty article were reviewed on October 2, 2026 [1,2]. The challenge prioritizes human pH-selective binding, then mouse cross-reactivity, then human affinity. Human binding is required at pH 6.5, with no detectable binding at pH 7.4. Mouse binding is required, but mouse pH switching is not. The requested 100 is a reserve pool because Track 3 permits at most 20 submitted designs. The deadline is October 4, 2026, 23:59 AoE, corresponding to October 5 at 06:59 in America/Chicago.

The updated antibody-specific novelty scale uses CDRH3 edit identity. A value below 70% is the novel-CDRH3 branch; a familiar framework can lead to antibody novelty level 3 rather than a generic-protein rejection. We therefore **do not mutate conserved framework residues solely to lower whole-chain identity**, and do not apply the generic structural novelty scale to reject these intended nanobodies.

**Adaptyv's submission assessment is authoritative.** Our CDR boundary extraction is not the organizer's IMGT annotation, our BLAST searches are not its MMseqs2 pipeline, and our databases are not guaranteed to match every patent or current sequence snapshot. Neither official antibody classification nor upload acceptance has been established. No modern organizer classifier or structure predictor was run locally.

The generic framework is retained, but its original three binding loops were removed. Accepted loops trace to independently generated MT19937/CCD construction, with subsequent classical sequence and geometry design. Refining our own untested computational sequences is distinct from taking an existing experimental binder as the initial seed. The final lineage audit checks this recorded construction; it is not an organizer eligibility ruling. Participant registration, rights, licensing, and other personal eligibility attestations remain the submitter's responsibility.

## 3. What this stage actually added

Stage 6 generated **{newseq} new unique sequences** in a 21-member contact-combination experiment: 16 new combinations and five unchanged-sequence controls. Every member was constructed on the same saved C00003-derived backbone, given matched refinement, and evaluated against both human and mouse targets. This removes unequal refinement histories as an explanation for favorable comparisons. The unchanged C00009 and C00017 sequence controls repair the earlier peptide-geometry problem on this shared backbone; the earlier failed coordinate models remain archived.

Ten existing designs from alternative docking poses received matched human and mouse calculations. Six additional existing sequences received missing mouse follow-up calculations. These K/L records are reassessments, **not additional unique sequences**. This stage has **{len(models)} completed detailed models**, of which **{sum(m['combined_geometry_pass'] for m in models)} pass the combined geometry/context gate**. The underlying campaign contains **{len(allhist):,} unique historical sequences**, including rejected designs; this is not a successful-binder count.

For the ten alternative poses, **{len(altpass)}** currently pass both the checked human geometry and acid-on direction across the independent histidine assumptions. The other outcomes are retained, including contrary pH directions and conformation clashes. Coarse-model positivity cannot overrule contrary detailed results.

### Matched expanded protonation comparison

Positive contrast means the modeled neutral-pH proxy is less favorable than the acidic-pH proxy. Every row below uses the same structural-construction procedure and the same expanded-analysis method, including explicit introduced acidic groups.

| Internal record | Role | Minimum conditional contrast, kcal/mol | Explicit titrating sites | Maximum held-out energy error, kcal/mol |
|---|---|---:|---:|---:|
{chr(10).join(table)}

The differences support modest **conditional-model improvements**, not measured affinity or guaranteed pH specificity. New J00012, J00017, and J00018 models failed the stricter mouse geometry gate and were not rescued by relaxing that gate. A human pH score alone cannot qualify them for the recommended panel.

## 4. Classical protonation calculation and numerical validation

The calculation uses Amber99SB/OBC molecular mechanics, explicit neutral histidine tautomers and protonated histidines, and protonated/deprotonated carboxylate alternatives. Human and mouse coordinates are evaluated separately. Human force-field calculations use the recorded 334–505 target crop, while full-reference context, alternate human conformations, and resolved glycans are checked separately. Prior target-capping controls remain part of the evidence but are not falsely represented as an exhaustive full-receptor calculation.

This stage adds a **third-order finite cluster expansion**, a pre-2011 protein-energy strategy [3,4]. Point, pair, and triplet microstate energies determine the expansion. Exact tensor contraction then sums the proton-binding polynomial for each independent assumed-pKa combination. The polynomial is exact for the supplied energy tensor; the tensor itself is approximate where higher-order physical interactions are omitted. This is a physical energy interpolation, not a learned protein sequence or structure predictor.

Before using it on the new combinations, the implementation was tested against seven archived exhaustive datasets. Pair-order approximations were not selected after their larger discrepancies; the third-order form was retained. Each new model also received exact held-out energy evaluations, which were not used to fit the expansion. All **{sum(d['numerical_validation_pass'] for d in clusters.values())} of {len(clusters)}** completed models pass the predeclared 0.05 kcal/mol held-out discrepancy check. The maximum observed held-out discrepancy is **{summary['cluster_max_held_out_energy_error_kcal']:.5f} kcal/mol**. These finite checks do not establish a rigorous bound on every unobserved microstate.

Histidine free-state pKas are varied independently over 5.8, 6.3, and 6.8, with alternative neutral-tautomer fractions. Acidic-residue baseline priors and upward shifts are also varied independently. Salt/dielectric alternatives remain explicit. Assumed priors are **not predicted pKas**. Receptor-wide protonation, glycan ensembles, conformational entropy, explicit-water networks, and kinetic/detection effects are incomplete. Counts of microstates or assumption combinations are not independent experiments.

The numerical suite passes **{tests['passed']}/{tests['count']} checks**, including an independent direct log-sum-exp calculation, exact low-order synthetic cluster recovery, detection of omitted higher-order terms, and archived-result recovery. These tests validate arithmetic and implementation behavior, not biological activity.

## 5. Geometry, glycans, and unbound behavior

A uniform audit checks side-chain/bond geometry, chirality, peptide torsions, severe clashes, full same-species receptor context, and the checked human 1IVO/1NQL conformations. Failed models remain attached to their sequence evidence. A repaired geometry can support a sequence, but a contrary human pH direction from a valid geometry still excludes it.

The same archived finite glycan torsion grid is used across the pool. The final selection excludes candidates overlapping more than half of the admissible grid arrangements at any checked glycan root. **Those unweighted fractions are not glycan populations or binding probabilities.** Unresolved sugars, induced fit, and full glycan thermodynamics are not claimed to be modeled. Each final sequence has an explicit passing record, and the selection is tested for repeatability so that an excluded glycan case cannot reappear through a missing record.

Ten matched control/variant sequences completed unbound local minimization and passed native and independent geometry checks. These runs had no positional restraints but retained chirality restraints. They are local-relaxation diagnostics, not folding simulations, melting-temperature estimates, or entropy corrections. Earlier stage-4/5 short dynamics did not establish stabilization from the proposed loop mutations; the campaign's loop-preorganization risk therefore remains unresolved. No new thermal-stability claim is made for this release.

Two implementation issues were corrected transparently. First, single-chain free structures were exported as chain A, while the initial postprocessing wrapper expected complex chain B. All ten molecular calculations had completed; the wrapper's nonzero exits were retained and the actual exported chain independently audited without changing coordinates or energies. Second, glycan rescreening originally risked omitting previously excluded records. The corrected code screens the population before applying the glycan gate, and the final release checks complete coverage and idempotent selection. Failed source versions and logs remain archived.

## 6. Screening and ranking the actual panel

Every final sequence is unique, uses the 20 standard amino acids, is within the length limit, has exactly two cysteines, and lacks an N-X-S/T sequon. The latter two are conservative campaign design choices, not additional competition rules or proof of expression. Other sequence liabilities are reported rather than presented as validated developability predictions.

All 100 occur in the reconciled final search records, including the 133-sequence main search and a targeted completion search after the final model comparisons. Swiss-Prot, PDB sequences, and the augmented antibody collection each completed BLAST searches with explicit query-ID reconciliation and an exact protein positive control. Whole-chain framework hits are retained in the private evidence, not treated as blanket antibody novelty failures. The maximum local CDRH3 edit identity among the final 100 is **{100*v['local_CDR3_maximum_edit_identity']:.2f}%**, below the local 70% safeguard. These searches do not certify organizer novelty or exhaustive database coverage.

The final pool contains **{v['distinct_CDR3s']} distinct designed CDRH3 sequences** and **{v['distinct_docking_pose_ids']} recorded docking poses**. No CDRH3 appears in more than six panel entries. The ranking first separates evidence strata, then uses the conservative pH proxy within each stratum, with small predeclared redundancy penalties. These penalties are panel-design preferences, not physical free-energy corrections. Mouse pH selectivity is not used as an exclusion criterion.

Evidence-stratum counts across the 100 are **{v['evidence_tiers_all100']}**; across the top 20 they are **{v['evidence_tiers_top20']}**:

- **Tier 0:** checked human/mouse all-atom geometry and expanded engineered-acid/histidine analysis.
- **Tier 1:** checked human/mouse all-atom geometry, with the narrower histidine pH analysis.
- **Tier 2:** detailed human analysis, with mouse geometric prefilter only.
- **Tier 3:** exploratory reserves supported by initial geometry and a coarse pH model only. The campaign has observed false positives from that coarse model; these entries do not have the same evidential support as the leaders.

### Private top-20 codebook

The following internal identifiers are excluded from the public CSV. This table and the full evidence JSON are private campaign artifacts unless the submitter chooses to publish them.

| Rank | Public codename | Internal evidence record | Tier | Ranking contrast proxy, kcal/mol |
|---:|---|---|---:|---:|
{chr(10).join(top)}

## 7. Historical-method scope, provenance, and reproducibility

The accompanying `METHODS_CERTIFICATE.md` certifies the scientific method classes in the accepted final lineage, using the user's allowed modern implementations and datasets. No modern binder generator or learned protein predictor was run. The supplied modern mouse structure is input data, not a new prediction generated by this campaign.

An early PCG64 exploratory pilot predates the regenerated accepted lineage and is explicitly excluded. Consequently, the certificate does **not** claim that literally every exploratory operation throughout the campaign complied. It identifies the accepted MT19937 lineage and records the excluded pilot rather than rewriting its history.

All 100 final design records pass the saved ancestry audit, including sequence/coordinate identity, CDR masks, saved MT19937 loop inputs, and source hashes. The 6,240 stage-5 checkpoint payload checksums passed before continuation. All new code, held-out energies, negative results, structures, commands, seeds, lineage records, and corrected/failed versions are preserved. Large original input archives, unpacked runtimes, BLAST indexes, and reproducible caches are referenced rather than duplicated in the campaign checkpoint.

The optional `make_submission.py` is only an offline formatter for choosing the next non-rejected reserves. It does not perform or bypass official novelty assessment. Its **{ft['passed']} format/error-path tests** passed. Final FASTA and CSV round trips verify every sequence, name, order, class, and row count.

Floating-point/native-library differences can alter refinement trajectories even with the same recorded seeds. Restoring archived outputs is different from reproducing them bit-for-bit through a fresh calculation. No process is intended to remain running after checkpoint packaging; the final archive verification records payload hashes independently.

## 8. What remains experimentally unknown

The release does not establish any dissociation constant, acidic/neutral binding ratio, assay nondetection threshold, mouse affinity, expression yield, oligomeric state, folding stability, or therapeutic safety/efficacy. The strongest modeled contrast remains compatible with undesired neutral binding. Common-framework preservation and geometric checks reduce some avoidable risks but do not demonstrate a functional binder.

The file is an experimentally testable, ranked candidate panel developed under the historical-method constraint. It is not a certificate of biological success, official novelty acceptance, or selection for synthesis.

## Sources

[1] Challenge and FAQ, checked October 2, 2026: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr

[2] Updated Proteinbase novelty definitions, checked October 2, 2026: https://www.adaptyvbio.com/blog/novelty

[3] Cluster-expansion protein-energy precedent (2005): https://doi.org/10.1103/PhysRevLett.95.148103

[4] Ultra-Fast Evaluation of Protein Energies Directly from Sequence (2006): https://doi.org/10.1371/journal.pcbi.0020063

Additional historical primary-source identifiers and actual-versus-considered method use are in `source_register.json` and the prior-stage method audits. Numerical tables in this report are campaign computations, not measurements extracted from these sources.
'''
 (O/'CAMPAIGN_REPORT.md').write_text(report)
 sources=read(R/'stage5/reference/source_register.json');sources['stage6_additions']=[{'method':'Finite cluster expansion of protein energies','year':2005,'doi':'10.1103/PhysRevLett.95.148103','use':'Third-order explicit proton-state energy interpolation with held-out exact checks'},{'method':'Protein energy cluster expansion','year':2006,'doi':'10.1371/journal.pcbi.0020063','use':'Historical precedent, not a downloaded modern model'},{'method':'Exact finite tensor contraction / binding polynomial sum','year':'pre-2011 mathematical operation','use':'Modern numerical implementation, no learned inference'}];sources['current_rules']=read(S/'reference/rules_novelty_update.json');(O/'source_register.json').write_text(json.dumps(sources,indent=2))
 cert='''# Historical-method certificate for the delivered EGFR panel

## Certification scope

For the **100 delivered sequences and their accepted design/evaluation lineage**, I certify that the scientific method classes run or implemented are established in **2010 or earlier**, using the modern programming languages, implementations, and datasets expressly permitted by the user. No modern protein-binder generator, protein-specific learned sequence model, learned structure predictor, or learned interaction potential was run to generate or evaluate these designs.

Modern ChatGPT assisted planning, code authoring, and reporting. The amino-acid sequences were produced by the recorded classical construction and search code, not by running a modern protein generator.

This certification is about the recorded scientific methods and accepted ancestry. It is **not** a statement that the current chat system, operating system, compiler, package versions, file-format utilities, or data releases existed before 2011. It is not a certification of experimental binding, folding, safety, novelty acceptance, or competition eligibility.

## Methods actually used in the accepted lineage

| Scientific method class | Historical basis | Actual role |
|---|---|---|
| Generic VHH framework, native binding loops removed | Framework study/3EAK, 2008 | Structural scaffold only; not a known EGFR-binding loop seed |
| Internal-coordinate construction, Ramachandran/steric filters | 1963 and earlier geometry | Backbone construction and exclusion checks |
| Cyclic coordinate descent | Protein-loop implementation, 2003 | Independent loop closure |
| Finite side-chain rotamers, classical packing, simulated annealing | Pre-2011 methods | Sequence and side-chain searches |
| MT19937 | 1998 | Accepted stochastic sampling; seeds preserved |
| Powell optimization | 1964 | Selected side-chain geometry searches |
| Amber99SB | 2006 | All-atom molecular mechanics |
| OBC generalized Born | 2004 | Implicit solvent electrostatics |
| Explicit proton states and binding polynomials | Classical proton linkage; Wyman, 1964 | pH-direction diagnostics across stated assumptions |
| Finite cluster expansion | Protein-energy precedents, 2005/2006 | Third-order microstate energy approximation, independently tested |
| Exact finite sums/tensor contractions | Classical mathematics | Partition-polynomial evaluation, not a learned method |
| Kabsch alignment | 1976 | Structural comparisons |
| Verlet/Andersen/RATTLE approaches | 1982/1980/1983 | Earlier short classical dynamics controls; not folding proof |
| Finite glycan torsion grids and steric distances | Classical geometry | Resolved-glycan stress checks, not population inference |
| Gapped BLAST/BLOSUM62 | 1997/1992 | Local sequence search; modern BLAST implementation allowed |
| Edit-distance identity | Levenshtein, 1965/1966 | Local CDRH3 comparison; not modern IMGT annotation |
| Deterministic threshold filters and greedy diversity selection | Elementary pre-2011 operations | Explicit experimental-priority ordering |

The full source register distinguishes executed methods from ideas considered but not executed. Current data include the supplied mouse predicted structure, modern sequence databases, and downloaded software implementations. Their use follows the user's express allowance; no mouse structure prediction was run by this campaign.

## Explicit historical exception: excluded exploratory pilot

The early campaign used a modern default **PCG64** random-number generator in an exploratory pilot. That operation does not meet a literal pre-2011 algorithm restriction. It was disclosed, excluded, and superseded by a regenerated MT19937 lineage. The delivered sequences pass the saved audit linking them to the regenerated inputs.

Accordingly, I **cannot truthfully certify that every exploratory operation in the entire historical campaign complied**. The affirmative certificate is scoped to the accepted final design/evaluation lineage. The exception is preserved in the archive, not hidden or retrospectively relabeled.

## Remaining boundaries

ProteinTyper, ESMFold/ESMFold2, ANARCI, MMseqs2, Foldseek, and modern binder generators were not run locally for this panel. Adaptyv's external submission assessment is separate from this campaign's restricted scientific workflow. All final candidates remain experimentally unvalidated, and official novelty/classification remain pending.
'''
 (O/'METHODS_CERTIFICATE.md').write_text(cert)
 readme='''# EGFR ranked candidate files

`egfr_ranked_100.fasta` contains the requested 100 ranked, unique designs. These are unvalidated experimental hypotheses, not established binders.

**For Track 3, use `egfr_track3_top20.csv`.** It contains exactly 20 entries in rank order and only the required columns. Do not upload the entire 100-entry reserve CSV as one Track 3 submission. Official novelty/classification and other eligibility checks remain with Adaptyv; none has been performed through the submission portal here.

`egfr_track3_top20.fasta` contains the same top 20 in FASTA format. `egfr_100_RESERVE_POOL_NOT_SINGLE_UPLOAD.csv` is the full reserve pool for convenient replacement after official assessment.

The CSV and FASTA contain only public codenames and sequences. The codebook, full evidence JSON, report, source register, and method certificate are separate artifacts. Submitted methodology may be made public; opaque names are not a guarantee of secrecy.

## Replacing an officially rejected entry

The optional formatter preserves rank order, excludes listed public identifiers, and limits output to 20. It does not assess or bypass novelty. Put rejected identifiers, one per line, in `rejected_names.txt`, then use:

```bash
python make_submission.py --fasta egfr_ranked_100.fasta --reject-file rejected_names.txt --output egfr_track3_revised_top20.csv
```

The original CSV is ready without running this script. Always keep the total submitted count within the Track 3 limit; local formatting cannot track the account's existing submissions.

## Interpreting evidence

Read `CAMPAIGN_REPORT.md` and the evidence-tier column in `private_codebook_and_evidence.tsv`. Detailed-model leaders and coarse-only reserves do not have equivalent support. Human acidic binding, human neutral-pH nondetection, mouse binding, affinity, folding, and expression remain unverified for every sequence. Numerical scores are conditional proxies, not measured or calibrated binding free energies.

`METHODS_CERTIFICATE.md` gives the accepted-lineage historical-method certification and explicitly discloses the excluded early PCG64 pilot. Do not characterize the entire exploratory history as uniformly compliant.
'''
 (O/'README.md').write_text(readme);shutil.copy2(S/'code/make_submission.py',O/'make_submission.py');print('Wrote report, source registry, certificate, and README',len(report.split()),'report words')
if __name__=='__main__':main()
