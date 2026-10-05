# EGFR campaign revision 2: evidence-led sequence refinement

Completed UTC: 2026-10-04T20:38:52.953713+00:00. Track 3. Two Phone-a-Friend requests used, one unused.

## 1. Deliverables and decision

The revised file contains **100 unique, ranked candidate sequences**, 119-127 residues long. **9 sequences in the pool are new in this revision; 0 are in the top 20.** The remaining retained sequences keep their original public identifiers and exact amino-acid sequences. File order defines rank; any numeric prefix in a retained identifier is an immutable name from the earlier release, not a revised rank.

No new sequence met the complete descriptive promotion rule. The original top-20 ordering is therefore retained; the revised reserve pool contains additional computationally characterized candidates rather than an invented new lead.

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
|---|---:|---:|
| Matched unchanged parent (M00000) [M00000] | 0.799256 | 0.005097 |
| BondVoyage_Cipher [M00018] | 0.762038 | 0.005090 |
| SilentCatalyst_Umbra [M00021] | 0.780473 | 0.004594 |

Positive contrast means an acid-on direction in this model. It is **not** an experimental binding free energy, dissociation constant, binding ratio or demonstration of neutral-pH nondetection. The matched unchanged control is the proper comparison; small differences from the prior J00008 score are additional-refinement effects, not mutation gains.

Both variants have slightly lower conservative contrast than the matched parent. Their interaction proxies also become more favorable at neutral pH, so stronger modeled attraction cannot be presented as a specificity improvement. The exact paired-parameter comparisons are:

| Record | Change in minimum contrast | Range of paired contrast changes, kcal/mol | Parameter combinations with larger contrast |
|---|---:|---:|---:|
| BondVoyage_Cipher | -0.037218 | -0.192660 to -0.001139 | 0 / 354,294 |
| SilentCatalyst_Umbra | -0.018783 | -0.132451 to 0.008982 | 10,049 / 354,294 |

The intended gain, if supported, is improved loop preorganization while retaining an acid-on model direction, **not** a stronger pH switch. No candidate is marked as satisfying the biological no-detectable-binding requirement.

## 5. Thermal tests and the artificial-restraint confound

The primary comparison consists of two separately seeded 20-picosecond production trajectories for each of the unchanged parent, M00018, M00021 and native fold control. Every run has a 2-picosecond warmup, 300 K target temperature, a 2-femtosecond timestep, Amber99SB/OBC, synchronized velocity-Verlet/RATTLE integration and Andersen-style constraint-group collisions generated by MT19937. Warmup position restraints are removed in production.

The inherited protocol adds artificial chirality restraints. Glycine replacements create additional stereocenters and therefore additional restraints, which could falsely appear as improved rigidity. Before any mutant thermal runs, the protocol was amended so the **primary comparison has no added chirality restraints**. Standard force-field improper terms remain. The four earlier added-restraint parent/native runs are retained as controls and are not substituted for the primary matched parent.

All saved CA-trace RMSDs were independently reconstructed by framework alignment and checked against the simulation records. The complete, unrestrained-production results are:

| Record | Seed | Mean loop RMSD, A | Late-half loop RMSD, A | Mean framework RMSD, A | Mean temperature, K |
|---|---:|---:|---:|---:|---:|
| Matched unchanged parent (M00000) | 198005 | 2.126 | 2.549 | 0.703 | 299.12 |
| Matched unchanged parent (M00000) | 198006 | 1.869 | 2.156 | 0.749 | 300.19 |
| BondVoyage_Cipher | 198005 | 1.734 | 2.113 | 0.719 | 298.17 |
| BondVoyage_Cipher | 198006 | 1.799 | 2.177 | 0.705 | 300.94 |
| SilentCatalyst_Umbra | 198005 | 1.762 | 2.207 | 0.796 | 299.40 |
| SilentCatalyst_Umbra | 198006 | 1.743 | 2.076 | 0.746 | 300.65 |
| NATIVE_3EAK | 198005 | 1.406 | 1.621 | 0.802 | 300.00 |
| NATIVE_3EAK | 198006 | 1.189 | 1.255 | 0.737 | 299.11 |

The predeclared descriptive rule requires at least 0.15 A lower mean loop displacement in **each** seed, no worsening in either late half, temperatures within 10 K of the matched parent, retained expanded pH contrast within 0.10 kcal/mol, and geometry/novelty safeguards. This rule is not a statistical significance test. Equal seed labels do not create identical atomic noise in different proteins.

| Variant | Mean loop reduction, seed 198005 / 198006, A | Late-half reduction, A | Descriptive thermal rule |
|---|---:|---:|---|
| BondVoyage_Cipher | 0.392 / 0.070 | 0.435 / -0.021 | Not met |
| SilentCatalyst_Umbra | 0.364 / 0.126 | 0.342 / 0.081 | Not met |

An additional descriptive check examines backbone movement at the six intended pH-contact positions. It was added during execution, before the complete comparison, and is not misrepresented as the original primary endpoint. Whole-loop RMSD can conceal worse movement at functional contacts.

| Variant | Contact-backbone reduction, seed 198005 / 198006, A |
|---|---:|
| BondVoyage_Cipher | 0.547 / 0.583 |
| SilentCatalyst_Umbra | 0.519 / 0.470 |

These are exceptionally short, nonequilibrated stress tests, not folding simulations, melting temperatures, thermodynamic stability estimates or evidence of expression. Changed amino-acid masses and the limited trajectory duration also limit interpretation. CA motion is not a direct measurement of side-chain hydrogen-bond occupancy or preservation of the binding mode.

### Phase-matched endpoint geometry

Applying the minimized-structure 25-degree peptide-angle cutoff directly to single thermal endpoints rejected all four initial control endpoints, including both native fold controls. This exposed a phase mismatch in the proposed audit, not a reason to call the native protein unfolded. The raw outcomes were preserved. Before any complete unrestrained trajectory was available, we recorded a correction: quench **all twelve** endpoints identically with 500 minimizer iterations and no added positional/chirality restraints, then apply the unchanged strict numerical geometry thresholds to those minimized endpoints. No loop-motion threshold or pH criterion was relaxed.

Raw endpoint strict passes: **0 / 12**. After phase-matched quenching, basic geometry passes: **12 / 12**; independent strict passes: **9 / 12**. All failed endpoints remain in the archive. These checks concern endpoints only; only CA traces, not full all-atom intermediate trajectories, were retained.
Some control endpoints still cross the heuristic peptide threshold after quenching. This remains a limitation of this force-field/geometry screen; it is not evidence that the experimentally determined native protein fails to fold. No thresholds were relaxed to remove these warnings. The affected records are:
- NATIVE_3EAK seed 198005 (source extra chirality restraints=True): 95-96:154.75deg
- NATIVE_3EAK seed 198005 (source extra chirality restraints=False): 55-56:152.81deg, 80-81:-149.90deg
- NATIVE_3EAK seed 198006 (source extra chirality restraints=False): 55-56:151.49deg

### The unbound-minimization warning remains

Before these thermal tests, local unbound minimization rearranged the two variants more than the parent:

| Record | Loop CA displacement on unbound minimization, A |
|---|---:|
| Matched unchanged parent (M00000) | 0.292952 |
| BondVoyage_Cipher | 0.565323 |
| SilentCatalyst_Umbra | 0.415182 |

Those local-reorganization results are retained even when another diagnostic is favorable. They do not by themselves measure folding free energy, but they prevent describing every assessment as an improvement.

## 6. Novelty and competition rules

All 23 distinct query sequences, including the unchanged control, completed searches against local Swiss-Prot, PDB sequence and augmented antibody databases with exact-reference positive controls and reconciled query identifiers. All 22 new sequences are below the local 70% CDRH3 edit-identity safeguard; the maximum is 50%. The local CDRH3 reference set contains 203,818 segments. Search inputs, results, positive controls, source mappings and hashes are retained.

The current Adaptyv policy has a separate antibody novelty scale: a conserved framework is not automatically disqualifying, and CDRH3 similarity is central once a sequence is classified as antibody-like [2]. We did not run the modern organizer classifier, IMGT numbering implementation, learned predictor, MMseqs2 or Foldseek. Local CDR extraction/edit-distance and BLAST are safeguards, not a reproduction of Adaptyv's exact pipeline. Official classification, novelty level and acceptance remain pending.

Track 3 permits **at most 20 designs**. Upload format is a ranked CSV with `name,sequence,molecule_class`; all delivered entries are intended nanobodies. The original de novo binding-loop lineage is preserved, but organizer eligibility is not self-certified. No public identifier contains an instruction, internal design label or methodological claim. Participant registration, rights, licensing and other personal attestations remain the submitter's responsibility.

The live challenge deadline has been extended to **October 6, 2026 at 23:59 AoE (UTC-12), which is October 7 at 06:59 in America/Chicago** [1]. Historical reports retain their earlier deadline observations; this revision corrects the current instruction. Human binding at pH6.5 with no detectable binding at pH7.4 remains the central criterion, with mouse cross-reactivity and human affinity also required. Mouse pH switching is not imposed by this campaign as an extra success criterion.

## 7. Selection, testing and provenance

The v2 pool contains **84 distinct designed CDRH3 sequences**, with at most **6** entries sharing a CDRH3. This is a count, not a claim that all structural hypotheses are independent. Evidence tiers in the 100 are: {'0': 13, '1': 21, '2': 3, '3': 63}. Tier 0 is expanded protonation plus both-species geometry, tier 1 narrower protonation plus both-species geometry, tier 2 detailed human evidence with more limited mouse checks, and tier 3 exploratory reserves. No tier certifies biological function.

Independent release validation passed **1022 checks**, including **51 coordinate/sequence identity checks**. Earlier in this revision, 369 generation/lineage/sequence-search checks, 104 finite-polynomial tests, six analytic alignment cases and 13 formatter/error-path tests passed. These counts describe implementation and consistency checks, not independent biological validation.

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
