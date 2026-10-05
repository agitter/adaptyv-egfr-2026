# EGFR campaign: stage 3

**Date:** October 1, 2026. **Status:** computational improvement checkpoint, not a
final submission. Two Phone-a-Friend requests used; the third remains unused.

The output is a larger, better-audited candidate library and a narrower set of
promising physical hypotheses. No candidate has experimentally established human
binding, mouse binding, or undetectable binding at pH 7.4. The final ranked
100-sequence FASTA, secret public codenames, and top-20 competition CSV are still
pending. The interim FASTA in this archive is explicitly **NOT FOR SUBMISSION**.

## 1. Rules rechecked

The live challenge and FAQ were successfully retrieved on October 1. Human EGFR
binding, mouse cross-reactivity, and acidic-selective human binding remain the
three objectives. The explicit pH comparison is 6.5 versus 7.4, with no detectable
human binding at 7.4. The ranking order emphasizes pH, then mouse, then human
affinity. No assay detection threshold or required numerical dissociation
constant is published on the page. The FAQ also allows multidimensional winner
categories; this does not remove the user's requirement to pursue all three.

Track 3 allows at most 20 submitted designs. The requested 100 will be a ranked
pool with reserves. The required upload is a ranked CSV with name, sequence,
and molecule_class; FASTA remains the user's requested primary artifact. Current
lengths, 119-127 residues, fit the 10-250 limit. Generic nanobody frameworks are
allowed, but modifying a known binder is not. The current loops were generated
from scratch on the retained generic framework. The public deadline is
October 4, 2026, 23:59 Anywhere on Earth.

Sources: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
and https://www.adaptyvbio.com/blog/novelty . Retrieval records are preserved.

## 2. What changed in the search

The earlier coarse score omitted interactions between newly selected loop side
chains. The revised search includes intrabinder charge interactions, complete
framework atoms, target-associated glycan exclusions during rotamer selection,
and a mouse objective that does not incorrectly penalize stronger mouse binding.

A receptor-histidine docking search examined 71,640 orientations, accepted 51,
and retained 49 new poses in addition to 474 regenerated stage-2 poses. The new
sequence search then explored receptor-His recognition, mixed binder-His
interfaces, and explicit geometric contacts to both histidine ring nitrogens.

The two-ring-nitrogen requirement matters: a single designed hydrogen bond can
remain intact in one neutral histidine tautomer, undermining a putative pH
switch. The detailed calculation therefore includes HID, HIE and HIP, not just
one neutral form versus the protonated form.

Targeted A/B/C branches are mutations of the campaign's own untested designs,
not of an existing experimental binder. The A branch removes adverse introduced
histidines. B explores acidic contacts and side-chain conformations. C removes
nonessential positive charges and tests nearby acidic contacts while preserving
the leading two target-His contacts. A proposed additional binder-His contact
had zero acceptable geometric options; no speculative sequence was invented to
stand in for that failed search.

## 3. Library and audit counts

The machine-readable inventory is `output/stage3_summary.json`.

| Item | Result |
| --- | ---: |
| Sequence/structure records | 1,114 |
| Distinct sequences | 1,097 |
| Distinct CDR3 sequences | 1,044 |
| Length range | 119-127 residues |
| Initial human/mouse/glycan geometry prefilter pass | 1,114 |
| Records without an N-X-S/T sequon | 1,113 |
| Native-geometry detailed human MM/proton calculations | 49 |
| Those with acid-on direction across all 54 assumptions | 27 |
| Detailed mouse sequence evaluations | 7 |
| CDR3 records below the local 70% edit-identity threshold | 1,114 |
| Completed whole-chain searches, each of 3 databases | 61 records / 60 distinct sequences |

The 17 duplicate sequence records represent alternative structural hypotheses,
not additional candidates for the final 100. All five sequence/structure checks
per record passed: amino-acid alphabet, allowed length, coordinate-residue versus
sequence agreement, CDR slicing, and new-loop masks. The interim unique FASTA
collapses duplicates but deliberately retains poor designs for provenance.

All 49 detailed human refinements passed independent severe-clash, bond-length,
chirality and whole-receptor/glycan checks. An additional peptide-omega audit
flagged G00012_R00: one peptide bond deviates 31.34 degrees from trans. This is
retained as a warning, not silently discarded. The leading C00003/C00009 pair
has no such omega flag. A00008 has a glycosylation sequon and is excluded from
final candidacy regardless of its pH score.

The automatic exploratory ranking was also hardened: transformed designs with
unrecomputed parent scores or placeholder zero scores cannot enter it. Seven
regression checks passed. This changes no retained candidate or physical result;
the actual detailed batches were independently selected. It is not the final
competition ranking.

## 4. Numerical and mechanistic checks

The new direct Amber99SB/OBC interaction calculation agrees with OpenMM Reference
matched-partner subtraction to 7.83e-7 kcal/mol on the saved test complex.
Across six salt/dielectric models, comparison against OpenMM CustomGB gives a
maximum error of 2.07e-6 kcal/mol. Translation and partner-exchange invariance
hold to 3.8e-11 kJ/mol. At 1,000 nm partner separation, the residual interaction
is -0.0148 kcal/mol, consistent with the expected long-range approach to zero.
These are arithmetic checks, not proof of biological accuracy.

The proton ensemble passed 20 synthetic and Wyman-derivative checks. The mixed
acid/histidine ensemble passed five additional mathematical checks, and every
enumerated carboxylate charge-increment test passed. Surface-area tests recover
the isolated-sphere result and the correct decrease on overlap.

The independent acid-on FcRn/Fc control, PDB 1I1A, has the acid-on direction in
all 54 tested assumptions, with contrast 1.92-3.40 kcal/mol. Its sequence and
binding loops are used only as a control. The historical experimental mechanism
is described by Martin et al. (2001), DOI 10.1016/S1097-2765(01)00230-1. Correct
control direction is useful but is not calibration of predicted KD values.

## 5. Leading redesigns: improvement that survived additional checks

The G00117_C04 family targets human H358/H383. The aligned mouse target retains
a compatible positive interaction opportunity, including the human H383 to
mouse arginine substitution. The mouse is not required by the rules to have the
same pH switch.

C00003 replaces two nonessential loop lysines with glutamines (K36Q/K105Q) on
the campaign's own B00013 geometry. C00009 adds N52E to those changes. B00013 is
the same sequence as G00117_C04 with a different side-chain arrangement and
must not be counted as another distinct submitted sequence.

The following values are **conditional model contrasts**, defined as neutral-
pH minus acidic-pH interaction proxy. Positive values mean an acid-on direction
within this model. They are not measured free energies, confidence intervals,
probabilities, KD values, or proof of absent neutral-pH binding.

| Diagnostic, kcal/mol | Parent G00117_C04 | C00003 | C00009 |
| --- | ---: | ---: | ---: |
| Worst of 54 histidine/salt/dielectric assumptions | 0.596 | 0.869 | 0.937 |
| Worst after target-terminal capping | 0.594 | 0.854 | 0.943 |
| Worst with four receptor histidines, including sites outside the default cutoff | Not run | 0.853 | 0.947 |
| Worst of 162 histidine plus carboxylate assumptions | 0.577 | 0.766 | 0.855 |
| Worst on separate all-interface-HIP-refined geometry | Not run | 1.154 | 1.209 |

The broader site-cutoff audit explicitly includes H358, H370, H383 and H433
for each lead (81 microstates, 54 assumptions). Both retain the acid-on
direction. This and the carboxylate test are separate sensitivity tests, not a
fully joint enumeration of every titratable residue.

The carboxylate stress test includes both neutral oxygen-protonation forms.
C00003 uses 81 microstates; C00009 uses 243. Their pH trend survives these
additional states rather than relying on permanently charged designed acids.
The free pKa values and shifts remain declared assumptions, not computed or
measured pKa estimates.

At one explicitly fixed comparison scenario (solute dielectric 1, screening
parameter 1.25 nm^-1, unbound histidine pKa 6.3, neutral HIE fraction 0.5), the
human acidic interaction proxy improves from -5.46 for the parent to -7.00 for
C00003 and -6.42 for C00009. Mouse acidic proxies are -16.84, -19.34 and -19.46,
respectively. Neutral human proxies also remain favorable (-4.49, -5.23, -4.61),
so this is **not** evidence that neutral binding has been eliminated. The
increase in contrast is genuine within the calculation, but the competition's
no-detectable-neutral-binding criterion is unresolved.

Post-refinement whole-receptor surface calculations give approximately 491 A^2
interface area for C00003/C00009, using the conventional half-sum of buried
partner areas. C00003 changes by only 0.15% from 192 to 512 sphere points.
Neither lead has heavy-atom clashes below 2 A against the checked unliganded
1NQL or ligand-bound 1IVO receptor conformation. This resolves the specific
active-receptor collision seen in the earlier library; it does not establish
simultaneous EGF binding, and ligands were not part of that compatibility test.

## 6. Important failures and unresolved risk

The A-series histidine removals produce promising pH directions in some cases
but unfavorable acidic human interaction proxies of approximately +14 to +17
kcal/mol. They are not promoted just for looking pH selective. G00053_C05 has
better acidic interaction proxies, but its worst native-geometry contrast is
only 0.008 kcal/mol; it improves after protonated-geometry refinement and is
therefore geometry-sensitive. It also clashes extensively with active 1IVO.

Using a common interface-His set changes some results materially. For A00002,
including H383 just beyond the default 8-A cutoff reduces the minimum contrast
from approximately 0.644 to 0.345 kcal/mol. The old P00024_07 parent has a
negative minimum contrast (-0.589) when both its binder histidines and the three
common receptor sites are enumerated. Site-cutoff sensitivity is explicitly
retained, not hidden by comparing unmatched site sets.

**Loop preorganization is the clearest next problem.** Unbound position-
unrestrained, chirality-protected minimizations of C00003 and C00009 pass the
local geometry audit, with loop CA RMS displacements of 0.446 and 0.438 A.
However, the 10-ps unbound Verlet trajectory of their parent G00117_C04 shows
loop CA RMS displacement growing to 4.95 A, including 7.25 A in CDR2; the
framework remains at 0.98 A. This is substantial loop motion, despite passing
bond/clash/chirality checks. A local-minimization pass must not be called folding
or preorganization validation.

That trajectory has no thermostat or equilibration. Its initial kinetic energy
falls to approximately half, and it has a 26.32 kJ/mol initial total-energy
transient; after the first picosecond, sampled total-energy range is 2.53 kJ/mol.
It is not an equilibrated 300 K stability simulation. Nevertheless, the observed
motion is a warning that merits matched-control tests and loop stabilization.
The two charge-balanced mutants have not had equivalent dynamics runs; the
parent warning is a family risk, not a measured failure of each mutant.

## 7. Novelty and eligibility boundary

All 1,114 recorded CDR3s have nearest local edit identity below 70% against
203,818 reference segments; the maximum is 64.706%. Five exact-distance tests
and an exact-match positive control passed. The reference normalization retains
446,966 unique antibody sequences and original source pointers.

All three whole-sequence BLAST searches completed for the physically selected
61-record pool, including every charge-balance design. The Swiss-Prot, PDB and
augmented-antibody outputs each contain all 61 processed-query confirmations,
and their hashes match the execution record. None of the reported top hits is
an exact full-query match. Expected antibody framework homology is present.

A broader whole-antibody search was intentionally stopped and its output is
labeled PARTIAL. It must not be represented as completed. Full CDR3 screening
covers the entire 1,114-record library; whole-sequence completed searches cover
only the specified 61-record pool. The organizer uses different software,
numbering and database snapshots. Local motif-based annotation, incomplete
patent coverage, and the lack of the organizer's exact classifier prevent an
upload-acceptance guarantee. Generic-protein novelty thresholds must not be
substituted for the antibody-specific rules.

## 8. Methods, provenance and continuation

Scientific method classes in the current lineage satisfy the <=2010 constraint,
using the user's permitted modern implementations and current datasets. Read
`METHODS_AGE_AUDIT.md` for the exact scope, references, and the disclosed excluded
earlier PCG64 pilot. No modern learned generator, predictor or potential was run.

The checkpoint retains source code, accepted and rejected sequences, coordinates,
state energies, numerical controls, sensitivity analyses, logs and input/output
hash manifests. Original uploaded databases/software and regenerable indexes
are referenced rather than duplicated. Stage 1's large checkpoint was not
available in this runtime; inherited reports and reconstructed stage-2 methods
remain explicit, not falsely described as bit-identical recovery.

The next defined stage is to improve loop preorganization while preserving the
validated contact hypothesis, compare neutral and protonated alternative
geometries and unbound controls, expand distinct-sequence physical coverage,
and assemble the final diversity-aware ranked 100. Apparent pH selectivity alone
is not enough: both species' favorable contacts, neutral discrimination,
structural quality, sequence liabilities and novelty must remain separate
selection gates. No work is left running after the checkpoint is finalized.
