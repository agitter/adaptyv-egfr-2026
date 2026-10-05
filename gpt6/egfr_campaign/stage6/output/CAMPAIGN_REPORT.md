# EGFR classical-design campaign: ranked 100-sequence release

**Release date:** October 2, 2026. **Track:** 3. **Status:** computational designs for experimental testing; no experimental binding or expression data. **Phone-a-Friend:** two requests used, one unused.

## 1. Deliverables and decision

The release contains **100 distinct, ranked amino-acid sequences**, with public codenames and separate private evidence. Their lengths are **119–127 residues**. All are intended single-domain nanobodies on a common VHH framework. The public upload files contain no internal design identifiers, instructions, or embedded prompts.

`egfr_ranked_100.fasta` is the requested complete reserve pool. `egfr_track3_top20.csv` is the ordered, 20-entry submission file. Its columns are exactly `name,sequence,molecule_class`, with `nanobody` as the class. `egfr_100_RESERVE_POOL_NOT_SINGLE_UPLOAD.csv` is a convenience export, **not a valid recommendation to upload all 100 in Track 3**. The top-20 FASTA is also supplied.

The ranking is an experimental-priority ordering, **not a prediction that any candidate wins or meets all three biological criteria**. Every candidate retains explicit UNVERIFIED fields for acidic human binding, neutral-pH nondetection, mouse binding, affinity, folding, and expression. Neutral-pH interaction-energy proxies remain favorable; no numerical score is relabeled as evidence of neutral-pH nondetection.

The first-ranked candidate is **001_SilentLigand_Onyx** (private source `stage6:J00008`). Its conservative ranking contrast is **0.8065 kcal/mol within the stated conditional model**, not a measured binding free energy or a calibrated dissociation constant.

## 2. Current rules and novelty policy

The live challenge and updated novelty article were reviewed on October 2, 2026 [1,2]. The challenge prioritizes human pH-selective binding, then mouse cross-reactivity, then human affinity. Human binding is required at pH 6.5, with no detectable binding at pH 7.4. Mouse binding is required, but mouse pH switching is not. The requested 100 is a reserve pool because Track 3 permits at most 20 submitted designs. The deadline is October 4, 2026, 23:59 AoE, corresponding to October 5 at 06:59 in America/Chicago.

The updated antibody-specific novelty scale uses CDRH3 edit identity. A value below 70% is the novel-CDRH3 branch; a familiar framework can lead to antibody novelty level 3 rather than a generic-protein rejection. We therefore **do not mutate conserved framework residues solely to lower whole-chain identity**, and do not apply the generic structural novelty scale to reject these intended nanobodies.

**Adaptyv's submission assessment is authoritative.** Our CDR boundary extraction is not the organizer's IMGT annotation, our BLAST searches are not its MMseqs2 pipeline, and our databases are not guaranteed to match every patent or current sequence snapshot. Neither official antibody classification nor upload acceptance has been established. No modern organizer classifier or structure predictor was run locally.

The generic framework is retained, but its original three binding loops were removed. Accepted loops trace to independently generated MT19937/CCD construction, with subsequent classical sequence and geometry design. Refining our own untested computational sequences is distinct from taking an existing experimental binder as the initial seed. The final lineage audit checks this recorded construction; it is not an organizer eligibility ruling. Participant registration, rights, licensing, and other personal eligibility attestations remain the submitter's responsibility.

## 3. What this stage actually added

Stage 6 generated **16 new unique sequences** in a 21-member contact-combination experiment: 16 new combinations and five unchanged-sequence controls. Every member was constructed on the same saved C00003-derived backbone, given matched refinement, and evaluated against both human and mouse targets. This removes unequal refinement histories as an explanation for favorable comparisons. The unchanged C00009 and C00017 sequence controls repair the earlier peptide-geometry problem on this shared backbone; the earlier failed coordinate models remain archived.

Ten existing designs from alternative docking poses received matched human and mouse calculations. Six additional existing sequences received missing mouse follow-up calculations. These K/L records are reassessments, **not additional unique sequences**. This stage has **68 completed detailed models**, of which **59 pass the combined geometry/context gate**. The underlying campaign contains **1,818 unique historical sequences**, including rejected designs; this is not a successful-binder count.

For the ten alternative poses, **1** currently pass both the checked human geometry and acid-on direction across the independent histidine assumptions. The other outcomes are retained, including contrary pH directions and conformation clashes. Coarse-model positivity cannot overrule contrary detailed results.

### Matched expanded protonation comparison

Positive contrast means the modeled neutral-pH proxy is less favorable than the acidic-pH proxy. Every row below uses the same structural-construction procedure and the same expanded-analysis method, including explicit introduced acidic groups.

| Internal record | Role | Minimum conditional contrast, kcal/mol | Explicit titrating sites | Maximum held-out energy error, kcal/mol |
|---|---|---:|---:|---:|
| J00000 | Unchanged C00003 parent | 0.5954 | 5 | 0.00006 |
| J00001 | Unchanged C00009 sequence | 0.6101 | 6 | 0.00008 |
| J00002 | Unchanged C00017 sequence | 0.7156 | 7 | 0.00012 |
| J00003 | Unchanged H00011 sequence | 0.6963 | 6 | 0.00531 |
| J00004 | Unchanged H00021 sequence | 0.6409 | 5 | 0.00002 |
| J00005 | New contact combination | 0.7529 | 7 | 0.00539 |
| J00006 | New contact combination | 0.7260 | 7 | 0.00739 |
| J00008 | New contact combination | 0.8065 | 8 | 0.00552 |
| J00009 | New contact combination | 0.7628 | 8 | 0.00568 |
| J00010 | New contact combination | 0.7402 | 7 | 0.00730 |
| J00014 | New contact combination | 0.8059 | 7 | 0.00556 |

The differences support modest **conditional-model improvements**, not measured affinity or guaranteed pH specificity. New J00012, J00017, and J00018 models failed the stricter mouse geometry gate and were not rescued by relaxing that gate. A human pH score alone cannot qualify them for the recommended panel.

## 4. Classical protonation calculation and numerical validation

The calculation uses Amber99SB/OBC molecular mechanics, explicit neutral histidine tautomers and protonated histidines, and protonated/deprotonated carboxylate alternatives. Human and mouse coordinates are evaluated separately. Human force-field calculations use the recorded 334–505 target crop, while full-reference context, alternate human conformations, and resolved glycans are checked separately. Prior target-capping controls remain part of the evidence but are not falsely represented as an exhaustive full-receptor calculation.

This stage adds a **third-order finite cluster expansion**, a pre-2011 protein-energy strategy [3,4]. Point, pair, and triplet microstate energies determine the expansion. Exact tensor contraction then sums the proton-binding polynomial for each independent assumed-pKa combination. The polynomial is exact for the supplied energy tensor; the tensor itself is approximate where higher-order physical interactions are omitted. This is a physical energy interpolation, not a learned protein sequence or structure predictor.

Before using it on the new combinations, the implementation was tested against seven archived exhaustive datasets. Pair-order approximations were not selected after their larger discrepancies; the third-order form was retained. Each new model also received exact held-out energy evaluations, which were not used to fit the expansion. All **11 of 11** completed models pass the predeclared 0.05 kcal/mol held-out discrepancy check. The maximum observed held-out discrepancy is **0.00739 kcal/mol**. These finite checks do not establish a rigorous bound on every unobserved microstate.

Histidine free-state pKas are varied independently over 5.8, 6.3, and 6.8, with alternative neutral-tautomer fractions. Acidic-residue baseline priors and upward shifts are also varied independently. Salt/dielectric alternatives remain explicit. Assumed priors are **not predicted pKas**. Receptor-wide protonation, glycan ensembles, conformational entropy, explicit-water networks, and kinetic/detection effects are incomplete. Counts of microstates or assumption combinations are not independent experiments.

The numerical suite passes **148/148 checks**, including an independent direct log-sum-exp calculation, exact low-order synthetic cluster recovery, detection of omitted higher-order terms, and archived-result recovery. These tests validate arithmetic and implementation behavior, not biological activity.

## 5. Geometry, glycans, and unbound behavior

A uniform audit checks side-chain/bond geometry, chirality, peptide torsions, severe clashes, full same-species receptor context, and the checked human 1IVO/1NQL conformations. Failed models remain attached to their sequence evidence. A repaired geometry can support a sequence, but a contrary human pH direction from a valid geometry still excludes it.

The same archived finite glycan torsion grid is used across the pool. The final selection excludes candidates overlapping more than half of the admissible grid arrangements at any checked glycan root. **Those unweighted fractions are not glycan populations or binding probabilities.** Unresolved sugars, induced fit, and full glycan thermodynamics are not claimed to be modeled. Each final sequence has an explicit passing record, and the selection is tested for repeatability so that an excluded glycan case cannot reappear through a missing record.

Ten matched control/variant sequences completed unbound local minimization and passed native and independent geometry checks. These runs had no positional restraints but retained chirality restraints. They are local-relaxation diagnostics, not folding simulations, melting-temperature estimates, or entropy corrections. Earlier stage-4/5 short dynamics did not establish stabilization from the proposed loop mutations; the campaign's loop-preorganization risk therefore remains unresolved. No new thermal-stability claim is made for this release.

Two implementation issues were corrected transparently. First, single-chain free structures were exported as chain A, while the initial postprocessing wrapper expected complex chain B. All ten molecular calculations had completed; the wrapper's nonzero exits were retained and the actual exported chain independently audited without changing coordinates or energies. Second, glycan rescreening originally risked omitting previously excluded records. The corrected code screens the population before applying the glycan gate, and the final release checks complete coverage and idempotent selection. Failed source versions and logs remain archived.

## 6. Screening and ranking the actual panel

Every final sequence is unique, uses the 20 standard amino acids, is within the length limit, has exactly two cysteines, and lacks an N-X-S/T sequon. The latter two are conservative campaign design choices, not additional competition rules or proof of expression. Other sequence liabilities are reported rather than presented as validated developability predictions.

All 100 occur in the reconciled final search records, including the 133-sequence main search and a targeted completion search after the final model comparisons. Swiss-Prot, PDB sequences, and the augmented antibody collection each completed BLAST searches with explicit query-ID reconciliation and an exact protein positive control. Whole-chain framework hits are retained in the private evidence, not treated as blanket antibody novelty failures. The maximum local CDRH3 edit identity among the final 100 is **55.56%**, below the local 70% safeguard. These searches do not certify organizer novelty or exhaustive database coverage.

The final pool contains **85 distinct designed CDRH3 sequences** and **8 recorded docking poses**. No CDRH3 appears in more than six panel entries. The ranking first separates evidence strata, then uses the conservative pH proxy within each stratum, with small predeclared redundancy penalties. These penalties are panel-design preferences, not physical free-energy corrections. Mouse pH selectivity is not used as an exclusion criterion.

Evidence-stratum counts across the 100 are **{'0': 11, '1': 14, '2': 3, '3': 72}**; across the top 20 they are **{'0': 11, '1': 9}**:

- **Tier 0:** checked human/mouse all-atom geometry and expanded engineered-acid/histidine analysis.
- **Tier 1:** checked human/mouse all-atom geometry, with the narrower histidine pH analysis.
- **Tier 2:** detailed human analysis, with mouse geometric prefilter only.
- **Tier 3:** exploratory reserves supported by initial geometry and a coarse pH model only. The campaign has observed false positives from that coarse model; these entries do not have the same evidential support as the leaders.

### Private top-20 codebook

The following internal identifiers are excluded from the public CSV. This table and the full evidence JSON are private campaign artifacts unless the submitter chooses to publish them.

| Rank | Public codename | Internal evidence record | Tier | Ranking contrast proxy, kcal/mol |
|---:|---|---|---:|---:|
| 1 | 001_SilentLigand_Onyx | stage6:J00008 | 0 | 0.8065 |
| 2 | 002_Stereosecret_Cipher | stage6:J00014 | 0 | 0.8059 |
| 3 | 003_CatalyticCaper_Cinder | stage6:J00002 | 0 | 0.7156 |
| 4 | 004_ChiralCipher_Cipher | stage6:J00009 | 0 | 0.7628 |
| 5 | 005_PolypeptidePlot_Cinder | stage6:J00005 | 0 | 0.7529 |
| 6 | 006_BondVoyage_Nightfall | stage5:H00021 | 0 | 0.6409 |
| 7 | 007_IsoelectricAlias_Nightfall | stage6:J00010 | 0 | 0.7402 |
| 8 | 008_MotifMisdirect_Cinder | stage6:J00006 | 0 | 0.7260 |
| 9 | 009_IsoelectricAlias_Cipher | stage6:J00001 | 0 | 0.6101 |
| 10 | 010_SilentLigand_Cinder | stage5:N2070600 | 0 | 0.6013 |
| 11 | 011_IsoelectricAlias_Velvet | stage5:B00000 | 0 | 0.5954 |
| 12 | 012_PeptidePhantom_Nightfall | stage4:S00048 | 1 | 0.8391 |
| 13 | 013_SilentLigand_Velvet | stage5:H00018 | 1 | 0.8290 |
| 14 | 014_MolarMirage_Cinder | stage6:J00020 | 1 | 0.8504 |
| 15 | 015_MolarMirage_Velvet | stage6:J00015 | 1 | 0.8407 |
| 16 | 016_EnzymeEnigma_Onyx | stage6:J00013 | 1 | 0.8377 |
| 17 | 017_Covertase_Velvet | stage4:S00022 | 1 | 0.7909 |
| 18 | 018_AmideAlibi_Nightfall | stage5:T3080600 | 1 | 0.6632 |
| 19 | 019_MolarMirage_Cipher | stage5:T3070900 | 1 | 0.6483 |
| 20 | 020_CatalyticCaper_Cipher | stage6:J00019 | 1 | 0.7450 |

## 7. Historical-method scope, provenance, and reproducibility

The accompanying `METHODS_CERTIFICATE.md` certifies the scientific method classes in the accepted final lineage, using the user's allowed modern implementations and datasets. No modern binder generator or learned protein predictor was run. The supplied modern mouse structure is input data, not a new prediction generated by this campaign.

An early PCG64 exploratory pilot predates the regenerated accepted lineage and is explicitly excluded. Consequently, the certificate does **not** claim that literally every exploratory operation throughout the campaign complied. It identifies the accepted MT19937 lineage and records the excluded pilot rather than rewriting its history.

All 100 final design records pass the saved ancestry audit, including sequence/coordinate identity, CDR masks, saved MT19937 loop inputs, and source hashes. The 6,240 stage-5 checkpoint payload checksums passed before continuation. All new code, held-out energies, negative results, structures, commands, seeds, lineage records, and corrected/failed versions are preserved. Large original input archives, unpacked runtimes, BLAST indexes, and reproducible caches are referenced rather than duplicated in the campaign checkpoint.

The optional `make_submission.py` is only an offline formatter for choosing the next non-rejected reserves. It does not perform or bypass official novelty assessment. Its **13 format/error-path tests** passed. Final FASTA and CSV round trips verify every sequence, name, order, class, and row count.

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
