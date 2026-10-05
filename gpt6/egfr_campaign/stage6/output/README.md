# EGFR ranked candidate files

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
