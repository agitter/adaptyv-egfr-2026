# EGFR revised candidate files (v2)

100 unique ranked hypotheses, 119-127 residues. 9 new sequences in this pool; 0 new sequences in its top20. None is experimentally validated.

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
