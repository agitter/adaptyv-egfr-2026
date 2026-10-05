#!/bin/sh
set -eu
S=/mnt/data/egfr_campaign/stage2
for name in swissprot pdb antibodies; do
 "$S/runtime/bin/blastp" -query "$S/intermediate/novelty/candidates.fasta" -db "$S/intermediate/blast/$name" -out "$S/intermediate/novelty/blast_${name}.tsv" -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp' -evalue 0.01 -max_target_seqs 20 -num_threads 2 -seg no > "$S/logs/blast_${name}.log" 2>&1
done
