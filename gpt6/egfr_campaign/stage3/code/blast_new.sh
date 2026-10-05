#!/bin/sh
set -eu
R=/mnt/data/egfr_campaign
S="$R/stage3"
mkdir -p "$S/intermediate/novelty"
for db in antibodies_augmented swissprot pdb; do
 "$R/stage2/runtime/bin/blastp" -query "$S/intermediate/novelty/candidates.fasta" -db "$R/stage2/intermediate/blast/$db" -out "$S/intermediate/novelty/${db}_hits.tsv" -outfmt '6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send' -max_target_seqs 10 -num_threads 2 -seg no -evalue 0.001 > "$S/logs/blast_${db}.log" 2>&1
done
