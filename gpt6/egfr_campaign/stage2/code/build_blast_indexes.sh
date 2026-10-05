#!/bin/sh
set -eu
S=/mnt/data/egfr_campaign/stage2
mkdir -p "$S/intermediate/blast"
for name in swissprot pdb; do
 "$S/runtime/bin/makeblastdb" -in "$S/intermediate/$name.fasta" -dbtype prot -out "$S/intermediate/blast/$name" > "$S/logs/makeblastdb_${name}.log" 2>&1
done
