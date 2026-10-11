#!/usr/bin/env bash
# Step 02 - BLAST survey of the six contracted families against the B. crepitans transcriptomes.
#
#   tblastn : P. javanus query proteins -> Bc contigs (forward search, E-value <= 1e-5)
#   blastx  : best-hit Bc contigs     -> P. javanus proteome (reciprocal-best-hit check)
#
# Requires: BLAST+ (tested with v2.12.0), output of 01_extract_and_prepare.py.
set -euo pipefail

DIR="/mnt/agents/output/A1_Bc_presence_survey"
DATA="$DIR/data"
RES="$DIR/results"
PJA_PEP="/mnt/agents/temp/PJA_PEP.fas"
EVALUE=1e-5
THREADS=3

cat "$DATA/Bc_glands.fna" "$DATA/Bc_combined.fna" > "$DATA/Bc_all.fna"

echo "== building nucleotide DB (Bc assemblies) =="
makeblastdb -in "$DATA/Bc_all.fna" -dbtype nucl -out "$DATA/Bc_db"

echo "== building protein DB (P. javanus proteome) =="
makeblastdb -in "$PJA_PEP" -dbtype prot -out "$DATA/PJA_db"

echo "== tblastn: 27 PJA queries vs Bc assemblies =="
tblastn -query "$DATA/pja_27_queries.faa" -db "$DATA/Bc_db" \
        -evalue $EVALUE -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
        -num_threads $THREADS -out "$RES/tblastn_pja_vs_Bc.tsv"

echo "== extract best-hit contigs per query =="
awk '!seen[$1]++' "$RES/tblastn_pja_vs_Bc.tsv" | cut -f2 > /tmp/bh_ids.txt
awk -v ids=/tmp/bh_ids.txt 'NR==FNR{keep[$1]=1;next} /^>/{h=substr($1,2);p=keep[h]} p' ids "$DATA/Bc_all.fna" > "$RES/best_hit_contigs.fna"

echo "== blastx: reciprocal-best-hit check vs PJA proteome =="
blastx -query "$RES/best_hit_contigs.fna" -db "$DATA/PJA_db" \
       -evalue $EVALUE -outfmt "6 qseqid sseqid pident length evalue bitscore" \
       -max_target_seqs 1 -num_threads $THREADS -out "$RES/blastx_Bc_vs_PJA.tsv"

echo "done. Results in $RES/"
wc -l "$RES/tblastn_pja_vs_Bc.tsv" "$RES/blastx_Bc_vs_PJA.tsv"
