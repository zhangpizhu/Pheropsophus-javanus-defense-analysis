#!/bin/bash
# ============================================================
# GOx 基因组缺失检验流程（回应 Reviewer 1 Minor Comment:
# "Have you tested whether GOx is actually absent from the
#  P. javanus genome (e.g., by tBLASTn search with known
#  insect GOx sequences)?"）
#
# 用法: bash GOx_tblastn_运行脚本.sh <基因组FASTA路径> [蛋白集FASTA路径]
# 示例: bash GOx_tblastn_运行脚本.sh Pjavanus.genome.fa Pjavanus.proteins.fa
# ============================================================

set -euo pipefail

GENOME=${1:? "请提供基因组FASTA路径作为第1个参数"}
PROTEIN=${2:-""}          # 可选：注释蛋白集FASTA
QUERY=GOx_query_sequences.fasta   # 与本脚本放同一目录
THREADS=16

# ---- 0. 检查 BLAST+ 是否可用（需 tblastn/makeblastdb/blastp）----
# 若未安装: conda install -c bioconda blast   或   module load blast
for tool in makeblastdb tblastn; do
    command -v $tool >/dev/null || { echo "错误: 未找到 $tool，请先安装 NCBI BLAST+"; exit 1; }
done

# ---- 1. 建库 ----
makeblastdb -in "$GENOME" -dbtype nucl -out Pjavanus_genome_db

# ---- 2. tBLASTn: 16条昆虫GOx/GMC家族蛋白 -> 基因组六框翻译 ----
# outfmt 6 列: qseqid sseqid pident length qcovs evalue bitscore sstart send
tblastn -query "$QUERY" -db Pjavanus_genome_db \
    -evalue 1e-5 -num_threads $THREADS \
    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore sstart send" \
    -out GOx_tblastn_genome.tsv

echo "=== tBLASTn 完成，命中数: $(wc -l < GOx_tblastn_genome.tsv) ==="

# ---- 3. (可选) BLASTp: 查询注释蛋白集 ----
if [ -n "$PROTEIN" ]; then
    makeblastdb -in "$PROTEIN" -dbtype prot -out Pjavanus_protein_db
    blastp -query "$QUERY" -db Pjavanus_protein_db \
        -evalue 1e-5 -num_threads $THREADS \
        -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
        -out GOx_blastp_proteome.tsv
    echo "=== BLASTp 完成，命中数: $(wc -l < GOx_blastp_proteome.tsv) ==="
fi

# ---- 4. 汇总每个查询序列的最佳命中（按bitscore排序，便于判读）----
echo "=== 每条查询序列的 Top3 命中（tBLASTn）==="
sort -k1,1 -k7,7nr GOx_tblastn_genome.tsv | awk '!seen[$1]++ || count[$1]<2 {print; count[$1]++}' | head -50

echo ""
echo "全部完成。请将 GOx_tblastn_genome.tsv（和 GOx_blastp_proteome.tsv）发回给我解读。"
