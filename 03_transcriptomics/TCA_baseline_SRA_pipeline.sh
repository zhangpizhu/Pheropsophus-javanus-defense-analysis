#!/usr/bin/env bash
# ============================================================
# #11 验证版：T. castaneum 腺体 vs 腹部对照 的种内基线定量
# 目的：为 PJA_HEGs 提供"种内腺体特异性基线"，验证 21 个交集基因
#       的腺体上调不是物种间基线差异的假象
# 数据：PRJNA242780 (Li et al. 2013, PLoS Genet)
#   腺体组：  SRR1206008 (Male abdominal glands)
#            SRR1206010 (Female abdominal glands)
#   对照组：  SRR1206001 (Mid-abdomen，腹部非腺体组织)
# 注意：全部为 SINGLE-end 测序
# 重要：参考基因组与 GTF 必须使用您原分析处理 T. castaneum 时
#       所用的同一版本（保持流程完全一致，结果才可比）
# ============================================================
set -euo pipefail
THREADS=8
WORKDIR=~/tca_baseline          # 按需修改
REF_DIR=${WORKDIR}/ref          # 放您原来用的 T. castaneum 参考基因组
GENOME_FA=${REF_DIR}/tcastaneum.fa       # ← 改成您的实际文件名
GTF=${REF_DIR}/tcastaneum.gtf            # ← 改成您的实际文件名
HISAT2_INDEX=${REF_DIR}/tcastaneum       # ← 您原来建的 hisat2 index 前缀

mkdir -p ${WORKDIR}/{fastq,clean,bam,counts}
cd ${WORKDIR}

# ---------- 1. 下载 SRA（三选一） ----------
# 方案 A：ENA 直连（推荐，国内通常较快）
for R in SRR1206001 SRR1206008 SRR1206010; do
  wget -c -P fastq "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${R:0:6}/00${R: -1}/${R}/${R}.fastq.gz" \
    || wget -c -P fastq "https://ftp.sra.ebi.ac.uk/vol1/fastq/${R:0:6}/00${R: -1}/${R}/${R}.fastq.gz"
done
# 方案 B：sra-tools
# prefetch SRR1206001 SRR1206008 SRR1206010 && fasterq-dump --split-files -O fastq SRR12060*
# 方案 C：如果您本地还留着两个腺体样本的 clean 数据，只需下载 SRR1206001

# ---------- 2. Fastp 质控（与您原流程同版本同参数） ----------
for R in SRR1206001 SRR1206008 SRR1206010; do
  fastp -i fastq/${R}.fastq.gz -o clean/${R}.clean.fq.gz \
        -h clean/${R}.fastp.html -j clean/${R}.fastp.json -w ${THREADS}
done

# ---------- 3. HISAT2 比对（单端！用 -U 不是 -1/-2） ----------
for R in SRR1206001 SRR1206008 SRR1206010; do
  hisat2 -x ${HISAT2_INDEX} -U clean/${R}.clean.fq.gz -p ${THREADS} \
    | samtools sort -@ ${THREADS} -o bam/${R}.sorted.bam -
  samtools index bam/${R}.sorted.bam
done

# ---------- 4. featureCounts 计数（单端默认即可，参数同原流程） ----------
featureCounts -a ${GTF} -o counts/tca_counts.txt -T ${THREADS} \
  bam/SRR1206001.sorted.bam bam/SRR1206008.sorted.bam bam/SRR1206010.sorted.bam

echo "完成。输出：counts/tca_counts.txt —— 下一步运行 11_PJA_HEGs验证_比较分析.R"
