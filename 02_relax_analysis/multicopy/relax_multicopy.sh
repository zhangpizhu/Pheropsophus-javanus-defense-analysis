#!/bin/bash
# ============================================================================
# RELAX 批量流程：241 个多拷贝基因家族
# 流程：提取序列 -> muscle(AA) -> trimAL(AA) -> pal2nal(AA->CDS密码子比对)
#       -> IQ-TREE 建基因树 -> PJA 分支标注 {TEST} -> hyphy relax
# 用法：bash RELAX_多拷贝家族_流程脚本.sh
# ============================================================================
set -euo pipefail

# --------------------- 用户配置区（按实际情况修改） --------------------------
ORTHO_SEQ_DIR="/path/to/Orthogroup_Sequences"   # OrthoFinder 输出：每家族一个 OGxxxxxxx.fa（蛋白）
ORTHOGROUPS_TXT="/path/to/Orthogroups.txt"      # OrthoFinder 输出：家族-基因ID 对应表
CDS_DB="/path/to/all_species_cds.fasta"         # 所有物种 CDS 合并的单个大 fasta 文件
                                                # 序列头格式：物种名|转录本ID（与 Orthogroups.txt 中一致，如 PJA|Chr1.618）
ID_LIST="multicopy_241_family_IDs.txt"            # 本脚本同目录，两列：家族ID <tab> 组别
WORKDIR="RELAX_multicopy_241"                    # 输出工作目录
PJA_PREFIX="PJA|"                                # P. javanus 序列头前缀（精确到竖线，避免误标）
THREADS=8                                        # IQ-TREE / 并行用线程数
HYPHY_THREADS=4
# ----------------------------------------------------------------------------

mkdir -p ${WORKDIR}/{01_cds,02_muscle,03_trimAL,04_pal2nal,05_iqtree,06_tree_test,07_relax_json,logs}
export ORTHOGROUPS_TXT CDS_DB ID_LIST WORKDIR ORTHO_SEQ_DIR   # 传递给内嵌 Python

# ============ 步骤 0：一次性从各物种 CDS 库中提取 241 个家族的 CDS ============
# 只扫描 CDS 库一遍，按 Orthogroups.txt 的基因->家族映射分拣输出
python3 - <<'PYEOF'
import os, glob, re

ortho_txt = os.environ.get("ORTHOGROUPS_TXT", "Orthogroups.txt")
cds_db    = os.environ.get("CDS_DB", "all_species_cds.fasta")
id_list   = os.environ.get("ID_LIST", "multicopy_241_family_IDs.txt")
outdir    = os.environ.get("WORKDIR", "RELAX_multicopy_241") + "/01_cds"

# 需要的家族集合（含组别）
target = {}
for line in open(id_list):
    parts = line.split()
    if parts:
        target[parts[0]] = parts[1] if len(parts) > 1 else "NA"

# 基因 -> 家族 映射（仅保留目标家族）
# OrthoFinder 的 Orthogroups.txt 实际格式：每行 "OG0000000: 基因1 基因2 基因3 ..."（空格分隔，无表头）
gene2og = {}
with open(ortho_txt, errors="ignore") as f:
    for line in f:
        if ":" not in line:
            continue
        og, genes = line.split(":", 1)
        og = og.strip()
        if og not in target:
            continue
        for gene in genes.split():
            if gene:
                gene2og[gene] = og
print(f"[step0] 目标家族 {len(target)} 个，涉及基因 {len(gene2og)} 个")

# 扫描 CDS 库并分拣
buffers = {og: [] for og in target}
seen = set()
files = [cds_db]   # 单一大 fasta 文件
def flush(name, seq, og):
    buffers[og].append((name, seq))

for fp in files:
    if not os.path.isfile(fp):
        continue
    name, seq = None, []
    def emit():
        if name is None:
            return
        gid = name.split()[0]                    # 取序列头第一个词作为基因ID
        gid_clean = gid.lstrip(">")
        if gid_clean in gene2og:
            og = gene2og[gid_clean]
            buffers[og].append((gid_clean, "".join(seq)))
            seen.add(gid_clean)
    for line in open(fp, errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            emit()
            name, seq = line, []
        elif name is not None:
            seq.append(line.strip())
    emit()

for og, recs in buffers.items():
    if recs:
        with open(os.path.join(outdir, og + ".cds.fasta"), "w") as o:
            for gid, s in recs:
                o.write(f">{gid}\n{s}\n")
n_written = sum(1 for og in target if buffers[og])
print(f"[step0] 已输出 {n_written}/{len(target)} 个家族的 CDS fasta -> {outdir}")
missing = [g for g, s in gene2og.items() if g not in seen]
if missing:
    print(f"[step0][警告] {len(missing)} 个基因在 CDS 库中未找到序列（前5个）: {missing[:5]}")
PYEOF

# ====== 步骤 0.5：蛋白-CDS 一致性预检（剔除注释不一致的基因） ======
# 判定：CDS 去终止子后翻译 == 蛋白序列（逐残基严格一致）
# 不一致的基因从该家族的 AA 和 CDS 文件中同时剔除，并记入 logs/dropped_genes.log
python3 - <<'PYEOF'
import os

workdir = os.environ.get("WORKDIR", "RELAX_multicopy_241")
ortho_dir = os.environ.get("ORTHO_SEQ_DIR", ".")
id_list = os.environ.get("ID_LIST", "multicopy_241_family_IDs.txt")
cds_in = os.path.join(workdir, "01_cds")
aa_out = os.path.join(workdir, "01_aa_filtered")
cds_out = os.path.join(workdir, "01_cds_filtered")
os.makedirs(aa_out, exist_ok=True)
os.makedirs(cds_out, exist_ok=True)

CODON = {}
b = "TCAG"; a1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for x in b:
    for y in b:
        for z in b:
            CODON[x+y+z] = a1[i]; i += 1
def tr(s):
    s = s.upper()
    return "".join(CODON.get(s[j:j+3], "X") for j in range(0, len(s)-2, 3))

def read_fasta(fp):
    d = {}; name = None
    for line in open(fp, errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; d[name] = []
        elif name:
            d[name].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}

ogs = [l.split()[0] for l in open(id_list) if l.strip()]
log = open(os.path.join(workdir, "logs", "dropped_genes.log"), "w")
n_ok_fam = 0
for og in ogs:
    aa_fp = os.path.join(ortho_dir, og + ".fa")
    if not os.path.exists(aa_fp):
        aa_fp = os.path.join(ortho_dir, og + ".fasta")
    cds_fp = os.path.join(cds_in, og + ".cds.fasta")
    if not (os.path.exists(aa_fp) and os.path.exists(cds_fp)):
        continue
    aa = read_fasta(aa_fp); cds = read_fasta(cds_fp)
    keep = []
    for gid, prot in aa.items():
        nuc = cds.get(gid)
        if nuc is None:
            log.write(f"{og}\t{gid}\tno_CDS\n"); continue
        prot_tr = tr(nuc).rstrip("*")
        if prot_tr == prot.rstrip("*"):
            keep.append(gid)
        else:
            log.write(f"{og}\t{gid}\ttranslation_mismatch(pep={len(prot)},cds_aa={len(prot_tr)})\n")
    if len(keep) >= 10:
        with open(os.path.join(aa_out, og + ".fa"), "w") as o:
            for g in keep: o.write(f">{g}\n{aa[g]}\n")
        with open(os.path.join(cds_out, og + ".cds.fasta"), "w") as o:
            for g in keep: o.write(f">{g}\n{cds[g]}\n")
        n_ok_fam += 1
    else:
        log.write(f"{og}\tALL\ttoo_few_consistent_genes({len(keep)})\n")
log.close()
print(f"[step0.5] 一致性预检完成：{n_ok_fam} 个家族通过（>=10 条一致序列），详见 logs/dropped_genes.log")
PYEOF

# ============ 步骤 1-6：逐家族处理 ============
cut -f1 ${ID_LIST} | while read -r OG; do
    AA_RAW="${WORKDIR}/01_aa_filtered/${OG}.fa"
    CDS="${WORKDIR}/01_cds_filtered/${OG}.cds.fasta"
    if [ ! -f "${AA_RAW}" ] || [ ! -f "${CDS}" ]; then
        echo "[跳过] ${OG} 缺少蛋白或 CDS 文件" | tee -a ${WORKDIR}/logs/skip.log
        continue
    fi

    # 步骤 1：muscle 蛋白比对
    muscle -align "${AA_RAW}" -output "${WORKDIR}/02_muscle/${OG}.aa.aln" \
        > ${WORKDIR}/logs/${OG}.muscle.log 2>&1

    # 步骤 2：pal2nal 转换（必须先于 trimAL！pal2nal 校验 蛋白==CDS翻译，
    #          若先 trimAL 删列，蛋白序列残缺，必然校验失败）
    pal2nal.pl "${WORKDIR}/02_muscle/${OG}.aa.aln" "${CDS}" \
        -output fasta -nogap -nomismatch \
        > "${WORKDIR}/04_pal2nal/${OG}.codon.full.fasta" 2> ${WORKDIR}/logs/${OG}.pal2nal.log

    # 步骤 3：trimAL 在 AA 层面选择保留列（-automated1）
    trimal -in "${WORKDIR}/02_muscle/${OG}.aa.aln" \
           -out "${WORKDIR}/03_trimAL/${OG}.aa.trim.aln" \
           -automated1 >> ${WORKDIR}/logs/${OG}.trimal.log 2>&1

    # 步骤 4：把 trimAL 在 AA 层面的保留列映射到密码子比对
    #         （AA 第 i 列 <-> 密码子第 3i..3i+2 列；逐整 codon 修剪，不破坏阅读框）
    python3 - "${WORKDIR}/02_muscle/${OG}.aa.aln" \
              "${WORKDIR}/03_trimAL/${OG}.aa.trim.aln" \
              "${WORKDIR}/04_pal2nal/${OG}.codon.full.fasta" \
              "${WORKDIR}/04_pal2nal/${OG}.codon.fasta" <<'PYEOF'
import sys
aa_raw_fp, aa_trim_fp, codon_fp, out_fp = sys.argv[1:5]

def read_fasta(fp):
    d = {}; order = []; name = None
    for line in open(fp, errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; d[name] = []; order.append(name)
        elif name:
            d[name].append(line.strip())
    return order, {k: "".join(v) for k, v in d.items()}

order, aa_raw = read_fasta(aa_raw_fp)
_, aa_trim = read_fasta(aa_trim_fp)
_, codon = read_fasta(codon_fp)

raw_rows = [aa_raw[k] for k in order]
trim_rows = [aa_trim[k] for k in order]
n_raw = len(raw_rows[0]); n_trim = len(trim_rows[0])

# 1) 原 AA 比对中哪些列被 trimAL 保留（整列匹配，trim 列是 raw 列的子序列）
kept_raw = []; j = 0
for i in range(n_raw):
    if j < n_trim and all(raw_rows[k][i] == trim_rows[k][j] for k in range(len(order))):
        kept_raw.append(i); j += 1
assert j == n_trim, "trimAL 列匹配失败"

# 2) pal2nal -nogap 已删除含 gap 的列：raw 中无 gap 列 <=> codon 列（一一对应）
gapless = [i for i in range(n_raw) if all(raw_rows[k][i] != "-" for k in range(len(order)))]
kept_gapless = [i for i in kept_raw if i in set(gapless)]
codon_cols = [gapless.index(i) for i in kept_gapless]   # codon 比对中的列号

# 3) 逐 codon 修剪并输出
with open(out_fp, "w") as o:
    for k in order:
        seq = codon[k]
        trimmed = "".join(seq[3*c:3*c+3] for c in codon_cols)
        o.write(f">{k}\n{trimmed}\n")
print(f"[step4] {out_fp}: 保留 {len(codon_cols)} 个密码子列")
PYEOF

    # 步骤 4：IQ-TREE 建基因树（自动选模型）
    iqtree2 -s "${WORKDIR}/04_pal2nal/${OG}.codon.fasta" \
            -m MFP -T ${THREADS} --quiet \
            -pre "${WORKDIR}/05_iqtree/${OG}" \
            >> ${WORKDIR}/logs/${OG}.iqtree.log 2>&1

    # 步骤 5：PJA 末梢分支标注 {TEST}
    python3 - "${WORKDIR}/05_iqtree/${OG}.treefile" \
                "${WORKDIR}/06_tree_test/${OG}.test.tree" "${PJA_PREFIX}" <<'PYEOF'
import sys, re
src, dst, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
tree = open(src).read().strip()
# 在 "PJA...名字" 之后、"冒号+分支长度" 之前插入 {TEST}
labeled = re.sub(r'(' + re.escape(prefix) + r'[^,():\s]*)(?=:)', r'\1{TEST}', tree)
n = labeled.count('{TEST}')
assert n > 0, f"未标注到任何 {prefix} 分支，请检查序列头前缀"
open(dst, "w").write(labeled + "\n")
print(f"[step5] {src}: 标注 {n} 个 TEST 分支")
PYEOF

    # 步骤 6：hyphy relax
    hyphy relax --alignment "${WORKDIR}/04_pal2nal/${OG}.codon.fasta" \
                --tree "${WORKDIR}/06_tree_test/${OG}.test.tree" \
                --test TEST CPU=${HYPHY_THREADS} \
                --output "${WORKDIR}/07_relax_json/${OG}.RELAX.json" \
                > ${WORKDIR}/logs/${OG}.hyphy.log 2>&1

    echo "[完成] ${OG}"
done

echo "==================================================="
echo "全部完成。结果 JSON 在 ${WORKDIR}/07_relax_json/"
echo "下一步：把 07_relax_json 整个目录打包发回，做三组 k/P 汇总统计"
