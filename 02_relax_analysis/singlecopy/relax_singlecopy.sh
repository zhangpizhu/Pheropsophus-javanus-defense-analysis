#!/bin/bash
# ============================================================================
# RELAX 批量流程：42 个单拷贝基因家族（对照组 controlB）
# 与多拷贝流程的区别：不建基因树，直接复用物种树拓扑；
#                     按每个家族比对中实际出现的物种裁剪物种树，再标 PJA 为 {TEST}
# 目录假设（脚本所在目录）：
#   controlB_stable_123_singlecopy_IDs.txt      家族 ID 清单（每行一个 OG 号）
#   iqtree.treefile                             物种树（未标记）
#   rename_cds_process_data/OGxxxxxxx_rename_trimal_cds_processed.fasta
#                                               已比对好的密码子 CDS（序列头=物种名）
# 用法：bash RELAX_单拷贝家族_流程脚本.sh
# ============================================================================
set -euo pipefail

# --------------------- 用户配置区 ---------------------
ID_LIST="controlB_stable_123_singlecopy_IDs.txt"
SPECIES_TREE="iqtree.treefile"
ALN_DIR="rename_cds_process_data"          # 密码子比对目录
ALN_SUFFIX="_rename_trimal_cds_processed.fasta"
TEST_SPECIES="PJA"                          # 标注 {TEST} 的物种名
WORKDIR="RELAX_singlecopy_42"
HYPHY_THREADS=4
# ------------------------------------------------------

mkdir -p ${WORKDIR}/{01_tree_test,02_relax_json,logs}
export SPECIES_TREE ALN_DIR ALN_SUFFIX TEST_SPECIES WORKDIR

while read -r OG; do
    [ -z "${OG}" ] && continue
    ALN="${ALN_DIR}/${OG}${ALN_SUFFIX}"
    if [ ! -f "${ALN}" ]; then
        echo "[跳过] ${OG} 无比对文件" | tee -a ${WORKDIR}/logs/skip.log
        continue
    fi

    # 步骤 1：裁剪物种树到该家族实际包含的物种 + 标 {TEST}
    python3 - "${ALN}" "${WORKDIR}/01_tree_test/${OG}.test.tree" "${OG}" <<'PYEOF'
import sys, os, re

aln_fp, out_fp, og = sys.argv[1], sys.argv[2], sys.argv[3]
tree_str = open(os.environ["SPECIES_TREE"]).read().strip()
test_sp  = os.environ["TEST_SPECIES"]

# 比对中出现的物种（序列头）
species = set()
for line in open(aln_fp, errors="ignore"):
    if line.startswith(">"):
        species.add(line[1:].split()[0])

# ---- 极简 Newick 解析（支持内部节点自举值与分支长度） ----
def parse(s):
    s = s.strip().rstrip(";")
    pos = 0
    def parse_node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            children = []
            while True:
                children.append(parse_node())
                if s[pos] == ",":
                    pos += 1
                else:
                    break
            assert s[pos] == ")"
            pos += 1
            label, blen = parse_label_len()
            return {"children": children, "label": label, "len": blen}
        else:
            label, blen = parse_label_len()
            return {"children": None, "label": label, "len": blen}
    def parse_label_len():
        nonlocal pos
        label = ""; blen = ""
        while pos < len(s) and s[pos] not in ",:)":
            label += s[pos]; pos += 1
        if pos < len(s) and s[pos] == ":":
            pos += 1
            while pos < len(s) and s[pos] not in ",()":
                blen += s[pos]; pos += 1
        return label, blen
    return parse_node()

def tips(node):
    if node["children"] is None:
        return [node["label"]]
    out = []
    for c in node["children"]:
        out += tips(c)
    return out

def prune(node, keep):
    """返回 (节点或None, 累加长度)"""
    if node["children"] is None:
        return (node, 0.0) if node["label"] in keep else (None, 0.0)
    newch = []
    for c in node["children"]:
        cc, _ = prune(c, keep)
        if cc is not None:
            newch.append(cc)
    if not newch:
        return None, 0.0
    if len(newch) == 1:                       # 单子收缩：长度累加
        ch = newch[0]
        try:
            ch["len"] = str(float(ch["len"] or 0) + float(node["len"] or 0))
        except ValueError:
            pass
        return ch, 0.0
    node["children"] = newch
    return node, 0.0

def serialize(node):
    if node["children"] is None:
        lbl = node["label"]
        if lbl == test_sp:
            lbl = lbl + "{TEST}"
        return lbl + (":" + node["len"] if node["len"] else "")
    inner = ",".join(serialize(c) for c in node["children"])
    return "(" + inner + ")" + (node["label"] or "") + (":" + node["len"] if node["len"] else "")

tree = parse(tree_str)
tree_tips = set(tips(tree))
missing = species - tree_tips
if missing:
    sys.exit(f"[错误] {og}: 比对中的物种 {missing} 不在物种树里")
keep = tree_tips & species
tree2, _ = prune(tree, keep)
newick = serialize(tree2) + ";"
assert "{TEST}" in newick, f"{og}: 未标到 {test_sp} 分支（该家族缺少 {test_sp}？）"
open(out_fp, "w").write(newick + "\n")
print(f"[step1] {og}: 树裁剪到 {len(keep)} 物种, TEST 已标注")
PYEOF

    # 步骤 2：hyphy relax
    hyphy relax --alignment "${ALN}" \
                --tree "${WORKDIR}/01_tree_test/${OG}.test.tree" \
                --test TEST CPU=${HYPHY_THREADS} \
                --output "${WORKDIR}/02_relax_json/${OG}.RELAX.json" \
                > ${WORKDIR}/logs/${OG}.hyphy.log 2>&1

    echo "[完成] ${OG}"
done < "${ID_LIST}"

# ============ 结果汇总：JSON -> TSV ============
python3 - <<'PYEOF'
import os, json, glob

workdir = os.environ.get("WORKDIR", "RELAX_singlecopy_42")
rows = []
for fp in sorted(glob.glob(os.path.join(workdir, "02_relax_json", "*.RELAX.json"))):
    og = os.path.basename(fp).replace(".RELAX.json", "")
    try:
        j = json.load(open(fp))
        tr = j.get("test results", {})
        k = tr.get("relaxation or intensification parameter", "")
        p = tr.get("p-value", "")
        lr = tr.get("LR", "")
        rows.append((og, k, lr, p))
    except Exception as e:
        rows.append((og, "PARSE_ERROR", "", str(e)[:60]))
out = os.path.join(workdir, "RELAX_summary.tsv")
with open(out, "w") as o:
    o.write("OG\tk\tLR\tp_value\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
print(f"[汇总] {len(rows)} 个家族 -> {out}")
PYEOF

echo "==================================================="
echo "完成。把 ${WORKDIR}/RELAX_summary.tsv（或整个 02_relax_json/）发回即可"
