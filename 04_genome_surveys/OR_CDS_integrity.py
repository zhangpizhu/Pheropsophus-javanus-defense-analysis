# ============================================================
# #10 假基因检查 · 步骤2：CDS 完整性核查（内部终止子 / 移码）
# 用法: python OR_CDS完整性核查_步骤2.py 你的CDS文件.fasta
# 输入: 从基因组注释中提取的 CDS fasta（gffread -x 或类似工具生成）
# 输出: 终端表格 + OR_CDS_check.txt
# 判读: 内部终止子 >=1 或 长度%3!=0 → 假基因化/组装错误候选
# ============================================================
import sys

def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]; seqs[name] = ''
        elif line:
            seqs[name] += line.upper()
    return seqs

CODON = {}
bases = 'TCAG'
aa1 = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a+b+c] = aa1[i]; i += 1

def translate(s):
    return ''.join(CODON.get(s[i:i+3], 'X') for i in range(0, len(s)-2, 3))

TARGETS = ["Chr15.729","Chr7.150","Chr7.158","Chr7.159","Chr10.410","Chr7.77",
           "Chr6.790","Chr4.822","Chr4.829","Chr14.471","Chr14.469","Chr6.874",
           "Chr6.873","Chr17.166","Chr1.365","Chr1.366","Chr1.362","Chr1.364",
           "Chr6.872","Chr6.875"]

seqs = read_fasta(sys.argv[1])
hits = {n: s for n, s in seqs.items() if any(t in n for t in TARGETS)}
print(f"在 CDS 文件中匹配到 {len(hits)}/20 个目标基因\n")

lines = ["gene\tcds_len\tlen%3\tinternal_stops\tverdict"]
for n, s in sorted(hits.items()):
    prot = translate(s)
    internal = prot[:-1].count('*')          # 末尾终止子正常，只数内部的
    rem = len(s) % 3
    verdict = "OK" if (internal == 0 and rem == 0) else \
              ("PSEUDOGENE?" if internal > 0 else "PARTIAL_CDS")
    line = f"{n}\t{len(s)}\t{rem}\t{internal}\t{verdict}"
    print(line); lines.append(line)

missing = [t for t in TARGETS if not any(t in n for n in hits)]
if missing:
    print("\n未匹配（核对ID前缀）:", ", ".join(missing))

with open("OR_CDS_check.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n>> 结果已保存: OR_CDS_check.txt")
