import os

ortho_txt = os.environ.get("ORTHOGROUPS_TXT", "Orthogroups.txt")
cds_db    = os.environ.get("CDS_DB", "all_species_cds.fasta")
id_list   = os.environ.get("ID_LIST", "RELAX_multicopy_241_IDs.txt")
outdir    = os.environ.get("WORKDIR", "RELAX_multicopy_241") + "/01_cds"
os.makedirs(outdir, exist_ok=True)

target = {}
for line in open(id_list):
    parts = line.split()
    if parts:
        target[parts[0]] = parts[1] if len(parts) > 1 else "NA"

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

buffers = {og: [] for og in target}
seen = set()
files = [cds_db]
def emit(name, seq):
    if name is None:
        return
    gid = name.split()[0]
    gid_clean = gid.lstrip(">")
    if gid_clean in gene2og:
        og = gene2og[gid_clean]
        buffers[og].append((gid_clean, "".join(seq)))
        seen.add(gid_clean)

for fp in files:
    if not os.path.isfile(fp):
        continue
    name, seq = None, []
    for line in open(fp, errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            emit(name, seq)
            name, seq = line, []
        elif name is not None:
            seq.append(line.strip())
    emit(name, seq)

n_written = 0
for og, recs in buffers.items():
    out_fp = os.path.join(outdir, og + ".cds.fasta")
    if os.path.exists(out_fp) and os.path.getsize(out_fp) > 0:
        continue
    if recs:
        with open(out_fp, "w") as o:
            for gid, s in recs:
                o.write(f">{gid}\n{s}\n")
        n_written += 1
print(f"[step0] 新写出 {n_written} 个家族（其余已存在），输出目录 {outdir}")
missing = [g for g, s in gene2og.items() if g not in seen]
if missing:
    print(f"[step0][警告] {len(missing)} 个基因在 CDS 库中未找到序列（前5个）: {missing[:5]}")