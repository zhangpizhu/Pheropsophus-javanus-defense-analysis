import os

workdir = os.environ.get("WORKDIR", "RELAX_multicopy_241")
ortho_dir = os.environ.get("ORTHO_SEQ_DIR", ".")
id_list = os.environ.get("ID_LIST", "RELAX_multicopy_241_IDs.txt")
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
log = open(os.path.join(workdir, "logs", "dropped_genes.log"), "a")
n_ok_fam = 0
for og in ogs:
    aa_fp = os.path.join(ortho_dir, og + ".fa")
    if not os.path.exists(aa_fp):
        aa_fp = os.path.join(ortho_dir, og + ".fasta")
    cds_fp = os.path.join(cds_in, og + ".cds.fasta")
    if not (os.path.exists(aa_fp) and os.path.exists(cds_fp)):
        continue
    if os.path.exists(os.path.join(aa_out, og + ".fa")) and os.path.getsize(os.path.join(aa_out, og + ".fa")) > 0:
        n_ok_fam += 1
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