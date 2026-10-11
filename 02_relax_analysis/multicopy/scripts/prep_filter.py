import sys, os
# prep_filter.py v2.6
# NEW IN v2.6 (root-cause fix for hyphy relax failures):
#   hyphy relax hard-fails when the codon alignment contains stop codons
#   ("The input alignment must have the number of sites that is divisible
#   by 3 and must not contain stop codons"). Our input peptides end with
#   '*' and our CDS end with the terminal stop codon; after MUSCLE+trimAl
#   some '*' columns survive and map to TAA/TAG/TGA in the codon alignment.
#   Fix at the source:
#   - strip ONE terminal '*' from each peptide
#   - strip the matching terminal stop codon (TAA/TAG/TGA) from each CDS
#   - drop sequences still containing '*' afterwards (internal stop =
#     pseudogene / frameshift / misannotation)
#   - length filter uses true protein length (terminal '*' excluded)
#   - translation check runs on the stripped pair
#   - too-few families now get a skip marker (worker honours it)
#   - if only the CDS changed (terminal stop stripped, peptide identical)
#     the CDS file is rewritten without purging downstream (mapping only
#     covers 3 x AA-width columns, so downstream output is unchanged)
# Everything else identical to v2.5 (translation-consistency + MIN_AA_LEN
# filter, composition-change purge, PJA_absent skip).

workdir = os.environ.get("WORKDIR")
ortho_dir = os.environ.get("ORTHO_SEQ_DIR")
id_list = os.environ.get("ID_LIST")
min_aa = int(os.environ.get("MIN_AA_LEN", "100"))

aa_out = os.path.join(workdir, "01_aa_filtered")
cds_out = os.path.join(workdir, "01_cds_filtered")
cds_in = os.path.join(workdir, "01_cds")
os.makedirs(aa_out, exist_ok=True); os.makedirs(cds_out, exist_ok=True)

STOPS = ("TAA", "TAG", "TGA")

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

def strip_terminal_stop_cds(s):
    s = s.strip()
    if len(s) >= 3 and s[-3:].upper() in STOPS:
        return s[:-3], True
    return s, False

def strip_terminal_star(prot):
    prot = prot.strip()
    if prot.endswith("*"):
        return prot[:-1], True
    return prot, False

def read_fasta(fp):
    d = {}; name = None
    for line in open(fp, errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; d[name] = []
        elif name:
            d[name].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}

# downstream products to purge when the filtered peptide set/content changed
PURGE_PATS = [
    ("02_muscle",   ["{og}.aa.aln"]),
    ("03_trimAL",   ["{og}.aa.trim.aln"]),
    ("04_pal2nal",  ["{og}.codon.full.fasta", "{og}.codon.fasta"]),
    ("05_iqtree",   ["{og}.treefile", "{og}.log", "{og}.log.gz", "{og}.ckp", "{og}.ckp.gz",
                     "{og}.iqtree", "{og}.contree", "{og}.splits.nex", "{og}.ufboot",
                     "{og}.uuid"]),
    ("06_tree_test",["{og}.test.tree"]),
    ("07_relax_json",["{og}.RELAX.json"]),
]

def purge_downstream(workdir, og):
    for sub, pats in PURGE_PATS:
        for p in pats:
            fp = os.path.join(workdir, sub, p.format(og=og))
            if os.path.exists(fp): os.remove(fp)

def clear_markers(workdir, og, kinds=("done", "fail", "skip")):
    for m in kinds:
        mp = os.path.join(workdir, "logs", m, og)
        if os.path.exists(mp): os.remove(mp)

ogs = [l.split()[0] for l in open(id_list) if l.strip()]
log = open(os.path.join(workdir, "logs", "dropped_genes.log"), "a")
n_ok = 0; n_changed = 0; n_pja_absent = 0; n_too_few = 0
n_cds_only = 0
n_star_stripped = 0; n_stop_stripped = 0; n_internal_stop = 0
for og in ogs:
    aa_fp = os.path.join(ortho_dir, og + ".fa")
    if not os.path.exists(aa_fp): aa_fp = os.path.join(ortho_dir, og + ".fasta")
    cds_fp = os.path.join(cds_in, og + ".cds.fasta")
    if not (os.path.exists(aa_fp) and os.path.exists(cds_fp)):
        continue
    aa = read_fasta(aa_fp); cds = read_fasta(cds_fp)
    keep = []
    for gid, prot in aa.items():
        nuc = cds.get(gid)
        if nuc is None:
            log.write(f"{og}\t{gid}\tno_CDS\n"); continue
        prot2, had_star = strip_terminal_star(prot)
        nuc2, had_stop = strip_terminal_stop_cds(nuc)
        if had_star: n_star_stripped += 1
        if had_stop: n_stop_stripped += 1
        if "*" in prot2:
            n_internal_stop += 1
            log.write(f"{og}\t{gid}\tinternal_stop_residue\n"); continue
        prot_tr = tr(nuc2)
        if prot_tr != prot2:
            log.write(f"{og}\t{gid}\ttranslation_mismatch(pep={len(prot2)},cds_aa={len(prot_tr)})\n"); continue
        if len(prot2) < min_aa:
            log.write(f"{og}\t{gid}\tbelow_min_aa_len({len(prot2)}<{min_aa})\n"); continue
        keep.append((gid, prot2, nuc2))
    pja_left = [g for g, _, _ in keep if g.startswith("PJA|") or g.startswith("PJA")]
    if len(keep) < 10:
        n_too_few += 1
        log.write(f"{og}\tALL\ttoo_few_consistent_genes({len(keep)})\n")
        purge_downstream(workdir, og)
        for fp in (os.path.join(aa_out, og + ".fa"), os.path.join(cds_out, og + ".cds.fasta")):
            if os.path.exists(fp): os.remove(fp)
        clear_markers(workdir, og)
        with open(os.path.join(workdir, "logs", "skip", og), "w") as s:
            s.write(f"{og}\ttoo_few_consistent_genes({len(keep)} genes after filter)\n")
        continue
    if not pja_left:
        n_pja_absent += 1
        purge_downstream(workdir, og)
        for fp in (os.path.join(aa_out, og + ".fa"), os.path.join(cds_out, og + ".cds.fasta")):
            if os.path.exists(fp): os.remove(fp)
        clear_markers(workdir, og, ("done", "fail"))
        with open(os.path.join(workdir, "logs", "skip", og), "w") as s:
            s.write(f"{og}\tPJA_absent_after_filter(all PJA copies < {min_aa} aa)\n")
        continue
    # write only if content changed
    new_aa = "".join(f">{g}\n{p}\n" for g, p, _ in keep)
    new_cds = "".join(f">{g}\n{n}\n" for g, _, n in keep)
    aa_fp_out = os.path.join(aa_out, og + ".fa")
    cds_fp_out = os.path.join(cds_out, og + ".cds.fasta")
    old_aa = open(aa_fp_out, encoding="utf-8", errors="ignore").read() if os.path.exists(aa_fp_out) else None
    old_cds = open(cds_fp_out, encoding="utf-8", errors="ignore").read() if os.path.exists(cds_fp_out) else None
    aa_changed = (old_aa != new_aa)
    cds_changed = (old_cds != new_cds)
    if aa_changed:
        with open(aa_fp_out, "w") as o: o.write(new_aa)
        with open(cds_fp_out, "w") as o: o.write(new_cds)
        n_changed += 1
        purge_downstream(workdir, og)
        clear_markers(workdir, og)
    elif cds_changed:
        # peptide set identical; only terminal CDS stop codons were stripped.
        # codon mapping never reaches the trailing stop column, so existing
        # downstream products stay valid: rewrite CDS only.
        with open(cds_fp_out, "w") as o: o.write(new_cds)
        n_cds_only += 1
    else:
        pass  # unchanged: keep all existing downstream products
    n_ok += 1
log.close()
print(f"[step0.5v2.6] MIN_AA_LEN={min_aa}: {n_ok} families usable; "
      f"{n_changed} recomposed (downstream purged for rebuild); "
      f"{n_cds_only} cds-only restripped (downstream kept); "
      f"{n_pja_absent} excluded (no PJA left); {n_too_few} too-few-genes; "
      f"terminal '*' stripped on {n_star_stripped} peptides, "
      f"terminal stop codon stripped on {n_stop_stripped} CDS, "
      f"{n_internal_stop} internal-stop genes dropped")
