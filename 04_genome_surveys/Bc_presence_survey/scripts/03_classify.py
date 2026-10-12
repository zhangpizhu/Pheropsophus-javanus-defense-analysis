#!/usr/bin/env python3
"""
Step 03 - Filter tBLASTn hits, run the reciprocal-best-hit check, estimate
          B. crepitans copy numbers, and produce the family classification
          (Pangdun's categories I/II/III) reported in Supplementary Table S36, Part D.

Outputs:
  results/partD_table.csv      per-query survey table (all 27 queries)
  results/family_summary.csv   per-family copy numbers and category assignment
"""
import os
import pandas as pd

RES = os.environ.get("A1_RES_DIR", "/mnt/agents/output/A1_Bc_presence_survey/results")
os.makedirs(RES, exist_ok=True)

# per-family PJA gene sets (must match 01_extract_and_prepare.py)
TARGET_FAMILIES = {
    "OG0000018": ["Chr1.474", "Chr12.128", "Chr13.61", "Chr13.66", "Chr4.705", "Chr4.706", "Unplaced6.11"],
    "OG0000048": ["Chr12.542", "Chr12.699", "Chr12.700", "Chr5.537", "Chr5.538", "Chr5.539", "Chr5.540"],
    "OG0000494": ["Chr5.469"],
    "OG0000273": ["Chr16.191", "Chr9.96", "Chr9.97"],
    "OG0000230": ["Chr14.436", "Chr14.438", "Chr14.440", "Chr4.129"],
    "OG0000272": ["Chr1.58", "Chr1.800", "Chr4.452", "Chr4.697", "Chr7.48"],
}
CYP_FAMILIES = {"OG0000230", "OG0000272"}  # superfamily-fragmented: RBH only to superfamily level

# other-beetle copy-number ranges (from Orthogroups.GeneCount.txt)
OTHER_RANGE = {
    "OG0000018": "26-37", "OG0000048": "16-28", "OG0000494": "1-21",
    "OG0000273": "2-33", "OG0000230": "2-14", "OG0000272": "2-13",
}

hits = pd.read_csv(f"{RES}/tblastn_pja_vs_Bc.tsv", sep="\t", header=None,
                   names=["q", "s", "pident", "alen", "qlen", "slen", "evalue", "bits"])
hits["og"] = hits["q"].str.split("|").str[1]
hits["qgene"] = hits["q"].str.split("|").str[0]

# query lengths from the query fasta (self-contained, also covers no-hit queries)
qlen_map, _cur = {}, None
for line in open(f"{os.environ.get('A1_DATA_DIR', '/mnt/agents/output/A1_Bc_presence_survey/data')}/pja_27_queries.faa"):
    line = line.strip()
    if line.startswith(">"):
        _cur = line[1:].split("|")[0]
        qlen_map[_cur] = 0
    else:
        qlen_map[_cur] += len(line)

# best hit per query (by bit score) — used for the identity ranges quoted in the paper
best_per_query = hits.sort_values("bits", ascending=False).drop_duplicates("qgene")

rbh = pd.read_csv(f"{RES}/blastx_Bc_vs_PJA.tsv", sep="\t", header=None,
                  names=["s", "pja", "pident", "alen", "evalue", "bits"])
rbh_best = rbh.sort_values("bits", ascending=False).drop_duplicates("s")
rbh_map = {r["s"]: (r["pja"].replace("PJA|", ""), r["pident"])
           for _, r in rbh_best.iterrows()}

# ---- per-query table (Supplementary Table S36, Part D) ----------------------
rows = []
for og, genes in TARGET_FAMILIES.items():
    for g in genes:
        hh = hits[hits["qgene"] == g]
        row = {"gene": "evm.model." + g, "orthogroup": og, "query_aa": qlen_map[g]}
        for src, tag in [("BcGL", "best_hit_gland_S1"), ("BcCB", "best_hit_combined_S2")]:
            hs = hh[hh["s"].str.startswith(src)]
            if len(hs):
                b = hs.sort_values("bits", ascending=False).iloc[0]
                row[tag] = b["s"].split("|")[1]
                row[tag + "_identity"] = f"{b['pident']:.1f}%"
                row[tag + "_aln_aa"] = int(b["alen"])
            else:
                row[tag], row[tag + "_identity"], row[tag + "_aln_aa"] = "no hit", "-", "-"
        best_overall = hh.sort_values("bits", ascending=False).iloc[0]["s"] if len(hh) else None
        if best_overall and best_overall in rbh_map:
            g2, pid = rbh_map[best_overall]
            row["reciprocal_best_hit"] = f"{g2} ({pid:.0f}%)"
        else:
            row["reciprocal_best_hit"] = "-"
        rows.append(row)
partD = pd.DataFrame(rows)
partD.to_csv(f"{RES}/partD_table.csv", index=False)

# ---- family summary ----------------------------------------------------------
# Bc copy estimate: unique contigs with HSP covering >=25% of query OR >=40% identity.
# Not applicable to CYP families (superfamily cross-hits).
filt = hits[(hits["alen"] / hits["qlen"] >= 0.25) | (hits["pident"] >= 40)]
bestc = filt.sort_values("bits", ascending=False).drop_duplicates(["og", "s"])
summary = []
for og, genes in TARGET_FAMILIES.items():
    cc = bestc[bestc["og"] == og]
    n_hit_queries = hits[hits["og"] == og]["qgene"].nunique()
    bpq = best_per_query[best_per_query["og"] == og]
    s = {"orthogroup": og, "PJA_copies": len(genes),
         "other_beetles": OTHER_RANGE[og],
         "queries_with_Bc_hit": f"{n_hit_queries}/{len(genes)}",
         "identity_range": f"{bpq['pident'].min():.0f}-{bpq['pident'].max():.0f}%"}
    if og in CYP_FAMILIES:
        s["Bc_copy_estimate"] = "n.a. (CYP superfamily fragmented across 46 orthogroups)"
    else:
        s["Bc_copy_estimate"] = f"~{len(cc)} contigs"
    # category: I = Bc homolog present; III = no-hit short partial PJA model
    s["category"] = "I (all members with hits)" if n_hit_queries == len(genes) else \
                    f"I (homologs present; {len(genes)-n_hit_queries} short partial model(s) without hit)"
    summary.append(s)
pd.DataFrame(summary).to_csv(f"{RES}/family_summary.csv", index=False)

print(partD.to_string(index=False))
print("\n--- family summary ---")
print(pd.DataFrame(summary).to_string(index=False))
