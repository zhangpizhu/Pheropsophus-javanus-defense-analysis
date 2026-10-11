#!/usr/bin/env python3
"""
Step 01 - Extract B. crepitans assemblies from Vogel et al. (2025) Supplemental Files
          and build the P. javanus query protein set for the six contracted families.

Inputs (paths configurable below):
  - Vogel_et_al_Supplemental-File-S1_Transcriptome-GLANDS_Annotation-TPM-NOIseq.xlsx
  - Vogel_et_al_Supplemental-File-S2_Transcriptome-COMBINED_Annotation-TPM-NOIseq.xlsx
  - PJA_PEP.fas  (27,510-protein P. javanus proteome; headers ">PJA|<gene>")

Outputs:
  data/Bc_glands.fna        Bc gland part-assembly contigs (BcGL|...)
  data/Bc_combined.fna      Bc combined-tissue assembly contigs (BcCB|...)
  data/pja_27_queries.faa   27 query proteins (target family members, ">PJA|<gene>|<OG>")
"""
import os
import pandas as pd
import re

S1_XLSX = "/mnt/agents/temp/Vogel_et_al_Supplemental-File-S1_Transcriptome-GLANDS_Annotation-TPM-NOIseq.xlsx"
S2_XLSX = "/mnt/agents/temp/Vogel_et_al_Supplemental-File-S2_Transcriptome-COMBINED_Annotation-TPM-NOIseq.xlsx"
PJA_PEP = "/mnt/agents/temp/PJA_PEP.fas"
OUT_DIR = os.environ.get("A1_DATA_DIR", "/mnt/agents/output/A1_Bc_presence_survey/data")
os.makedirs(OUT_DIR, exist_ok=True)

# Target families: the three contracted cuticle-protein families, the contracted
# steroid-associated sulfotransferase family, and the two contracted CYP families.
TARGET_FAMILIES = {
    "OG0000018": ["Chr1.474", "Chr12.128", "Chr13.61", "Chr13.66", "Chr4.705", "Chr4.706", "Unplaced6.11"],
    "OG0000048": ["Chr12.542", "Chr12.699", "Chr12.700", "Chr5.537", "Chr5.538", "Chr5.539", "Chr5.540"],
    "OG0000494": ["Chr5.469"],
    "OG0000273": ["Chr16.191", "Chr9.96", "Chr9.97"],
    "OG0000230": ["Chr14.436", "Chr14.438", "Chr14.440", "Chr4.129"],
    "OG0000272": ["Chr1.58", "Chr1.800", "Chr4.452", "Chr4.697", "Chr7.48"],
}

MIN_CONTIG_LEN = 100  # nt; discard very short contigs


def extract_bc(path, sheet_keyword, out, source_tag):
    xl = pd.ExcelFile(path)
    sheet = [s for s in xl.sheet_names if sheet_keyword in s][0]
    df = xl.parse(sheet)
    name_col, seq_col = df.columns[0], df.columns[2]
    n = 0
    with open(out, "w") as fh:
        for _, r in df.iterrows():
            sq = re.sub(r"[^ACGTUNacgtun]", "", str(r[seq_col]))
            if len(sq) >= MIN_CONTIG_LEN:
                fh.write(f">{source_tag}|{r[name_col]}\n{sq.upper()}\n")
                n += 1
    print(f"{path}\n  sheet={sheet}: {n} contigs -> {out}")


def load_proteome(path):
    prot, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0].replace("PJA|", "")
            prot[cur] = ""
        else:
            prot[cur] += line
    return prot


if __name__ == "__main__":
    extract_bc(S1_XLSX, "GLAND", f"{OUT_DIR}/Bc_glands.fna", "BcGL")
    extract_bc(S2_XLSX, "COMB", f"{OUT_DIR}/Bc_combined.fna", "BcCB")

    prot = load_proteome(PJA_PEP)
    n = 0
    with open(f"{OUT_DIR}/pja_27_queries.faa", "w") as fh:
        for og, genes in TARGET_FAMILIES.items():
            for g in genes:
                fh.write(f">{g}|{og}\n{prot[g]}\n")
                n += 1
    print(f"{n} query proteins -> {OUT_DIR}/pja_27_queries.faa")
