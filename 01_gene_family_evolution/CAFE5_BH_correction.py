#!/usr/bin/env python3
"""
CAFE5_BH_correction.py
======================
Multiple-testing correction (Benjamini-Hochberg FDR) of CAFE5 gene-family
evolution results for the Pheropsophus javanus genome study.

Background
----------
CAFE5 reports a nominal family-wide P value per gene family
(Base_family_results.txt) and per-branch P values for families that are
significant at the family-wide level (Base_branch_probabilities.tab).
When ~10,000 families are tested, nominal P < 0.05 implies hundreds of
expected false positives, so we apply the Benjamini-Hochberg (BH)
procedure (Benjamini & Hochberg, 1995, JRSS B 57:289-300) at two levels:

  1. Family level  : BH across ALL tested families (family-wide P values).
  2. Branch level  : BH across the per-branch P values of the focal branch
                     (here PJA<13> = Pheropsophus javanus), restricted to
                     families significant at the nominal family-wide level
                     (the set for which CAFE5 reports branch P values).

A conservative subset is defined as families passing BOTH corrections
(q < 0.05 at family level AND q < 0.05 at branch level).

Inputs (standard CAFE5 Base-model outputs)
------------------------------------------
  Base_family_results.txt        FamilyID <tab> pvalue <tab> Significant
  Base_branch_probabilities.tab  FamilyID <tab> one P value per tree node
  Base_change.tab                FamilyID <tab> one signed change per node

Usage
-----
  python3 CAFE5_BH_correction.py \
      --family_results Base_family_results.txt \
      --branch_probs   Base_branch_probabilities.tab \
      --change         Base_change.tab \
      --branch         'PJA<13>' \
      --q              0.05 \
      --outdir         cafe_bh_out

Outputs (written to --outdir)
-----------------------------
  all_families_BH_qvalues.tsv    per-family P, family-level q, branch-level q,
                                 signed change on the focal branch
  familyBH_expanded_IDs.txt      families with family q < Q, change > 0
  familyBH_contracted_IDs.txt    families with family q < Q, change < 0
  branchBH_expanded_IDs.txt      nominal-sig families with branch q < Q, change > 0
  branchBH_contracted_IDs.txt    nominal-sig families with branch q < Q, change < 0
  strict_expanded_IDs.txt        both corrections significant, change > 0
  strict_contracted_IDs.txt      both corrections significant, change < 0
  summary.txt                    counts under each significance criterion

Dependencies: Python >= 3.8, numpy, pandas.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def bh_qvalues(pvals):
    """Benjamini-Hochberg adjusted q values (BH step-up procedure).

    Parameters
    ----------
    pvals : array-like of float, raw P values (NaN not allowed).

    Returns
    -------
    numpy array of q values aligned with the input order.
    """
    p = np.asarray(pvals, dtype=float)
    if np.isnan(p).any():
        raise ValueError("bh_qvalues: input contains NaN")
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    # enforce monotonicity: q(i) = min over j>=i of raw q(j)
    q_ordered = q[order]
    q_ordered = np.minimum.accumulate(q_ordered[::-1])[::-1]
    q[order] = np.minimum(q_ordered, 1.0)
    return q


def load_inputs(family_results, branch_probs, change, branch):
    fam = pd.read_csv(family_results, sep="\t")
    fam = fam.rename(columns={fam.columns[0]: "FamilyID"})

    br = pd.read_csv(branch_probs, sep="\t")
    br = br.rename(columns={br.columns[0]: "FamilyID"})
    if branch not in br.columns:
        sys.exit(f"ERROR: branch '{branch}' not found in {branch_probs}. "
                 f"Available: {list(br.columns)}")

    chg = pd.read_csv(change, sep="\t")
    chg = chg.rename(columns={chg.columns[0]: "FamilyID"})
    if branch not in chg.columns:
        sys.exit(f"ERROR: branch '{branch}' not found in {change}.")
    # CAFE5 writes changes with explicit '+' sign, e.g. '+34' -> parse to int
    chg["focal_change"] = pd.to_numeric(
        chg[branch].astype(str).str.replace("+", "", regex=False))

    return fam, br[["FamilyID", branch]], chg[["FamilyID", "focal_change"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_results", required=True)
    ap.add_argument("--branch_probs", required=True)
    ap.add_argument("--change", required=True)
    ap.add_argument("--branch", default="PJA<13>",
                    help="column name of the focal branch (default: PJA<13>)")
    ap.add_argument("--q", type=float, default=0.05, help="FDR threshold (default 0.05)")
    ap.add_argument("--outdir", default="cafe_bh_out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fam, br, chg = load_inputs(args.family_results, args.branch_probs,
                               args.change, args.branch)

    # ---- level 1: family-wide BH across all tested families ----
    fam["family_qvalue_BH"] = bh_qvalues(fam["pvalue"].values)

    # ---- level 2: branch-level BH within nominally significant families ----
    sig_nominal = fam.loc[fam["pvalue"] < args.q, "FamilyID"]
    br_sig = br[br["FamilyID"].isin(sig_nominal)].copy()
    br_sig["branch_qvalue_BH"] = bh_qvalues(br_sig[args.branch].values)

    out = (fam[["FamilyID", "pvalue", "family_qvalue_BH"]]
           .merge(br_sig[["FamilyID", "branch_qvalue_BH"]], on="FamilyID", how="left")
           .merge(chg, on="FamilyID"))
    out = out.rename(columns={"pvalue": "family_pvalue"})
    out.to_csv(os.path.join(args.outdir, "all_families_BH_qvalues.tsv"),
               sep="\t", index=False)

    crit = {
        "familyBH": out["family_qvalue_BH"] < args.q,
        "branchBH": (out["family_pvalue"] < args.q) & (out["branch_qvalue_BH"] < args.q),
        "strict":   (out["family_qvalue_BH"] < args.q) & (out["branch_qvalue_BH"] < args.q),
    }

    lines = [f"Focal branch: {args.branch}", f"FDR threshold: {args.q}",
             f"Families tested (family level): {len(out)}",
             f"Nominally significant families (branch-level test set): {len(br_sig)}", ""]
    for name, mask in crit.items():
        sel = out[mask]
        exp = sel.loc[sel["focal_change"] > 0, "FamilyID"]
        con = sel.loc[sel["focal_change"] < 0, "FamilyID"]
        exp.to_csv(os.path.join(args.outdir, f"{name}_expanded_IDs.txt"),
                   index=False, header=False)
        con.to_csv(os.path.join(args.outdir, f"{name}_contracted_IDs.txt"),
                   index=False, header=False)
        lines.append(f"{name}: {len(sel)} significant "
                     f"({len(exp)} expanded, {len(con)} contracted on focal branch)")

    summary = "\n".join(lines) + "\n"
    with open(os.path.join(args.outdir, "summary.txt"), "w") as fh:
        fh.write(summary)
    sys.stdout.write(summary)


if __name__ == "__main__":
    main()
