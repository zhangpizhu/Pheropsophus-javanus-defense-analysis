#!/usr/bin/env python3
"""
selection_BH_correction.py
==========================
Benjamini-Hochberg correction of PAML selection-test results for the
Pheropsophus javanus genome study (INZ-2026-07-520).

Workflow (exactly as used to produce Supplementary Table S20 from
Supplementary Tables S18 and S19)
------------------------------------------------------------------
1. PAML codeml branch-site tests (PSG arm: Model A vs. null) and
   branch-model tests (REG arm: two-ratio vs. one-ratio) were run for all
   2,568 single-copy orthologs, with P. javanus as the foreground branch.
2. The nominal-significant lists (LRT P <= 0.05; Tables S18 and S19) contain,
   per gene, the likelihood-ratio statistic. The LRT P value is
   P = chi2_sf(LRT, df = 1).
3. Because every ortholog absent from these lists has a larger P value than
   every listed gene, the ranks of the listed genes within the full set of
   2,568 tested orthologs are simply their ranks within the list. The
   Benjamini-Hochberg q values are therefore computed exactly as
       q_i = min_{j >= i} [ P_j * N_tested / j ],   N_tested = 2568
   (Benjamini & Hochberg, 1995, JRSS B 57:289-300).
4. Genes with q < 0.05 form the corrected primary sets
   (147 PSGs and 183 REGs in the manuscript; Table S20).

Input
-----
A tab- or comma-separated table exported from Table S18 or Table S19, with a
header line and columns:

  Gene_ID, Likelihood_Ratio, Gene_name, Significance

(only Gene_ID and Likelihood_Ratio are used; extra columns are ignored).
The file should contain ONLY the nominal-significant genes of ONE test arm;
the total number of tested orthologs is passed via --n_tested.

Usage
-----
  python3 selection_BH_correction.py --input Table_S18_PSGs.tsv \
      --arm PSG --n_tested 2568 --q 0.05 --output Table_S18_BH.tsv

Output
------
A TSV with columns: Gene_ID, Gene_name, LRT, LRT_P, BH_q, significant
(yes/no at q < --q), plus a printed summary line with the counts.

Dependencies: Python >= 3.8 (standard library only).
"""

import argparse
import csv
import math
import sys


def chi2_sf_df1(x):
    """Survival function of the chi-square distribution with df = 1."""
    return math.erfc(math.sqrt(x / 2.0))


def sniff_delimiter(line):
    return '\t' if '\t' in line else ','


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--input', required=True,
                    help='Nominal-significant gene table (one test arm)')
    ap.add_argument('--arm', default='PSG', choices=['PSG', 'REG'],
                    help='Test arm label, used in output naming/summary')
    ap.add_argument('--n_tested', type=int, required=True,
                    help='Total orthologs tested in this arm (2568 in this study)')
    ap.add_argument('--q', type=float, default=0.05, help='FDR threshold (default 0.05)')
    ap.add_argument('--output', required=True, help='Output TSV path')
    args = ap.parse_args()

    with open(args.input, newline='') as fh:
        first = fh.readline()
        delim = sniff_delimiter(first)
        fh.seek(0)
        rdr = csv.reader(fh, delimiter=delim)
        rows = [r for r in rdr if r and r[0].strip()]
    # drop header
    if rows and rows[0][0].strip().lower().startswith('gene'):
        rows = rows[1:]

    recs = []
    for r in rows:
        gid = r[0].strip()
        lrt = float(r[1])
        name = r[2].strip() if len(r) > 2 else ''
        recs.append((gid, name, lrt))
    recs.sort(key=lambda x: -x[2])

    m = len(recs)
    if m > args.n_tested:
        sys.exit(f'ERROR: {m} input genes exceed --n_tested {args.n_tested}')

    qvals = [0.0] * m
    prev = 1.0
    for i in range(m - 1, -1, -1):
        p = chi2_sf_df1(recs[i][2])
        q = min(p * args.n_tested / (i + 1), prev, 1.0)
        prev = q
        qvals[i] = q

    n_sig = sum(1 for q in qvals if q < args.q)
    with open(args.output, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['Gene_ID', 'Gene_name', 'LRT', 'LRT_P', 'BH_q', 'significant'])
        for (gid, name, lrt), q in zip(recs, qvals):
            p = chi2_sf_df1(lrt)
            w.writerow([gid, name, f'{lrt:.6g}', f'{p:.6g}', f'{q:.6g}',
                        'yes' if q < args.q else 'no'])

    print(f'{args.arm}: nominal genes = {m}, N_tested = {args.n_tested}, '
          f'BH q < {args.q} -> {n_sig} genes')
    print(f'written: {args.output}')


if __name__ == '__main__':
    main()
