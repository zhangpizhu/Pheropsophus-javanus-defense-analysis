#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Figure S7 — RELAX analysis of selection intensity on P. javanus branches
across expanded, contracted and stable gene-family groups (multicopy and single-copy arms).

Input : data/FigS7_RELAX_input.csv  (one row per family; identical to Supplementary Table S22)
        columns: orthogroup, arm, group, k, LRT, p_value, q_BH, direction, significant
Output: FigS7_RELAX三组比较.png (300 dpi) + FigS7_RELAX三组比较.pdf (vector)

Panel A: distribution of log10(k) per group (boxplot + per-family points).
         Point color: blue  = significant relaxation   (BH q < 0.05, k < 1)
                      red   = significant intensification (BH q < 0.05, k > 1)
                      grey  = not significant
         Dashed line: k = 1 (no selection shift).
Panel B: percentage of families with a significant selection shift per group
         (stacked: relaxed + intensified). Fisher's exact tests vs the stable
         multicopy control; ** P < 0.01, ns = not significant.

Usage:
    python FigS7_RELAX_plot.py [input.csv] [output_prefix]
Defaults: input.csv      = data/FigS7_RELAX_input.csv (same folder as script)
          output_prefix  = FigS7_RELAX三组比较
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

# ---------------------------------------------------------------- input ----
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(script_dir, 'data/FigS7_RELAX_input.csv')
out_prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.join(script_dir, 'FigS7_RELAX三组比较')

df = pd.read_csv(csv_path)
df['log10k'] = np.log10(df['k'])
df['sig'] = df['significant'].astype(str).str.lower().eq('yes')
df['sig_relaxed'] = df['sig'] & (df['k'] < 1)
df['sig_intensified'] = df['sig'] & (df['k'] > 1)

# ------------------------------------------------------- group definition ----
# (key, axis label, row selector)
GROUPS = [
    ('expanded',  'Expanded',                    (df['group'] == 'expanded')),
    ('contracted', 'Contracted',                 (df['group'] == 'contracted')),
    ('stable_mc',  'Stable,\nmulticopy',         (df['group'] == 'stable') & (df['arm'] == 'multicopy')),
    ('stable_sc',  'Stable,\nsingle-copy',       (df['group'] == 'stable') & (df['arm'] == 'single-copy')),
]
REF = 'stable_mc'  # reference group for Fisher's exact tests

C_REL = '#2166ac'   # relaxed, blue
C_INT = '#b2182b'   # intensified, red
C_NS  = '#c7c7c7'   # not significant, grey

data = {key: df.loc[mask].copy() for key, _, mask in GROUPS}
labels = [f'{lab}\n(n = {len(data[key])})' for key, lab, _ in GROUPS]

# ------------------------------------------------------------ statistics ----
def fisher_vs_ref(key):
    a = int(data[key]['sig'].sum()); na = len(data[key]) - a
    b = int(data[REF]['sig'].sum()); nb = len(data[REF]) - b
    return fisher_exact([[a, na], [b, nb]])[1]

def stars(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'ns'

pvals = {key: fisher_vs_ref(key) for key, _, _ in GROUPS if key != REF}
for key, p in pvals.items():
    print(f'Fisher exact, {key} vs {REF}: P = {p:.4g}  ({stars(p)})')

# ---------------------------------------------------------------- figure ----
rng = np.random.default_rng(42)
fig = plt.figure(figsize=(10.8, 6.0084))  # aspect = 1846:1027 (identical to the embedded version)
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.22,
                      left=0.075, right=0.985, top=0.90, bottom=0.16)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])

# ---- Panel A: log10(k) distributions --------------------------------------
x = np.arange(len(GROUPS))
for i, (key, _, _) in enumerate(GROUPS):
    d = data[key].reset_index(drop=True)
    jit = rng.uniform(-0.13, 0.13, size=len(d))
    for mask, col, s in [(~d['sig'], C_NS, 9), (d['sig_relaxed'], C_REL, 26), (d['sig_intensified'], C_INT, 26)]:
        idx = np.where(mask.to_numpy())[0]
        axA.scatter(i + jit[idx], d['log10k'].to_numpy()[idx],
                    s=s, c=col, alpha=0.9 if s > 10 else 0.85, edgecolors='none', zorder=3)
bp = axA.boxplot([data[key]['log10k'] for key, _, _ in GROUPS], positions=x, widths=0.5,
                 showfliers=False, patch_artist=True,
                 boxprops=dict(facecolor='white', edgecolor='#4d4d4d', lw=1.2),
                 medianprops=dict(color='black', lw=1.6),
                 whiskerprops=dict(color='#4d4d4d', lw=1.0),
                 capprops=dict(color='#4d4d4d', lw=1.0))
axA.axhline(0, ls='--', c='#333333', lw=1.2, zorder=2)
axA.set_xticks(x)
axA.set_xticklabels(labels, fontsize=9)
axA.set_ylabel('log$_{10}$(k)', fontsize=12)
axA.grid(True, alpha=0.25, zorder=0)
axA.set_axisbelow(True)
legA = axA.legend(handles=[
    plt.Line2D([], [], marker='o', ls='', c=C_REL, label='Significant, relaxed (k < 1)'),
    plt.Line2D([], [], marker='o', ls='', c=C_INT, label='Significant, intensified (k > 1)'),
    plt.Line2D([], [], marker='o', ls='', c=C_NS, label='Not significant')],
    loc='lower left', fontsize=8.5, frameon=False)

# ---- Panel B: % families with significant shift ---------------------------
totals, rel_pct, int_pct = [], [], []
for key, _, _ in GROUPS:
    d = data[key]
    totals.append(100.0 * d['sig'].sum() / len(d))
    rel_pct.append(100.0 * d['sig_relaxed'].sum() / len(d))
    int_pct.append(100.0 * d['sig_intensified'].sum() / len(d))

axB.bar(x, rel_pct, width=0.62, color=C_REL, label='Relaxed (k < 1)', zorder=3)
axB.bar(x, int_pct, width=0.62, bottom=rel_pct, color=C_INT, label='Intensified (k > 1)', zorder=3)
for i, t in enumerate(totals):
    axB.text(i, t + 1.3, f'{t:.0f}%', ha='center', fontsize=11)
    key = GROUPS[i][0]
    if key != REF:
        axB.text(i, t + 5.6, stars(pvals[key]), ha='center', fontsize=10, color='#333333')
axB.set_xticks(x)
axB.set_xticklabels(labels, fontsize=9)
axB.set_ylabel('Families with significant\nselection shift (%)', fontsize=12)
axB.set_ylim(0, 52)
axB.grid(True, alpha=0.25, zorder=0)
axB.set_axisbelow(True)
axB.legend(loc='upper right', fontsize=9, frameon=False)

# ---- panel letters ----------------------------------------------------------
for ax, letter in [(axA, 'A'), (axB, 'B')]:
    ax.text(-0.10, 1.04, letter, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')

fig.savefig(out_prefix + '.png', dpi=300)
fig.savefig(out_prefix + '.pdf')
print('saved:', out_prefix + '.png /', out_prefix + '.pdf')
