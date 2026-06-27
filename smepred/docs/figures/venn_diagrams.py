#!/usr/bin/env python3
"""Regenerate data composition and feature engineering figure (data_composition.png)."""

import os, sys
from pathlib import Path
os.chdir(Path(__file__).parent.parent.parent)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

FIGS_DIR = Path("docs/figures")
FIGS_DIR.mkdir(exist_ok=True)

C_BLUE   = '#1f77b4'
C_RED    = '#d62728'
C_GREEN  = '#2ca02c'
C_PURPLE = '#9467bd'
C_ORANGE = '#ff7f0e'
C_TEAL   = '#17becf'

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 0.8, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.2))

labels = ['Position-Aware\n55,730 (66.7%)', 'Hetero Patent\n23,187 (27.8%)', 'CMsiRNAdb\n4,618 (5.5%)']
sizes = [55730, 23187, 4618]
colors_pie = [C_BLUE, C_ORANGE, C_GREEN]
wedges, texts = ax1.pie(sizes, labels=labels, colors=colors_pie, startangle=90, textprops={'fontsize': 7})
for w in wedges:
    w.set_edgecolor('white'); w.set_linewidth(1.5)
ax1.set_title('a) Training Data Composition (83,535 rows)', fontsize=10, fontweight='bold', pad=8)

categories = ['Per-Position\nFlags (33x42)', 'Global\nCounts (31x2)', 'Summary\nStats (9x2)', 'Log\nConc.']
dims = [1386, 62, 18, 1]
bars = ax2.bar(categories, dims, color=[C_BLUE, C_ORANGE, C_GREEN, C_TEAL], width=0.55, edgecolor='black', linewidth=0.8)
for bar, d in zip(bars, dims):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 15, f'{d}', ha='center', fontsize=8, fontweight='bold')
ax2.set_ylabel('Feature Dimensions')
ax2.set_title('b) Feature Vector Composition (1,467 total)', fontsize=10, fontweight='bold', pad=8)
ax2.set_ylim(0, 1550)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

fig.suptitle('Figure 2: Data Curation and Feature Engineering', fontsize=11, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig(FIGS_DIR / 'data_composition.png')
plt.close(fig)
print('  venn_diagrams.py -> data_composition.png')
