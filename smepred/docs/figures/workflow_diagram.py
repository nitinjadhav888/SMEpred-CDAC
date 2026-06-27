#!/usr/bin/env python3
"""Regenerate pipeline workflow diagram (workflow.png)."""

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
C_GRAY   = '#7f7f7f'

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 0.8, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})

styles = {
    'input':  dict(boxstyle='round,pad=0.35', facecolor='#e8f5e9', edgecolor=C_GREEN, linewidth=1.5),
    'step':   dict(boxstyle='round,pad=0.35', facecolor='#e3f2fd', edgecolor=C_BLUE, linewidth=1.5),
    'output': dict(boxstyle='round,pad=0.35', facecolor='#fce4ec', edgecolor=C_RED, linewidth=1.5),
    'data':   dict(boxstyle='round,pad=0.3',  facecolor='#fff3e0', edgecolor=C_ORANGE, linewidth=1.5),
}

def rbox(ax, x, y, w, h, text, style, fontsize=8, color='k', ha='center', va='center'):
    ax.text(x+w/2, y+h/2, text, ha=ha, va=va, fontsize=fontsize, color=color, bbox=style, zorder=5)

def rarrow(ax, x1, y1, x2, y2, color='#888', lw=1.2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='->', color=color, lw=lw), zorder=3)

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
ax.set_xlim(0, 10); ax.set_ylim(0, 3.5)
ax.axis('off')

# Step 1: Input
rbox(ax, 0.1, 2.5, 2.0, 0.7, 'Step 1: Input\n21-nt sense + antisense\nor mRNA sequence', styles['input'], fontsize=7)
rarrow(ax, 2.1, 2.85, 3.0, 2.85)

# Step 2: Feature extraction
rbox(ax, 3.0, 2.5, 2.0, 0.7, 'Step 2: Feature Extraction\n214-d (naked) / 1,467-d (mod)\nOne-hot + TNC + positional', styles['step'], fontsize=7)
rarrow(ax, 5.0, 2.85, 5.8, 2.85)

# Step 3: Scoring
rbox(ax, 5.8, 2.5, 1.8, 0.7, 'Step 3: Model Scoring\nLightGBM inference\nRaw efficacy (0-100)', styles['step'], fontsize=7)
rarrow(ax, 7.6, 2.85, 8.3, 2.85)

# Step 4: Output
rbox(ax, 8.3, 2.5, 1.5, 0.7, 'Step 4: Results\nRanked predictions\nExport JSON', styles['output'], fontsize=7)

# Bottom row - detail
# Detail boxes for Single-mod scan
rbox(ax, 1.8, 0.5, 2.6, 0.7, 'Single-Mod Scan\n1,302 variants (31 mods x 42 positions)\nEnumerates all single substitutions', styles['step'], fontsize=7)
rarrow(ax, 3.1, 1.2, 3.3, 1.9)

# Multi-mod beam search
rbox(ax, 5.2, 0.5, 2.6, 0.7, 'Multi-Mod Beam Search\nBeam width = 30, up to 14 rounds\nPlateau-based early stopping, ~20s', styles['step'], fontsize=7)
rarrow(ax, 6.5, 1.2, 6.7, 1.9)

# Biophysical validation
rbox(ax, 7.8, 0.5, 2.0, 0.7, 'Adjusted Score\nRaw - 0.70 x Penalty\n5 penalty domains', styles['data'], fontsize=7)

# Arrow from step 3 to scan
rarrow(ax, 6.7, 2.5, 6.7, 2.1)
rarrow(ax, 3.1, 2.5, 3.1, 2.1)
rarrow(ax, 3.1, 1.9, 6.5, 1.2)

# Arrow from beam search to biophysics
rarrow(ax, 7.8, 1.2, 7.8, 1.9)

ax.text(0.4, 3.3, 'Figure 7: HelixZero-CMS Prediction Workflow', fontsize=11, fontweight='bold')
fig.savefig(FIGS_DIR / 'workflow.png')
plt.close(fig)
print('  workflow_diagram.py -> workflow.png')
