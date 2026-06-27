#!/usr/bin/env python3
"""Regenerate system architecture diagram (architecture.png)."""

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

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 0.8, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})

styles = {
    'ui':   dict(boxstyle='round,pad=0.3', facecolor='#e3f2fd', edgecolor=C_BLUE, linewidth=1.5),
    'api':  dict(boxstyle='round,pad=0.3', facecolor='#fff3e0', edgecolor=C_ORANGE, linewidth=1.5),
    'core': dict(boxstyle='round,pad=0.3', facecolor='#e8f5e9', edgecolor=C_GREEN, linewidth=1.5),
    'data': dict(boxstyle='round,pad=0.3', facecolor='#f3e5f5', edgecolor=C_PURPLE, linewidth=1.5),
    'out':  dict(boxstyle='round,pad=0.3', facecolor='#fce4ec', edgecolor=C_RED, linewidth=1.5),
}

def rbox(ax, x, y, w, h, text, style, fontsize=8, color='k', ha='center', va='center'):
    ax.text(x+w/2, y+h/2, text, ha=ha, va=va, fontsize=fontsize, color=color, bbox=style, zorder=5)

def rarrow(ax, x1, y1, x2, y2, color='#888', lw=1.2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='->', color=color, lw=lw), zorder=3)

fig, ax = plt.subplots(1, 1, figsize=(8.5, 6))
ax.set_xlim(0, 10); ax.set_ylim(0, 7.5)
ax.axis('off')

ax.text(0.2, 7.0, 'Interface Layer', fontsize=8, fontweight='bold', color=C_BLUE)
ax.text(0.2, 5.3, 'API Layer', fontsize=8, fontweight='bold', color=C_ORANGE)
ax.text(0.2, 3.7, 'Core Engine', fontsize=8, fontweight='bold', color=C_GREEN)
ax.text(0.2, 1.9, 'Data & Models', fontsize=8, fontweight='bold', color=C_PURPLE)

rbox(ax, 0.8, 6.2, 2.5, 0.6, 'Web UI (app.html)\nSingle-Page Application', styles['ui'], fontsize=7)
rbox(ax, 4.0, 6.2, 2.0, 0.6, 'CLI (cli/run.py)\nClick Commands', styles['ui'], fontsize=7)
rbox(ax, 6.7, 6.2, 2.5, 0.6, 'REST API (api/main.py)\nFastAPI + Swagger', styles['api'], fontsize=7)

rbox(ax, 0.3, 4.8, 2.2, 0.5, '/rank - Gene Scan\n/single-mod - Enumeration\n/multi-mod - Beam Search\n'
     '/multi-mod-scan - Batch\n/rank - Score & Rank', styles['api'], fontsize=6)
rbox(ax, 3.2, 4.8, 2.2, 0.5, 'Prediction Orchestrator\npredictor.py', styles['core'], fontsize=7)
rbox(ax, 6.1, 4.8, 2.4, 0.5, 'Modification Engine\nmodification_engine.py\nBeam Search + Scoring', styles['core'], fontsize=7)
rbox(ax, 9.1, 4.8, 0.7, 0.5, 'Parser\nparser.py', styles['core'], fontsize=6)

rbox(ax, 0.5, 3.0, 2.0, 0.5, 'Feature Extraction\nfeatures.py\n214-d / 1,467-d', styles['core'], fontsize=7)
rbox(ax, 3.2, 3.0, 2.0, 0.5, 'LightGBM Model B\nmodel_b.pkl\n1,115 trees', styles['data'], fontsize=7)
rbox(ax, 5.9, 3.0, 2.0, 0.5, 'Biophysical Penalty\nbiophysics.py\n5 domains', styles['data'], fontsize=7)
rbox(ax, 8.6, 3.0, 1.2, 0.5, 'Filters\nfilters.py', styles['data'], fontsize=7)

rbox(ax, 0.5, 1.2, 2.0, 0.5, 'Training Data\n83,535 rows\n3 sources', styles['data'], fontsize=7)
rbox(ax, 3.2, 1.2, 2.0, 0.5, 'Naked Model V4\ncalibrator_naked.pkl\n214-d features', styles['data'], fontsize=7)
rbox(ax, 5.9, 1.2, 2.0, 0.5, 'Mod Codes JSON\nmodification_codes.json\n31 symbols', styles['data'], fontsize=7)
rbox(ax, 8.6, 1.2, 1.2, 0.5, 'Metadata\nmodel_b_meta.json', styles['data'], fontsize=7)

rbox(ax, 7.0, 6.2, 0.7, 0.6, 'JSON\nOutput', styles['out'], fontsize=7, color=C_RED)

rarrow(ax, 3.3, 6.2, 4.0, 5.5)
rarrow(ax, 6.0, 6.2, 6.7, 5.5)
rarrow(ax, 2.0, 4.8, 2.0, 3.6)
rarrow(ax, 4.3, 4.8, 4.2, 3.6)
rarrow(ax, 7.3, 4.8, 6.9, 3.6)
rarrow(ax, 1.5, 3.0, 1.5, 1.8)
rarrow(ax, 4.2, 3.0, 4.2, 1.8)
rarrow(ax, 6.9, 3.0, 6.9, 1.8)
rarrow(ax, 9.2, 4.8, 9.2, 3.6)
rarrow(ax, 9.2, 3.6, 7.9, 3.6)
rarrow(ax, 7.9, 3.6, 7.5, 5.5)
rarrow(ax, 9.2, 1.8, 9.2, 1.2)

ax.text(0.5, 7.3, 'Figure 1: HelixZero-CMS System Architecture', fontsize=11, fontweight='bold')
fig.savefig(FIGS_DIR / 'architecture.png')
plt.close(fig)
print('  schematic_diagram.py -> architecture.png')
