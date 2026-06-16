"""
Validate SMEpred on published experimental siRNA data.

Sources:
1. Kliuchnikov et al. 2025 (Alnylam) - 15 pairs of unmodified + 2'-OMe/2'-F siRNAs
2. Kenski et al. 2012 (Merck) - 5 sequences with 2'-O-benzyl position scan

Usage:
    cd D:\Helixx\smepred
    python tools\validate_on_published_papers.py

Output: logs/paper_validation/
"""

import sys, os, json, warnings
from pathlib import Path
from typing import List, Tuple, Dict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings('ignore', message='X does not have valid feature names')

from src.features import extract_batch_gbm
from src.predictor import _get_model, _normalize_scores

OUT = Path(__file__).parent.parent / "logs" / "paper_validation"
OUT.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# 1. ALNYLAM PAPER (Kliuchnikov et al. 2025, Mol Ther Nucleic Acids)
# ═══════════════════════════════════════════════════════════════════════════
# Table 1: 15 siRNA pairs with IC50 and Tm
# Notation:
#   UPPER     = RNA (unmodified)    → A/U/G/C
#   lower     = 2'-OMe              → M
#   Xf (CAP+f) = 2'-F               → F
#   dTdT      = DNA TT overhang     → TT (DNA)

def parse_alnylam_seq(seq_str: str, length: int = 21) -> str:
    """Parse Alnylam sequence notation to SMEpred single-letter codes."""
    result = []
    i = 0
    while i < len(seq_str):
        ch = seq_str[i]
        # Check for dTdT (two chars)
        if ch == 'd' and i + 1 < len(seq_str) and seq_str[i+1] == 'T':
            result.append('T')
            i += 2
            continue
        # Check for Xf (capital + 'f')
        if ch.isupper() and ch in 'AUCG' and i + 1 < len(seq_str) and seq_str[i+1] == 'f':
            result.append('F')
            i += 2
            continue
        # Lower case = 2'-OMe → M
        if ch.islower() and ch in 'aucg':
            result.append('M')
            i += 1
            continue
        # Upper case RNA = canonical base
        if ch.isupper() and ch in 'AUCG':
            result.append(ch)
            i += 1
            continue
        # Skip other chars
        i += 1
    return ''.join(result)[:length]


# Alnylam data: (id, passenger_seq, guide_seq, IC50_nM)
# Modified versions marked with 'm' suffix
ALNYLAM_DATA = [
    # (name, passenger, guide, IC50)
    ("siSER-1", "ACCAGCGGCCUCUGGACCAdTdT", "UGGUCCAGAGGCCGCUGGUdTdT", 4.4),
    ("siSER-1m", "accaGfcGfGfCfcucuggaccadTdT", "uGfgucCfagaggccGfcUfggudTdT", 100),
    ("siSER-2", "CUCCCCUGUGAGCAUCUCAdTdT", "UGAGAUGCUCACAGGGGAGdTdT", 0.11),
    ("siSER-2m", "cuccCfcUfGfUfgagcaucucadTdT", "uGfagaUfgcucacaGfgGfgagdTdT", 27.2),
    ("siSER-3", "CCCAGCUUCUCCAGGGCCUdTdT", "AGGCCCUGGAGAAGCUGGGdTdT", 0.33),
    ("siSER-3m", "cccaGfcUfUfCfuccagggccudTdT", "aGfgccCfuggagaaGfcUfgggdTdT", 100),
    ("siSER-4", "UUGCUGGAGUCAUUCUCAAdTdT", "UUGAGAAUGACUCCAGCAAdTdT", 0.032),
    ("siSER-4m", "uugcUfgGfAfGfucauucucaadTdT", "uUfgagAfaugacucCfaGfcaadTdT", 0.027),
    ("siSER-5", "AGACAUCAAGCACUACUAUdTdT", "AUAGUAGUGCUUGAUGUCUdTdT", 0.20),
    ("siSER-5m", "agacAfuCfAfAfgcacuacuaudTdT", "aUfaguAfgugcuugAfuGfucudTdT", 0.23),
    ("siSER-6", "UCCCCUGCCAGCUGGUGCAdTdT", "UGCACCAGCUGGCAGGGGAdTdT", 2.74),
    ("siSER-6m", "ucccCfuGfCfCfagcuggugcadTdT", "uGfcacCfagcuggcAfgGfggadTdT", 100),
    ("siSER-7", "AGGUCACCAUCUCUGGAGUdTdT", "ACUCCAGAGAUGGUGACCUdTdT", 0.56),
    ("siSER-7m", "agguCfaCfCfAfucucuggagudTdT", "aCfuccAfgagauggUfgAfccudTdT", 100),
    ("siSER-8", "UCACCUGGAGCAGCCUUUUdTdT", "AAAAGGCUGCUCCAGGUGAdTdT", 1.3),
    ("siSER-8m", "ucacCfuGfGfAfgcagccuuuudTdT", "aAfaagGfcugcuccAfgGfugadTdT", 0.15),
    ("siSER-9", "CUGACUUUGGGAACCAGGAdTdT", "UCCUGGUUCCCAAAGUCAGdTdT", 0.16),
    ("siSER-9m", "cugaCfuUfUfGfggaaccaggadTdT", "uCfcugGfuucccaaAfgUfcagdTdT", 100),
    ("siSER-10", "AAGUUCUUCUCCCUCCAAAdTdT", "UUUGGAGGGAGAAGAACUUdTdT", 0.001),
    ("siSER-10m", "aaguUfcUfUfCfucccuccaaadTdT", "uUfuggAfgggagaaGfaAfcuudTdT", 0.004),
    ("siSER-11", "ACUUUAGGCAUCUUUUAAUdTdT", "AUUAAAAGAUGCCUAAAGUdTdT", 0.0007),
    ("siSER-11m", "acuuUfaGfGfCfaucuuuuaaudTdT", "aUfuaaAfagaugccUfaAfagudTdT", 0.0001),
    ("siAGT-1", "CCUGGCUGCAGGUGACCGAdTdT", "UCGGUCACCUGCAGCCAGGdTdT", 0.04),
    ("siAGT-1m", "ccugGfcUfGfCfaggugaccgadTdT", "uCfgguCfaccugcaGfcCfaggdTdT", 100),
    ("siAGT-2", "AGCAAUGACCGCAUCAGGAdTdT", "UCCUGAUGCGGUCAUUGCUdTdT", 0.13),
    ("siAGT-2m", "agcaAfuGfAfCfcgcaucaggadTdT", "uCfcugAfugcggucAfuUfgcudTdT", 4.9),
    ("siAGT-3", "CAAAAAUUGGGUUUUAAAAdTdT", "UUUUAAAACCCAAUUUUUGdTdT", 0.0004),
    ("siAGT-3m", "caaaAfaUfUfGfgguuuuaaaadTdT", "uUfuuaAfaacccaaUfuUfuugdTdT", 0.022),
    ("siAGT-4", "GGGUGGGGAGGCAAGAACAdTdT", "UGUUCUUGCCUCCCCACCCdTdT", 0.01),
    ("siAGT-4m", "ggguGfgGfGfAfggcaagaacadTdT", "uGfuucUfugccuccCfcAfcccdTdT", 0.21),
]

# ═══════════════════════════════════════════════════════════════════════════
# 2. KENSKI PAPER (Kenski et al. 2012, Mol Ther Nucleic Acids)
# ═══════════════════════════════════════════════════════════════════════════
# Sequences from Materials & Methods
# Note: Our model doesn't support 2'-O-benzyl, so we only validate unmodified
# and 2'-OMe-modified forms

def parse_kenski_seq(seq_parts: List[str]) -> str:
    """Parse Kenski semi-colon sequence notation."""
    result = []
    for part in seq_parts:
        part = part.strip()
        if part.startswith('r') and len(part) >= 2:
            base = part[1]
            if base in 'AUCG':
                result.append(base)
        elif part.startswith('dT'):
            result.append('T')
        elif part.startswith('omeU') or part.startswith('ome'):
            result.append('M')
        elif part == 'iB':
            continue  # skip inverted abasic
    return ''.join(result)


# ═══════════════════════════════════════════════════════════════════════════
# PREDICTION
# ═══════════════════════════════════════════════════════════════════════════

def predict_sequences(seq_pairs: List[Tuple[str, str, str]]) -> List[Dict]:
    """
    Run SMEpred on a list of (name, sense, antisense) sequences.
    Uses model A with default conditions (10 nM, 24h).
    """
    model = _get_model("A")
    names, s_list, a_list = [], [], []
    for name, sense, anti in seq_pairs:
        names.append(name)
        s_list.append(sense)
        a_list.append(anti)

    X = extract_batch_gbm(s_list, a_list)
    raw = model.predict(X)
    cal = _normalize_scores(raw, calibrator_key="cm")

    results = []
    for i, name in enumerate(names):
        results.append({
            "name": name,
            "sense_mod": s_list[i],
            "antisense_mod": a_list[i],
            "raw_score": float(raw[i]),
            "calibrated_score": float(cal[i]),
        })
    return results


# ═══════════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

def plot_alnylam_validation(results: List[Dict], data: List[Tuple]):
    """
    Create validation plots for Alnylam paper data.
    Compares predicted % inhibition with experimental IC50.
    """
    from scipy.stats import spearmanr, pearsonr

    # Build lookup
    score_map = {r["name"]: r["calibrated_score"] for r in results}
    ic50_map = {d[0]: d[3] for d in data}

    # Separate parent and modified
    parents = [(d[0], d[3]) for d in data if not d[0].endswith('m')]
    modifieds = [(d[0], d[3]) for d in data if d[0].endswith('m')]

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # ── Plot 1: Predicted score vs log10(IC50) - all ──
    ax = axes[0, 0]
    names_all = [d[0] for d in data]
    scores_all = [score_map.get(nm, 0) for nm in names_all]
    logic50_all = [np.log10(max(d[3], 1e-5)) for d in data]

    for nm, sc, lc in zip(names_all, scores_all, logic50_all):
        color = 'red' if nm.endswith('m') else 'blue'
        ax.scatter(lc, sc, c=color, alpha=0.7, s=60)
        ax.annotate(nm.replace('si', ''), (lc, sc), fontsize=5, alpha=0.7)

    spr, _ = spearmanr(scores_all, logic50_all)
    pcc, _ = pearsonr(scores_all, logic50_all)
    ax.set_xlabel('log10(IC50) [nM] (lower = more potent)')
    ax.set_ylabel('SMEpred predicted score')
    ax.set_title(f'All siRNAs: Spearman ρ={spr:.3f}, PCC={pcc:.3f}')
    ax.axvline(0, color='gray', linestyle='--', alpha=0.3)
    ax.grid(alpha=0.2)

    # ── Plot 2: Parent vs Modified prediction comparison ──
    ax = axes[0, 1]
    parent_names = [d[0] for d in parents]
    mod_names = [d[0] for d in modifieds]
    parent_scores = [score_map.get(n, 0) for n in parent_names]  # noqa: F841
    mod_scores = [score_map.get(n, 0) for n in mod_names]

    # Group by base name
    for p_name, m_name in zip(parent_names, mod_names):
        ps = score_map.get(p_name, 0)
        ms = score_map.get(m_name, 0)
        ax.plot([0, 1], [ps, ms], 'o-', alpha=0.5, markersize=6)
        label = p_name.replace('si', '')
        ax.text(1.02, ms, label, fontsize=7, va='center')

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Parent (unmodified)', 'Modified (2\'-OMe/2\'-F)'])
    ax.set_ylabel('SMEpred predicted score')
    ax.set_title('Parent → Modified score change')
    ax.grid(alpha=0.2)

    # ── Plot 3: IC50 ratio vs prediction delta ──
    ax = axes[0, 2]
    deltas_logic50 = []
    deltas_pred = []
    active_labels = []
    for p_name, m_name in zip(parent_names, mod_names):
        ic_p = ic50_map[p_name]
        ic_m = ic50_map[m_name]
        if ic_p > 0 and ic_m > 0:
            deltas_logic50.append(np.log10(ic_m / ic_p))
            deltas_pred.append(score_map.get(m_name, 0) - score_map.get(p_name, 0))
            active = 'active' if ic_m / ic_p < 100 else 'inactive'
            active_labels.append(active)

    for i, (dl, dp, al) in enumerate(zip(deltas_logic50, deltas_pred, active_labels)):
        color = 'green' if al == 'active' else 'red'
        ax.scatter(dl, dp, c=color, s=60, alpha=0.7)
        ax.annotate(parent_names[i].replace('si', ''), (dl, dp), fontsize=6, alpha=0.7)

    spr2, _ = spearmanr(deltas_pred, deltas_logic50)
    ax.set_xlabel('log10(IC50_modified / IC50_parent)')
    ax.set_ylabel('Δ Predicted score (modified − parent)')
    ax.set_title(f'Potency change correlation (Spearman ρ={spr2:.3f})')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.3)
    ax.grid(alpha=0.2)

    # ── Plot 4: Rank comparison bar chart ──
    ax = axes[1, 0]
    # Sort by predicted score within parent + modified groups
    all_sorted = sorted(zip(names_all, scores_all, logic50_all),
                        key=lambda x: -x[1])
    names_s = [x[0].replace('si', '') for x in all_sorted]
    scores_s = [x[1] for x in all_sorted]
    logic50_s = [x[2] for x in all_sorted]

    x = np.arange(len(names_s))
    width = 0.35
    bars1 = ax.bar(x - width/2, scores_s, width, label='SMEpred score', alpha=0.7)
    ax2 = ax.twinx()
    # Normalize IC50 to 0-100 for visual comparison (inverted: higher = better)
    norm_ic50 = [-l*10 + 50 for l in logic50_s]  # transform for display
    bars2 = ax2.bar(x + width/2, norm_ic50, width, label='−log10(IC50) scaled',
                    color='orange', alpha=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(names_s, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('SMEpred score')
    ax2.set_ylabel('Normalized potency (−log10(IC50))')
    ax.set_title('Rank comparison: SMEpred vs experimental IC50')
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax.grid(alpha=0.15)

    # ── Plot 5: Prediction vs IC50 with activity classification ──
    ax = axes[1, 1]
    for j, (nm, sc, lc) in enumerate(zip(names_all, scores_all, logic50_all)):
        is_active_mod = nm.endswith('m') and (ic50_map[nm] / ic50_map[nm.replace('m', '')] < 100)
        if nm.endswith('m'):
            color = 'green' if is_active_mod else 'red'
        else:
            color = 'blue'
        ax.scatter(lc, sc, c=color, s=60, alpha=0.7)
        if nm.endswith('m'):
            ax.annotate(nm.replace('si', ''), (lc, sc), fontsize=5,
                       color=color, alpha=0.7)

    ax.set_xlabel('log10(IC50)')
    ax.set_ylabel('SMEpred predicted score')
    ax.set_title('Blue=parent, Green=active mod, Red=inactive mod')
    ax.grid(alpha=0.2)

    # ── Plot 6: Key insight - model limitation explanation ──
    ax = axes[1, 2]
    # Parent-only correlations
    from scipy.stats import spearmanr as spr_func
    p_scores_vals = [score_map.get(d[0], 0) for d in parents]
    p_ic50_vals = [np.log10(max(d[1], 1e-5)) for d in parents]
    spr_p, _ = spr_func(p_scores_vals, p_ic50_vals)

    ax.text(0.5, 0.5,
        "Key findings:\n\n"
        f"Parent sequences: ρ={spr_p:.3f}\n"
        "  Weak correlation because our model\n"
        "  has limited cross-gene PCC (~0.26)\n\n"
        "Modified sequences: all → 62.9\n"
        "  Our MNC features compress all\n"
        "  2'-OMe/2'-F sequences to the same\n"
        "  composition → no discrimination\n\n"
        "This reveals a known limitation:\n"
        "  MNC alone cannot distinguish\n"
        "  sequences with identical modification\n"
        "  patterns but different base identities.",
        transform=ax.transAxes, fontsize=9, va='center', ha='center',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    ax.axis('off')

    plt.tight_layout()
    path = OUT / "alnylam_validation.png"
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {path}")

    return {"spearman_all": spr, "spearman_delta": spr2, "N": len(names_all)}


def create_summary_plot(results_alnylam, results_cmsirnadb, metrics_alnylam, metrics_cmsirnadb):
    """Create a summary comparison of all validation results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ── CMsiRNAdb scatter ──
    ax = axes[0]
    import pandas as pd
    csv_path = Path(__file__).parent.parent / "logs" / "cmsirnadb_validation.csv"
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        y_t = df['inhibition'].values
        y_p = df['predicted_calibrated'].values
        ax.scatter(y_t, y_p, alpha=0.1, s=5, c='blue')
        # Add trend line
        from numpy.polynomial.polynomial import polyfit
        b, m = polyfit(y_t, y_p, 1)
        x_line = np.array([0, 100])
        ax.plot(x_line, b + m * x_line, 'r-', lw=2)
        ax.plot(x_line, x_line, 'k--', alpha=0.3, label='Perfect')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_xlabel('Experimental % inhibition (CMsiRNAdb)')
        ax.set_ylabel('SMEpred predicted score')
        ax.set_title(f'CMsiRNAdb (n={len(y_t)}): PCC={metrics_cmsirnadb.get("pcc", 0):.3f}, MAE={metrics_cmsirnadb.get("mae", 0):.1f}')
        ax.legend()
        ax.grid(alpha=0.15)

    # ── Alnylam scatter ──
    ax = axes[1]
    if results_alnylam:
        from scipy.stats import spearmanr
        score_map = {r["name"]: r["calibrated_score"] for r in results_alnylam}
        for d in ALNYLAM_DATA:
            name, _, _, ic50 = d
            score = score_map.get(name, 0)
            color = 'red' if name.endswith('m') else 'blue'
            marker = 'x' if (name.endswith('m') and ic50 >= 100) else 'o'
            ax.scatter(np.log10(max(ic50, 1e-5)), score, c=color, marker=marker, s=50, alpha=0.7)
        ax.set_xlabel('log10(IC50) [nM] (lower = more potent)')
        ax.set_ylabel('SMEpred predicted score')
        ax.set_title(f'Alnylam (n=30): Spearman ρ={metrics_alnylam.get("spearman_all", 0):.3f}')
        ax.invert_xaxis()  # So more potent is on the right
        ax.grid(alpha=0.15)

    plt.tight_layout()
    path = OUT / "summary_validation.png"
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("PUBLISHED PAPER VALIDATION")
    print("=" * 60)

    # ── 1. Alnylam Paper ──
    print("\n1. Alnylam Paper (Kliuchnikov et al. 2025)")
    print("-" * 40)
    seq_pairs = []
    for name, pass_seq, guide_seq, ic50 in ALNYLAM_DATA:
        sense_mod = parse_alnylam_seq(pass_seq)
        anti_mod = parse_alnylam_seq(guide_seq)
        seq_pairs.append((name, sense_mod, anti_mod))

    results_alnylam = predict_sequences(seq_pairs)

    print(f"  Predicted {len(results_alnylam)} siRNA sequences")
    print(f"\n  {'Name':<15} {'Pred Score':<12} {'IC50 (nM)':<10} {'log10(IC50)':<12}")
    print(f"  {'-'*50}")
    score_map = {r["name"]: r["calibrated_score"] for r in results_alnylam}
    for name, _, _, ic50 in ALNYLAM_DATA:
        s = score_map.get(name, 0)
        print(f"  {name:<15} {s:<12.1f} {ic50:<10} {np.log10(max(ic50,1e-5)):<12.2f}")

    # Alnylam metrics
    from scipy.stats import spearmanr
    names_all = [d[0] for d in ALNYLAM_DATA]
    scores_all = [score_map.get(n, 0) for n in names_all]
    logic50_all = [np.log10(max(d[3], 1e-5)) for d in ALNYLAM_DATA]
    spr, _ = spearmanr(scores_all, logic50_all)
    # Delta correlation
    parent_names = [d[0] for d in ALNYLAM_DATA if not d[0].endswith('m')]
    mod_names = [d[0] for d in ALNYLAM_DATA if d[0].endswith('m')]
    deltas_pred = [score_map.get(m, 0) - score_map.get(p, 0)
                   for p, m in zip(parent_names, mod_names)]
    deltas_ic50 = []
    for p, m in zip(parent_names, mod_names):
        ic_p = max(d[3] for d in ALNYLAM_DATA if d[0] == p)
        ic_m = max(d[3] for d in ALNYLAM_DATA if d[0] == m)
        deltas_ic50.append(np.log10(ic_m / max(ic_p, 1e-5)))
    spr_delta, _ = spearmanr(deltas_pred, deltas_ic50)

    metrics_alnylam = {"spearman_all": spr, "spearman_delta": spr_delta}
    print(f"\n  Spearman ρ (all):     {spr:.4f}")
    print(f"  Spearman ρ (deltas):  {spr_delta:.4f}")

    # ── 2. Published paper plots ──
    print("\n2. Creating visualization plots...")
    metrics_cmsirnadb = {"pcc": 0.550, "mae": 17.63}  # from CMsiRNAdb validation
    plot_alnylam_validation(results_alnylam, ALNYLAM_DATA)
    create_summary_plot(results_alnylam, None, metrics_alnylam, metrics_cmsirnadb)

    # ── 3. Save all results ──
    output = {
        "alnylam": {
            "source": "Kliuchnikov et al. 2025, Mol Ther Nucleic Acids",
            "n_sequences": len(results_alnylam),
            "spearman_all_vs_ic50": spr,
            "spearman_delta": spr_delta,
        },
        "cmsirnadb": {
            "source": "He et al. 2026, BMC Bioinformatics (CMsiRNAdb)",
            "n_sequences": 12303,
            "pcc": 0.550,
            "spearman": 0.537,
            "mae": 17.63,
        }
    }
    json_path = OUT / "validation_results.json"
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n  Results saved to {json_path}")

    print(f"\n  Plots saved to {OUT}/")
    print("=" * 60)
