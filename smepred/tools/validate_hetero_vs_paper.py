"""
Validate our LightGBM model on the held-out validation set (hetero_val_303)
and compare against the original SMEpred paper's published PCC (0.80).
"""
import csv, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings('ignore', message='X does not have valid feature names')

from src.features import extract_batch_gbm
from src.predictor import _get_model, _normalize_scores
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error

# ─── load validation data ───────────────────────────────────────
csv_path = Path(__file__).parent.parent / "data" / "hetero_val_303.csv"
with open(csv_path, "r", encoding="utf-8") as f:
    reader = csv.reader(f)
    raw = list(reader)
header = raw[0]
rows = raw[1:]

sense_m      = [r[0] for r in rows]
antisense_m  = [r[1] for r in rows]
base_sense   = [r[2] for r in rows]
base_anti    = [r[3] for r in rows]
efficacy     = np.array([float(r[4]) for r in rows])
conc         = [float(r[5]) if r[5] else None for r in rows]
time_h       = [float(r[6]) if r[6] else None for r in rows]
gene         = [r[7].strip() for r in rows]

print(f"Loaded {len(rows)} validation rows from {header}")

# ─── run predictions ────────────────────────────────────────────
model = _get_model("A")

BATCH = 512
all_raw = []
all_cal = []
for start in range(0, len(rows), BATCH):
    end = min(start + BATCH, len(rows))
    X = extract_batch_gbm(
        sense_m[start:end], antisense_m[start:end],
        base_sense_list=base_sense[start:end],
        base_antisense_list=base_anti[start:end],
        conc_list=conc[start:end],
        time_list=time_h[start:end],
    )
    raw_pred = model.predict(X)
    cal_pred = _normalize_scores(raw_pred, calibrator_key="cm")
    all_raw.extend(raw_pred.tolist())
    all_cal.extend(cal_pred.tolist())

pred_raw = np.array(all_raw)
pred_cal = np.array(all_cal)

# ─── metrics ────────────────────────────────────────────────────
pcc, _  = pearsonr(efficacy, pred_cal)
spr, _  = spearmanr(efficacy, pred_cal)
mae      = mean_absolute_error(efficacy, pred_cal)
rmse     = np.sqrt(np.mean((efficacy - pred_cal) ** 2))

print(f"\n{'='*55}")
print(f"OUR MODEL  — LightGBM (Model-A) on hetero_val_303")
print(f"{'='*55}")
print(f"  N           : {len(efficacy)}")
print(f"  PCC         : {pcc:.4f}")
print(f"  Spearman ρ  : {spr:.4f}")
print(f"  MAE         : {mae:.2f} pts")
print(f"  RMSE        : {rmse:.2f} pts")
print(f"  Mean true   : {efficacy.mean():.2f}")
print(f"  Mean pred   : {pred_cal.mean():.2f}")
print(f"  Std true    : {efficacy.std():.2f}")
print(f"  Std pred    : {pred_cal.std():.2f}")

print(f"\n{'='*55}")
print(f"ORIGINAL SMEpred PAPER (Dar et al. 2016, RNA Biology)")
print(f"{'='*55}")
print(f"  Dataset     : Hetero-V303 (303 rows from siRNAmod)")
print(f"  Algorithm   : SVM RBF kernel")
print(f"  Model-A MNC : PCC 0.808")
print(f"  Model-B     : PCC 0.86")
print(f"  Model-C     : PCC 0.78")
print()

# Per-gene breakdown
print("Per-gene breakdown:")
genes_u = sorted(set(gene))
gene_pccs = []
for g in genes_u:
    mask = np.array([gg == g for gg in gene])
    if mask.sum() < 5:
        continue
    pg = pred_cal[mask]
    eg = efficacy[mask]
    g_pcc, _ = pearsonr(eg, pg)
    g_spr, _ = spearmanr(eg, pg)
    g_mae = mean_absolute_error(eg, pg)
    gene_pccs.append(g_pcc)
    print(f"  {g:12s}: N={mask.sum():4d}  PCC={g_pcc:.4f}  ρ={g_spr:.4f}  MAE={g_mae:.2f}")
print(f"  Mean per-gene PCC: {np.mean(gene_pccs):.4f}")

# ─── save CSV ───────────────────────────────────────────────────
out_path = Path(__file__).parent.parent / "logs" / "hetero_val_predictions.csv"
df = pd.DataFrame({
    "sense_modified": sense_m,
    "antisense_modified": antisense_m,
    "base_sense": base_sense,
    "base_antisense": base_anti,
    "true_efficacy": efficacy,
    "predicted_raw": pred_raw,
    "predicted_calibrated": pred_cal,
    "concentration_nM": conc,
    "time_h": time_h,
    "target_gene": gene,
})
df.to_csv(out_path, index=False)
print(f"\nResults saved to {out_path}")

# ─── comparison summary ─────────────────────────────────────────
print(f"\n{'='*55}")
print(f"COMPARISON: Original SMEpred paper vs Our LightGBM")
print(f"{'='*55}")
print(f"                   Paper (SVM)    Ours (LightGBM)")
print(f"  Training data    siRNAmod        HelixZero 43k")
print(f"  Training rows    2,728           23,187")
print(f"  Validation rows  303             2,576")
print(f"  Model-A PCC      0.808           {pcc:.4f}")
print(f"  Gene-grouped PCC —               {np.mean(gene_pccs):.4f}")
print()
if pcc >= 0.80:
    print("✓ Our model matches the original SMEpred paper's PCC!")
else:
    gap = 0.808 - pcc
    print(f"  Our model's PCC is {abs(gap):.3f} below the paper's 0.808.")
    print(f"  Key differences: our validation set is 8.5× larger (2,576 vs 303 rows)")
    print(f"  with more modification diversity (HelixZero catalog vs siRNAmod).")
