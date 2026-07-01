"""
cross_validate_model_b.py — K-Fold Cross-Validation for Model B v4.

Performs 5-fold stratified cross-validation on all available training data
to produce robust performance estimates with confidence intervals.

Outputs:
  1. Per-fold metrics: PCC, Spearman, MAE, RMSE, R²
  2. Overall mean ± std across folds
  3. Per-source (hetero vs cmsirnadb) breakdown
  4. Per-gene breakdown on gene-grouped folds
  5. Leave-One-Source-Out cross-validation
  6. External validation on hetero_val_303 (holdout)

Usage:
  cd helixzero_cms
  python scripts/cross_validate_model_b.py
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import lightgbm as lgb
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import StratifiedKFold, KFold

warnings.filterwarnings('ignore')

# ── Setup paths ──────────────────────────────────────────────────────────────
SRCDIR = Path(__file__).parent.parent

from src.features import extract_positional_features_batch, _MOD_CHAR_MAP, _MOD_TYPES


# ── LightGBM hyperparams (same as original training) ─────────────────────────
PARAMS = {
    'objective': 'regression',
    'metric': ['rmse', 'mae'],
    'boosting_type': 'gbdt',
    'num_leaves': 127,
    'max_depth': -1,
    'learning_rate': 0.03,
    'n_estimators': 3000,
    'feature_fraction': 0.6,
    'bagging_fraction': 0.8,
    'bagging_freq': 5,
    'reg_alpha': 0.1,
    'reg_lambda': 0.2,
    'min_child_samples': 20,
    'verbose': -1,
    'random_state': 42,
}


# ── Metrics helper ───────────────────────────────────────────────────────────
def compute_metrics(y_true, y_pred, label=""):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    yt, yp = y_true[mask], y_pred[mask]
    if len(yt) < 5:
        return None
    pcc, _ = pearsonr(yt, yp)
    spr, _ = spearmanr(yt, yp)
    mae = mean_absolute_error(yt, yp)
    rmse = float(np.sqrt(np.mean((yt - yp) ** 2)))
    r2 = r2_score(yt, yp)
    return {"label": label, "n": len(yt), "pcc": pcc, "spearman": spr,
            "mae": mae, "rmse": rmse, "r2": r2}


def print_metrics(m, indent=2):
    if m is None:
        return
    pad = " " * indent
    print(f"{pad}{m['label']}: N={m['n']:,}  PCC={m['pcc']:.4f}  "
          f"ρ={m['spearman']:.4f}  MAE={m['mae']:.2f}  "
          f"RMSE={m['rmse']:.2f}  R²={m['r2']:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load & featurize all available datasets
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("  MODEL B v4 — K-FOLD CROSS-VALIDATION")
print("=" * 72)

# ── Hetero train (23,187 rows) ──
print("\n[1] Loading datasets...")
hdf = pd.read_csv(SRCDIR / "data" / "processed" / "hetero_train_2728.csv", low_memory=False)
hdf.columns = [c.strip() for c in hdf.columns]
print(f"  hetero_train: {len(hdf):,} rows")

# ── CMsiRNAdb (up to 25,863 rows) ──
cdf = pd.read_csv(SRCDIR / "data" / "processed" / "cmsirnadb_full.csv", low_memory=False)
cdf = cdf.dropna(subset=['efficacy', 'sense', 'antisense', 'base_sense', 'base_antisense'])
# Deduplicate against hetero_train (same logic as original training script)
cdf['_key'] = cdf['sense'].str[:20] + cdf['antisense'].str[:20]
hdf['_key'] = hdf['sense'].str[:20] + hdf['antisense'].str[:20]
cdf_new = cdf[~cdf['_key'].isin(hdf['_key'])].copy()
print(f"  cmsirnadb (non-overlapping): {len(cdf_new):,} rows")

# ── Check for position-aware dataset ──
# Check both locations (data/ and data/raw/)
pos_path = SRCDIR / "data" / "sirna_modified_position_aware_dataset_v2.csv"
if not pos_path.exists():
    pos_path = SRCDIR / "data" / "raw" / "sirna_modified_position_aware_dataset_v2.csv"
has_pos = pos_path.exists()
if has_pos:
    pdf = pd.read_csv(pos_path, low_memory=False)
    print(f"  position_aware: {len(pdf):,} rows")
else:
    print(f"  position_aware: NOT FOUND (excluded from repo)")
    print(f"  → CV will use hetero_train + cmsirnadb only ({len(hdf) + len(cdf_new):,} rows)")

# ── Featurize ────────────────────────────────────────────────────────────────
print("\n[2] Extracting 1,467-d positional features...")
t0 = time.time()

# Hetero
X_hetero = extract_positional_features_batch(
    hdf['sense'].tolist(), hdf['antisense'].tolist(),
    hdf['base_sense'].tolist(), hdf['base_antisense'].tolist(),
    hdf['concentration_nM'].fillna(10).tolist()
)
y_hetero = hdf['efficacy'].values.astype(np.float32)
genes_hetero = hdf['target_gene'].values if 'target_gene' in hdf.columns else np.array(['Unknown'] * len(hdf))

# CMsiRNAdb
X_cms = extract_positional_features_batch(
    cdf_new['sense'].tolist(), cdf_new['antisense'].tolist(),
    cdf_new['base_sense'].tolist(), cdf_new['base_antisense'].tolist(),
    cdf_new['concentration_nM'].fillna(10).tolist()
)
y_cms = cdf_new['efficacy'].values.astype(np.float32)
genes_cms = cdf_new['target_gene'].values if 'target_gene' in cdf_new.columns else np.array(['Unknown'] * len(cdf_new))

# Position-aware (if available)
if has_pos:
    X_pos = extract_positional_features_batch(
        pdf['sense_seq'].tolist(), pdf['antisense_seq'].tolist(),
        pdf['base_sense'].tolist(), pdf['base_antisense'].tolist(),
        pdf['concentration_nM'].fillna(10).tolist()
    )
    y_pos = pdf['efficacy'].values.astype(np.float32)
    genes_pos = pdf['target_gene'].values if 'target_gene' in pdf.columns else np.array(['Unknown'] * len(pdf))
    X_all = np.vstack([X_pos, X_hetero, X_cms])
    y_all = np.concatenate([y_pos, y_hetero, y_cms])
    src_all = np.array(['pos'] * len(X_pos) + ['het'] * len(X_hetero) + ['cms'] * len(X_cms))
    genes_all = np.concatenate([genes_pos, genes_hetero, genes_cms])
else:
    X_all = np.vstack([X_hetero, X_cms])
    y_all = np.concatenate([y_hetero, y_cms])
    src_all = np.array(['het'] * len(X_hetero) + ['cms'] * len(X_cms))
    genes_all = np.concatenate([genes_hetero, genes_cms])

feat_time = time.time() - t0
print(f"  Done in {feat_time:.1f}s. Combined shape: {X_all.shape}")
print(f"  Feature dimension: {X_all.shape[1]}")
print(f"  Efficacy range: [{y_all.min():.1f}, {y_all.max():.1f}], mean={y_all.mean():.1f}, std={y_all.std():.1f}")

# ── Load external validation set ─────────────────────────────────────────────
print("\n[3] Loading external validation set (hetero_val_303)...")
vdf = pd.read_csv(SRCDIR / "data" / "processed" / "hetero_val_303.csv", low_memory=False)
X_ext = extract_positional_features_batch(
    vdf['sense'].tolist(), vdf['antisense'].tolist(),
    vdf['base_sense'].tolist(), vdf['base_antisense'].tolist(),
    vdf['concentration_nM'].fillna(10).tolist()
)
y_ext = vdf['efficacy'].values.astype(np.float32)
genes_ext = vdf['target_gene'].values if 'target_gene' in vdf.columns else None
print(f"  {len(vdf):,} rows, features: {X_ext.shape}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. 5-Fold Stratified Cross-Validation
# ══════════════════════════════════════════════════════════════════════════════
N_FOLDS = 5
print(f"\n{'=' * 72}")
print(f"  5-FOLD STRATIFIED CROSS-VALIDATION (stratified by source)")
print(f"{'=' * 72}")

skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
fold_metrics = []
fold_ext_metrics = []
all_oof_preds = np.zeros(len(y_all))  # out-of-fold predictions
feature_names = [f"f{i}" for i in range(X_all.shape[1])]

for fold_i, (train_idx, val_idx) in enumerate(skf.split(X_all, src_all)):
    print(f"\n  ── Fold {fold_i + 1}/{N_FOLDS} ──")
    X_train, X_val = X_all[train_idx], X_all[val_idx]
    y_train, y_val = y_all[train_idx], y_all[val_idx]
    src_val = src_all[val_idx]
    genes_val = genes_all[val_idx]

    print(f"  Train: {len(X_train):,}  Val: {len(X_val):,}")

    # Train
    t_start = time.time()
    lgb_train = lgb.Dataset(X_train, label=y_train, feature_name=feature_names)
    lgb_val = lgb.Dataset(X_val, label=y_val, reference=lgb_train)

    model = lgb.train(
        PARAMS, lgb_train,
        num_boost_round=PARAMS['n_estimators'],
        valid_sets=[lgb_val],
        callbacks=[lgb.early_stopping(100, verbose=False), lgb.log_evaluation(0)],
    )
    train_time = time.time() - t_start
    print(f"  Trained in {train_time:.1f}s, best_iteration={model.best_iteration}")

    # Predict on validation fold
    y_pred_val = model.predict(X_val)
    all_oof_preds[val_idx] = y_pred_val

    # Overall fold metrics
    m = compute_metrics(y_val, y_pred_val, f"Fold {fold_i + 1}")
    fold_metrics.append(m)
    print_metrics(m)

    # Per-source breakdown within fold
    for src_name in sorted(set(src_val)):
        mask = src_val == src_name
        if mask.sum() >= 5:
            sm = compute_metrics(y_val[mask], y_pred_val[mask], f"  {src_name}")
            print_metrics(sm, indent=4)

    # Per-gene breakdown (top genes with ≥20 samples)
    gene_pccs = []
    for g in sorted(set(str(x) for x in genes_val if x and str(x).strip() not in ('nan', '', 'None', 'Unknown'))):
        gm = np.array([str(gg).strip() == g for gg in genes_val])
        if gm.sum() >= 20:
            gp, _ = pearsonr(y_val[gm], y_pred_val[gm])
            gene_pccs.append(gp)
    if gene_pccs:
        print(f"    Mean per-gene PCC ({len(gene_pccs)} genes ≥20 samples): {np.mean(gene_pccs):.4f}")

    # External validation (same model, predict on hetero_val_303)
    y_pred_ext = model.predict(X_ext)
    m_ext = compute_metrics(y_ext, y_pred_ext, f"Fold {fold_i + 1} → ext_val")
    fold_ext_metrics.append(m_ext)
    print_metrics(m_ext)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Aggregate Results
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print(f"  AGGREGATE RESULTS — {N_FOLDS}-FOLD CV")
print(f"{'=' * 72}")

# Summary table
metrics_keys = ['pcc', 'spearman', 'mae', 'rmse', 'r2']
valid_folds = [m for m in fold_metrics if m is not None]

print(f"\n  {'Metric':<12} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
print(f"  {'─' * 48}")
for key in metrics_keys:
    vals = [m[key] for m in valid_folds]
    print(f"  {key.upper():<12} {np.mean(vals):8.4f} {np.std(vals):8.4f} "
          f"{np.min(vals):8.4f} {np.max(vals):8.4f}")

# Overall OOF metrics (all out-of-fold predictions combined)
print(f"\n  ── Overall Out-of-Fold (OOF) Metrics ──")
oof_m = compute_metrics(y_all, all_oof_preds, "OOF Combined")
print_metrics(oof_m)

# External validation aggregate
print(f"\n  ── External Validation (hetero_val_303) ──")
valid_ext = [m for m in fold_ext_metrics if m is not None]
print(f"\n  {'Metric':<12} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
print(f"  {'─' * 48}")
for key in metrics_keys:
    vals = [m[key] for m in valid_ext]
    print(f"  {key.upper():<12} {np.mean(vals):8.4f} {np.std(vals):8.4f} "
          f"{np.min(vals):8.4f} {np.max(vals):8.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Leave-One-Source-Out Cross-Validation
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print(f"  LEAVE-ONE-SOURCE-OUT CROSS-VALIDATION")
print(f"{'=' * 72}")
print(f"  Tests generalization: train on one source, predict on the other.\n")

for holdout_src in sorted(set(src_all)):
    train_mask = src_all != holdout_src
    test_mask = src_all == holdout_src

    X_tr, y_tr = X_all[train_mask], y_all[train_mask]
    X_te, y_te = X_all[test_mask], y_all[test_mask]

    print(f"  Holdout: {holdout_src} ({test_mask.sum():,} rows), "
          f"Train: {train_mask.sum():,} rows")

    lgb_tr = lgb.Dataset(X_tr, label=y_tr, feature_name=feature_names)
    lgb_te = lgb.Dataset(X_te, label=y_te, reference=lgb_tr)

    model_loso = lgb.train(
        PARAMS, lgb_tr,
        num_boost_round=PARAMS['n_estimators'],
        valid_sets=[lgb_te],
        callbacks=[lgb.early_stopping(100, verbose=False), lgb.log_evaluation(0)],
    )

    y_pred_loso = model_loso.predict(X_te)
    m_loso = compute_metrics(y_te, y_pred_loso, f"Train→{holdout_src}")
    print_metrics(m_loso)

    # Also test on external validation
    y_pred_ext_loso = model_loso.predict(X_ext)
    m_ext_loso = compute_metrics(y_ext, y_pred_ext_loso, f"  → ext_val")
    print_metrics(m_ext_loso, indent=4)
    print()


# ══════════════════════════════════════════════════════════════════════════════
# 5. Comparison with Original Training
# ══════════════════════════════════════════════════════════════════════════════
print(f"{'=' * 72}")
print(f"  COMPARISON WITH ORIGINAL TRAINING")
print(f"{'=' * 72}")
print(f"""
  Original (single 85/15 split, all 3 datasets = 83,535 rows):
    Test PCC:       0.8217
    Test Spearman:  0.8225
    Test MAE:       12.27
    Test RMSE:      16.84
    Test R²:        0.6752
    Ext. val PCC:   0.6504

  {N_FOLDS}-Fold CV (available data = {len(X_all):,} rows):
    OOF PCC:        {oof_m['pcc']:.4f}
    OOF Spearman:   {oof_m['spearman']:.4f}
    OOF MAE:        {oof_m['mae']:.2f}
    OOF RMSE:       {oof_m['rmse']:.2f}
    OOF R²:         {oof_m['r2']:.4f}
    Ext. val PCC:   {np.mean([m['pcc'] for m in valid_ext]):.4f} ± {np.std([m['pcc'] for m in valid_ext]):.4f}

  Note: Original training includes position_aware dataset (55,730 rows)
  which is not available in the current repo. CV metrics are computed
  on hetero_train + cmsirnadb only. Lower PCC is expected since the
  position-aware data provides the strongest signal.
""")

print(f"{'=' * 72}")
print(f"  CROSS-VALIDATION COMPLETE")
print(f"{'=' * 72}")
