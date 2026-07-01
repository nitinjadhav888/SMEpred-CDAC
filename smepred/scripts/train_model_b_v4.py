"""
Train Model B v4 — Expanded feature space + all data.
Now uses extract_positional_features_batch for ALL rows (not CSV pre-computed).
_MODIFICATION_MAP expanded: F, M, L, S, D + E(MOE) + 2,3,4,6,8,Q,U,X (HelixZero codes).
"""

import pandas as pd, numpy as np, lightgbm as lgb, json, joblib, warnings, sys
from pathlib import Path
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr, spearmanr
warnings.filterwarnings('ignore')


from src.features import extract_positional_features_batch, _MODIFICATION_MAP, _MOD_CATEGORIES

MODELS_DIR = Path(__file__).parent
SRCDIR = Path(__file__).parent.parent

print("Model B v4 — Expanded feature extractor")
print("  _MODIFICATION_MAP:", dict(_MODIFICATION_MAP))
print("  _MOD_CATEGORIES:", _MOD_CATEGORIES)
n_flags = len(_MODIFICATION_MAP) + 2
n_pos = n_flags * 21 * 2
n_global = (len(_MOD_CATEGORIES) + 9) * 2
n_total = n_pos + n_global + 1
print(f"  Features per position: {n_flags}  Total dim: {n_total}")
print()

# ── 1. Load & featurize position-aware ──
print("[1] Loading position-aware dataset...")
# Note: this CSV was moved to data/raw/ during cleanup.
# Symlink or copy it back from backup if re-training.
pdf = pd.read_csv(SRCDIR / "data" / "raw" / "sirna_modified_position_aware_dataset_v2.csv", low_memory=False)
print(f"  {len(pdf):,} rows")

# Use extract_positional_features_batch for ALL data (consistent feature space)
sense_list = pdf['sense_seq'].tolist()
anti_list = pdf['antisense_seq'].tolist()
bs_list = pdf['base_sense'].tolist()
ba_list = pdf['base_antisense'].tolist()
conc_list = pdf['concentration_nM'].fillna(10).tolist()
y_pos = pdf['efficacy'].values.astype(np.float32)

print(f"  Extracting features for {len(sense_list):,} rows...")
X_pos = extract_positional_features_batch(sense_list, anti_list, bs_list, ba_list, conc_list)
print(f"  Done. Shape: {X_pos.shape}")

# ── 2. Load & featurize hetero_train ──
print("\n[2] Loading hetero_train_2728...")
hdf = pd.read_csv(SRCDIR / "data" / "processed" / "hetero_train_2728.csv", low_memory=False)
hdf.columns = [c.strip() for c in hdf.columns]
print(f"  {len(hdf):,} rows")
X_hetero = extract_positional_features_batch(
    hdf['sense'].tolist(), hdf['antisense'].tolist(),
    hdf['base_sense'].tolist(), hdf['base_antisense'].tolist(),
    hdf['concentration_nM'].fillna(10).tolist()
)
y_hetero = hdf['efficacy'].values.astype(np.float32)
print(f"  Shape: {X_hetero.shape}")

# ── 3. Load & featurize CMsiRNAdb ──
print("\n[3] Loading CMsiRNAdb...")
cdf = pd.read_csv(SRCDIR / "data" / "processed" / "cmsirnadb_full.csv", low_memory=False)
cdf = cdf.dropna(subset=['efficacy', 'sense', 'antisense', 'base_sense', 'base_antisense'])
cdf['_key'] = cdf['sense'].str[:20] + cdf['antisense'].str[:20]
pdf['_key'] = pdf['sense_seq'].str[:20] + pdf['antisense_seq'].str[:20]
cdf_new = cdf[~cdf['_key'].isin(pdf['_key'])].copy()
if len(cdf_new) > 5000:
    cdf_new = cdf_new.sample(n=5000, random_state=42)
print(f"  {len(cdf_new):,} rows (non-overlapping, sampled)")
X_cms = extract_positional_features_batch(
    cdf_new['sense'].tolist(), cdf_new['antisense'].tolist(),
    cdf_new['base_sense'].tolist(), cdf_new['base_antisense'].tolist(),
    cdf_new['concentration_nM'].fillna(10).tolist()
)
y_cms = cdf_new['efficacy'].values.astype(np.float32)
print(f"  Shape: {X_cms.shape}")

# ── 4. Combine ──
print("\n[4] Combining datasets...")
X_all = np.vstack([X_pos, X_hetero, X_cms])
y_all = np.concatenate([y_pos, y_hetero, y_cms])
src_labels = (['pos'] * len(X_pos) + ['het'] * len(X_hetero) + ['cms'] * len(X_cms))

sss = StratifiedShuffleSplit(n_splits=1, test_size=0.15, random_state=42)
trainval_idx, test_idx = next(sss.split(X_all, src_labels))
X_trainval, X_test = X_all[trainval_idx], X_all[test_idx]
y_trainval, y_test = y_all[trainval_idx], y_all[test_idx]
src_test = [src_labels[i] for i in test_idx]

sss2 = StratifiedShuffleSplit(n_splits=1, test_size=0.176, random_state=42)
src_trainval = [src_labels[i] for i in trainval_idx]
train_idx, val_idx = next(sss2.split(X_trainval, src_trainval))
X_train, X_val = X_trainval[train_idx], X_trainval[val_idx]
y_train, y_val = y_trainval[train_idx], y_trainval[val_idx]

print(f"\n  Train: {len(X_train):,}  Val: {len(X_val):,}  Test: {len(X_test):,}")

# ── 5. Load hetero_val_303 ──
print("\n[5] Loading hetero_val_303...")
vdf = pd.read_csv(SRCDIR / "data" / "processed" / "hetero_val_303.csv", low_memory=False)
X_val_ext = extract_positional_features_batch(
    vdf['sense'].tolist(), vdf['antisense'].tolist(),
    vdf['base_sense'].tolist(), vdf['base_antisense'].tolist(),
    vdf['concentration_nM'].fillna(10).tolist()
)
y_val_ext = vdf['efficacy'].values.astype(np.float32)
v_genes = vdf['target_gene'].tolist() if 'target_gene' in vdf.columns else None
print(f"  {len(vdf)} rows, features: {X_val_ext.shape}")

# ── 6. Train ──
print("\n[6] Training LightGBM...")
feature_names = [f"f{i}" for i in range(X_train.shape[1])]
lgb_train = lgb.Dataset(X_train, label=y_train, feature_name=feature_names)
lgb_val_ds = lgb.Dataset(X_val, label=y_val, reference=lgb_train)

params = {
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

model = lgb.train(
    params, lgb_train,
    num_boost_round=params['n_estimators'],
    valid_sets=[lgb_val_ds],
    callbacks=[lgb.early_stopping(100, verbose=True), lgb.log_evaluation(200)],
)

# ── 7. Evaluate ──
def print_metrics(y_true, y_pred, label, genes=None):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_t, y_p = y_true[mask], y_pred[mask]
    if len(y_t) < 5: return None
    pcc, _ = pearsonr(y_t, y_p)
    spr, _ = spearmanr(y_t, y_p)
    mae = mean_absolute_error(y_t, y_p)
    rmse = float(np.sqrt(np.mean((y_t - y_p) ** 2)))
    r2 = r2_score(y_t, y_p)
    print(f"  {label}: N={len(y_t):,}  PCC={pcc:.4f}  ρ={spr:.4f}  MAE={mae:.2f}  RMSE={rmse:.2f}  R²={r2:.4f}")
    if genes is not None:
        g_pccs = []
        for g in sorted(set(str(x) for x in genes if x and str(x).strip() not in ('nan','','None'))):
            gm = np.array([str(gg).strip() == g for gg in genes])
            if gm.sum() < 5: continue
            gp,_ = pearsonr(y_t[gm], y_p[gm])
            g_pccs.append(gp)
            gmae = mean_absolute_error(y_t[gm], y_p[gm])
            print(f"    {g:12s}: N={gm.sum():4d}  PCC={gp:.4f}  MAE={gmae:.2f}")
        if g_pccs:
            print(f"    Mean per-gene PCC: {np.mean(g_pccs):.4f}")
    return {"pcc": pcc, "spearman": spr, "mae": mae, "rmse": rmse, "r2": r2}

print(f"\n  --- Internal evaluation ---")
y_pred_test = model.predict(X_test)
test_metrics = print_metrics(y_test, y_pred_test, "Test (combined holdout)")
for src in sorted(set(src_test)):
    sm = np.array([s == src for s in src_test])
    if sm.sum() < 5: continue
    print_metrics(y_test[sm], y_pred_test[sm], f"  {src}")

print(f"\n  --- External validation: hetero_val_303 ---")
y_pred_v = model.predict(X_val_ext)
val_metrics = print_metrics(y_val_ext, y_pred_v, "Model B v4", v_genes)

print(f"\n  Historical:")
print(f"  Model A:      PCC=0.773  MAE=14.31")
print(f"  Model B v1:   PCC=0.573  MAE=18.56")
print(f"  Model B v2:   PCC=0.636  MAE=17.31")
if val_metrics:
    print(f"  Model B v4:   PCC={val_metrics['pcc']:.4f}  MAE={val_metrics['mae']:.2f}")

# ── 8. Save ──
model.save_model(str(MODELS_DIR / 'model_b.txt'))
joblib.dump(model, MODELS_DIR / 'model_b.pkl')
meta = {
    'version': 4,
    'date': pd.Timestamp.now().isoformat(),
    'mod_char_map': dict(_MODIFICATION_MAP),
    'mod_types': _MOD_CATEGORIES,
    'n_features': X_train.shape[1],
    'n_flags_per_pos': n_flags,
    'training_rows': {'position_aware': len(X_pos), 'hetero_patent': len(X_hetero), 'cmsirnadb': len(X_cms), 'total': len(X_all)},
    'best_iteration': model.best_iteration,
    'test_metrics': {k: float(v) for k,v in (test_metrics or {}).items()},
    'hetero_val_303': {k: float(v) for k,v in (val_metrics or {}).items()},
    'params': {k: str(v) if isinstance(v, (list, dict)) else v for k, v in params.items()},
}
with open(MODELS_DIR / 'model_b_meta.json', 'w') as f:
    json.dump(meta, f, indent=2)
print(f"\nSaved: {MODELS_DIR / 'model_b.pkl'}")
