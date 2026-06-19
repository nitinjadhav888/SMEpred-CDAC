"""
Train Model B on the position-aware modification dataset (55,730 rows).
Saves model_b.pkl and model_b_meta.json.
"""
import pandas as pd, numpy as np, lightgbm as lgb, json, joblib, warnings
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score, ndcg_score
warnings.filterwarnings('ignore')

MODELS_DIR = Path(__file__).parent
DATA_PATH = r"D:\Helixx\sirna_modified_position_aware_dataset_v2.csv"

df = pd.read_csv(DATA_PATH)
print(f"Loaded {len(df):,} rows × {df.shape[1]} cols")

# ── Feature columns (same selection as original training script) ──────
pos_binary_cols = [c for c in df.columns if
                   any(c.endswith(f) for f in ['_is_2F','_is_2OMe','_is_LNA',
                                                '_is_PS','_is_DNA',
                                                '_is_canonical','_is_modified'])]
global_cols = [c for c in df.columns if c.startswith(('ss_n_','as_n_',
               'ss_frac','as_frac','ss_seed','as_seed',
               'ss_gc','as_gc','ss_cleavage','as_cleavage',
               'ss_5term','as_5term','ss_3term','as_3term',
               'ss_alt','as_alt'))]
FEATURE_COLS = pos_binary_cols + global_cols + ['log_conc']
print(f"Features: {len(FEATURE_COLS)} ({len(pos_binary_cols)} pos + {len(global_cols)} global + 1 conc)")

# Prepare data
df['log_conc'] = np.log1p(df['concentration_nM'].fillna(10))
X = df[FEATURE_COLS].fillna(0).astype(np.float32)
y = df['efficacy'].values.astype(np.float32)

# ── Train/val/test split (80/10/10) ──────────────────────────────────
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.111, random_state=42)
print(f"Train: {len(X_train):,}  Val: {len(X_val):,}  Test: {len(X_test):,}")

# ── Data source stratified eval ──────────────────────────────────────
train_idx, test_idx = train_test_split(range(len(df)), test_size=0.1, random_state=42)
_, val_idx = train_test_split(train_idx, test_size=0.111, random_state=42)
test_sources = df.iloc[test_idx]['data_source'].value_counts()

# ── LightGBM training ────────────────────────────────────────────────
lgb_train = lgb.Dataset(X_train, label=y_train, feature_name=FEATURE_COLS)
lgb_val = lgb.Dataset(X_val, label=y_val, reference=lgb_train)

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

print(f"\nTraining LightGBM...")
model = lgb.train(
    params,
    lgb_train,
    num_boost_round=params['n_estimators'],
    valid_sets=[lgb_val],
    callbacks=[lgb.early_stopping(100, verbose=True), lgb.log_evaluation(200)],
)
best_iter = model.best_iteration

# ── Evaluate ─────────────────────────────────────────────────────────
y_pred = model.predict(X_test)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)
rmse = float(np.sqrt(np.mean((y_test - y_pred) ** 2)))
ndcg = ndcg_score([y_test], [y_pred], k=10)

print(f"\n{'='*50}")
print(f"Test RMSE: {rmse:.2f}")
print(f"Test MAE : {mae:.2f}")
print(f"Test R²  : {r2:.4f}")
print(f"NDCG@10  : {ndcg:.4f}")
print(f"{'='*50}")

# Per-source evaluation
print(f"\nPer-source test MAE:")
test_df = df.iloc[test_idx].copy()
test_df['pred'] = model.predict(X_test)
for src in test_df['data_source'].unique():
    sub = test_df[test_df['data_source'] == src]
    src_mae = mean_absolute_error(sub['efficacy'], sub['pred'])
    print(f"  {src:<45s}  MAE={src_mae:.2f}  n={len(sub):,}")

# ── Feature importance ───────────────────────────────────────────────
fi = pd.Series(model.feature_importance(importance_type='gain'), index=FEATURE_COLS)
print(f"\nTop 20 features by gain:")
print(fi.sort_values(ascending=False).head(20).to_string())

# ── Save model & metadata ────────────────────────────────────────────
model.save_model(str(MODELS_DIR / 'model_b.txt'))
meta = {
    'date': pd.Timestamp.now().isoformat(),
    'version': 1,
    'rows': len(df),
    'features': len(FEATURE_COLS),
    'feature_cols': FEATURE_COLS,
    'best_iteration': best_iter,
    'test_rmse': rmse,
    'test_mae': mae,
    'test_r2': r2,
    'test_ndcg10': ndcg,
    'params': {k: str(v) if isinstance(v, (list, dict)) else v for k, v in params.items()},
}
with open(MODELS_DIR / 'model_b_meta.json', 'w') as f:
    json.dump(meta, f, indent=2)

print(f"\nModel saved to {MODELS_DIR / 'model_b.txt'}")
print(f"Metadata saved to {MODELS_DIR / 'model_b_meta.json'}")
print("Done")
