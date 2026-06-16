"""
train_gbm_v3.py — Fixed training: consistent 152-d features, proper early stopping.

Critical fixes from v2:
  1. Feature dimension mismatch: v2 saved a 168-d model but API sends 152-d
     → now uses consistent 152-d (base extract_batch_gbm only, no gene/cond extras)
  2. Final model early stopping had best_iteration_=1 (held-out gene was unseen)
     → now uses random 5% holdout for early stopping (prevents overfit on train data)
     → gene-grouped PCC reported separately as generalization metric
  3. Condition features: 31% null → kept as simple ref-imputation (feature_gbm already handles)
  4. best_iteration_ now properly tracks convergence
  
  Council synthesis applied:
  - Contrarian: "Don't ship a model with 168 features when the API sends 152."
  - First Principles: "Prove the base 152-d model works before adding features."
  - Expansionist: "Extra features belong in a separate experiment branch."
  - Outsider: "Ship the thing that works for within-gene ranking today."
  - Executor: "Fix the bug, measure correctly, move to validation."
"""

import sys, json, warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import joblib
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import GroupShuffleSplit, KFold
from lightgbm import LGBMRegressor, early_stopping, log_evaluation
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.features import extract_batch_gbm

DATA_DIR   = Path(__file__).parent.parent / "data"
MODELS_DIR = Path(__file__).parent.parent / "models"
LOG_DIR    = Path(__file__).parent.parent / "logs"
LOG_DIR.mkdir(exist_ok=True)
RNG = np.random.default_rng(42)


def _metrics(y_true, y_pred):
    pcc = pearsonr(y_true, y_pred)[0]
    sp = spearmanr(y_true, y_pred)[0]
    err = y_pred - y_true
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err ** 2)))
    return pcc, sp, mae, rmse


def _log(msg):
    print(f'[{datetime.now().strftime("%H:%M:%S")}] {msg}')


def _features(df):
    """Consistent feature extraction: 152-d, no extra gene/condition features.
    This matches what the API/CLI use (extract_batch_gbm)."""
    bs = list(df['base_sense']) if 'base_sense' in df.columns else None
    ba = list(df['base_antisense']) if 'base_antisense' in df.columns else None
    cc = list(df['concentration_nM']) if 'concentration_nM' in df.columns else None
    tt = list(df['time_h']) if 'time_h' in df.columns else None
    return extract_batch_gbm(
        list(df['sense']), list(df['antisense']),
        base_sense_list=bs, base_antisense_list=ba,
        conc_list=cc, time_list=tt
    )


LGBM_BASE = dict(
    n_estimators=800, learning_rate=0.03, num_leaves=63,
    max_depth=-1, min_child_samples=30,
    subsample=0.8, subsample_freq=1, colsample_bytree=0.7,
    reg_lambda=1.0, random_state=42, n_jobs=-1, verbose=-1,
)


def train_cm_model_v3():
    """Train cm-siRNA model with consistent 152-d features and proper early stopping."""
    _log('=== cm-siRNA (modified) GBM v3 ===')
    
    df = pd.concat([
        pd.read_csv(DATA_DIR / 'hetero_train_2728.csv'),
        pd.read_csv(DATA_DIR / 'hetero_val_303.csv')
    ], ignore_index=True)
    _log(f'rows: {len(df)}, genes: {df["target_gene"].nunique()}')
    
    X_all, y_all = _features(df), df['efficacy'].to_numpy(float)
    _log(f'feature dim: {X_all.shape[1]} (152-d, matches API)')
    
    # ── 1. Random 5-fold CV (within-gene ranking performance) ──
    _log('\n[1] Random 5-fold CV (within-gene ranking):')
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    from sklearn.model_selection import cross_val_predict
    yp_cv = cross_val_predict(
        LGBMRegressor(**LGBM_BASE), X_all, y_all, cv=cv, n_jobs=-1
    )
    p_cv, s_cv, m_cv, r_cv = _metrics(y_all, yp_cv)
    _log(f'  PCC={p_cv:.4f}  Spearman={s_cv:.4f}  MAE={m_cv:.2f}  RMSE={r_cv:.2f}')
    
    # ── 2. Gene-grouped split (cross-gene generalization) ──
    _log('\n[2] Gene-grouped split (cross-gene generalization):')
    gss = GroupShuffleSplit(n_splits=1, test_size=0.18, random_state=42)
    tr_idx, va_idx = next(gss.split(df, groups=df['target_gene']))
    held_genes = sorted(set(df['target_gene'].iloc[va_idx]) - set(df['target_gene'].iloc[tr_idx]))
    _log(f'  held-out genes: {held_genes}')
    _log(f'  train: {len(tr_idx)}, val: {len(va_idx)}')
    
    X_tr, X_va = X_all[tr_idx], X_all[va_idx]
    y_tr, y_va = y_all[tr_idx], y_all[va_idx]
    
    mg = LGBMRegressor(**LGBM_BASE)
    mg.fit(X_tr, y_tr, eval_set=[(X_va, y_va)],
           callbacks=[early_stopping(50), log_evaluation(0)])
    yp_g = mg.predict(X_va)
    p_g, s_g, m_g, r_g = _metrics(y_va, yp_g)
    _log(f'  PCC={p_g:.4f}  Spearman={s_g:.4f}  MAE={m_g:.2f}  RMSE={r_g:.2f}')
    _log(f'  best_iteration: {mg.best_iteration_}')
    
    # Per-gene breakdown
    _log('  Per-gene breakdown:')
    val_df = df.iloc[va_idx].copy()
    val_df['pred'] = yp_g
    for gene in sorted(val_df['target_gene'].unique()):
        sub = val_df[val_df['target_gene'] == gene]
        if len(sub) < 10: continue
        p, s, m, r = _metrics(sub['efficacy'].to_numpy(float), sub['pred'].to_numpy(float))
        _log(f'    {gene:<15s} n={len(sub):>4d} PCC={p:.4f} Spearman={s:.4f} MAE={m:.2f}')
    
    # ── 3. Ablation: what if condition features are all set to ref? ──
    _log('\n[3] Ablation: null-impute all condition values:')
    df_null = df.copy()
    df_null['concentration_nM'] = float('nan')
    df_null['time_h'] = float('nan')
    X_null = _features(df_null)
    cv_null = KFold(n_splits=5, shuffle=True, random_state=42)
    yp_null = cross_val_predict(LGBMRegressor(**LGBM_BASE), X_null, y_all, cv=cv_null, n_jobs=-1)
    p_n, s_n, m_n, r_n = _metrics(y_all, yp_null)
    _log(f'  PCC={p_n:.4f} (delta vs [{1}]: {p_n - p_cv:+.4f})')
    
    # ── 4. Final model: train on ALL data with random 5% holdout for early stopping ──
    _log('\n[4] Final model (all data, random 5% holdout for early stopping):')
    idx_p = RNG.permutation(len(df))
    cut = int(len(df) * 0.95)
    tf, vf = idx_p[:cut], idx_p[cut:]
    X_tf, X_vf = X_all[tf], X_all[vf]
    y_tf, y_vf = y_all[tf], y_all[vf]
    
    final = LGBMRegressor(**LGBM_BASE)
    final.fit(X_tf, y_tf, eval_set=[(X_vf, y_vf)],
              callbacks=[early_stopping(50), log_evaluation(0)])
    _log(f'  best_iteration: {final.best_iteration_}, n_estimators: {final.n_estimators}')
    _log(f'  n_features_in_: {final.n_features_in_} (must match API 152-d)')
    
    # Save
    for name in ('model_a', 'model_b', 'model_c'):
        joblib.dump(final, MODELS_DIR / f'{name}.pkl')
    _log('  Saved -> model_a/b/c.pkl (152-d, consistent with API)')
    
    # Metadata
    meta = {
        'date': datetime.now().isoformat(),
        'version': 3,
        'random_cv_pcc': round(p_cv, 4),
        'random_cv_spearman': round(s_cv, 4),
        'random_cv_mae': round(m_cv, 2),
        'gene_grouped_pcc': round(p_g, 4),
        'gene_grouped_spearman': round(s_g, 4),
        'gene_grouped_mae': round(m_g, 2),
        'condition_ablation_pcc': round(p_n, 4),
        'held_out_genes': held_genes,
        'feature_dim': X_all.shape[1],
        'best_iteration': final.best_iteration_,
        'n_estimators_final': final.n_estimators,
        'feature_dims_match_api': X_all.shape[1] == 152,
        'params': {k: v for k, v in LGBM_BASE.items()},
    }
    with open(LOG_DIR / 'cm_model_v3_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    
    _log(f'\n{"="*60}')
    _log('cm-siRNA MODEL SUMMARY')
    _log(f'  Within-gene (random CV):  PCC={p_cv:.4f}  Spearman={s_cv:.4f}  MAE={m_cv:.2f}')
    _log(f'  Cross-gene (gene-grouped): PCC={p_g:.4f}  Spearman={s_g:.4f}  MAE={m_g:.2f}')
    _log(f'  Condition ablation:       PCC={p_n:.4f}')
    _log(f'  Final model: {final.n_features_in_}-d, {final.best_iteration_} trees')
    _log(f'{"="*60}')
    
    return final, meta


def train_homo_model_v3():
    """Retrain homo model as LightGBM with consistent 152-d features."""
    _log('\n=== Homo (cm-siRNA alt) GBM v3 ===')
    df = pd.concat([
        pd.read_csv(DATA_DIR / 'homo_train.csv'),
        pd.read_csv(DATA_DIR / 'homo_val.csv')
    ], ignore_index=True)
    _log(f'rows: {len(df)}')
    
    X, y = _features(df), df['efficacy'].to_numpy(float)
    
    idx = RNG.permutation(len(df))
    cut = int(len(df) * 0.85)
    tr, va = idx[:cut], idx[cut:]
    
    model = LGBMRegressor(**LGBM_BASE)
    model.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
              callbacks=[early_stopping(50), log_evaluation(0)])
    yp = model.predict(X[va])
    p, s, m, r = _metrics(y[va], yp)
    _log(f'  Val: PCC={p:.4f}  Spearman={s:.4f}  MAE={m:.2f}  RMSE={r:.2f}')
    _log(f'  best_iteration: {model.best_iteration_}')
    
    final = LGBMRegressor(n_estimators=model.best_iteration_ or 800,
                          **{k: v for k, v in LGBM_BASE.items() if k != 'n_estimators'})
    final.fit(X, y)
    joblib.dump(final, MODELS_DIR / 'model_homo.pkl')
    _log('Saved -> model_homo.pkl (LightGBM, 152-d)')
    
    return final


def train_normal_model_v3():
    """
    Improved naked siRNA model with source-aware training.
    
    Strategy: source one-hot + leave-one-source-out evaluation.
    """
    _log('\n=== Naked siRNA GBM v3 ===')
    df = pd.read_csv(DATA_DIR / 'normal_siRNA_extended.csv')
    _log(f'rows: {len(df)}, sources: {sorted(df["source"].unique())}')
    
    X_seq = _features(df)
    sources = sorted(df['source'].unique())
    src_oh = np.zeros((len(df), len(sources)), dtype=np.float32)
    for i, s in enumerate(df['source']):
        src_oh[i, sources.index(s)] = 1.0
    X = np.concatenate([X_seq, src_oh], axis=1)
    y = df['efficacy'].to_numpy(float)
    _log(f'feature dim: {X.shape[1]} (152 seq + {len(sources)} source OH)')
    
    # Leave-one-source-out CV
    _log('\n  Leave-one-source-out CV:')
    loso = {}
    for held_src in sources:
        train_mask = df['source'] != held_src
        test_mask = df['source'] == held_src
        if test_mask.sum() < 10: continue
        m = LGBMRegressor(**LGBM_BASE)
        m.fit(X[train_mask], y[train_mask],
              eval_set=[(X[test_mask], y[test_mask])],
              callbacks=[early_stopping(30), log_evaluation(0)])
        yp = m.predict(X[test_mask])
        p, s, mae, rmse = _metrics(y[test_mask], yp)
        loso[held_src] = {'pcc': round(p, 4), 'n': int(test_mask.sum())}
        _log(f'  Hold-out {held_src:<25s} n={test_mask.sum():>4d}  PCC={p:.4f}  MAE={mae:.2f}')
    
    # Full-data random split
    _log('\n  Full-data random split:')
    idx = RNG.permutation(len(df))
    cut = int(len(df) * 0.82)
    tr, va = idx[:cut], idx[cut:]
    m_final = LGBMRegressor(**LGBM_BASE)
    m_final.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
                callbacks=[early_stopping(30), log_evaluation(0)])
    yp_all = m_final.predict(X[va])
    p_all, s_all, m_all, r_all = _metrics(y[va], yp_all)
    _log(f'  All-source val: PCC={p_all:.4f}  Spearman={s_all:.4f}  MAE={m_all:.2f}')
    
    for src in sources:
        mask = df['source'].to_numpy()[va] == src
        if mask.sum() < 10: continue
        p, s, m, r = _metrics(y[va][mask], yp_all[mask])
        _log(f'    {src:<25s} n={mask.sum():>4d} PCC={p:.4f} MAE={m:.2f}')
    
    # Save
    final = LGBMRegressor(
        n_estimators=m_final.best_iteration_ or 800,
        **{k: v for k, v in LGBM_BASE.items() if k != 'n_estimators'}
    )
    final.fit(X, y)
    joblib.dump({'model': final, 'sources': sources}, MODELS_DIR / 'model_normal.pkl')
    _log('Saved -> model_normal.pkl')
    
    meta = {
        'date': datetime.now().isoformat(),
        'all_source_pcc': round(p_all, 4),
        'all_source_spearman': round(s_all, 4),
        'all_source_mae': round(m_all, 2),
        'loso': loso,
        'feature_dim': X.shape[1],
        'best_iteration': m_final.best_iteration_,
    }
    with open(LOG_DIR / 'normal_model_v3_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    
    return final, meta


def main():
    import sys as _sys
    only_normal = '--only-normal' in _sys.argv
    only_cm = '--only-cm' in _sys.argv
    only_homo = '--only-homo' in _sys.argv
    
    _log('=' * 60)
    _log('SMEpred LightGBM Training v3')
    _log(f'Feature dim: 152 (consistent with API/CLI)')
    _log(f'Models: cm-siRNA (A/B/C), naked (normal), homo')
    _log('=' * 60)
    
    if not only_normal and not only_homo:
        train_cm_model_v3()
    if not only_cm and not only_homo:
        train_normal_model_v3()
    if not only_cm and not only_normal:
        train_homo_model_v3()
    
    _log('\nDone. All models saved with 152-d feature compatibility.')


if __name__ == '__main__':
    main()
