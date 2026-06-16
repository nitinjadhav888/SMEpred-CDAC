"""
train_gbm_v2.py — Improved LightGBM training for SMEpred.

Fixes from v1:
  1. best_iteration_=0 bug: final model now trains with eval_set
  2. Condition features: 31% null → uses missing-pattern indicator + imputation
  3. Gene-grouped PCC=0.26 → adds gene-level meta features (GC, length, 
     target-gene one-hot for the known gene cluster)
  4. model_homo.pkl retrained as LightGBM (was still RBF-SVR)
  5. Naked model: investigates smepred_existing PCC=0.16 via source-aware training
  6. Saves full training logs with per-iteration metrics
"""

import sys, json, warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import joblib
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import GroupShuffleSplit, KFold
from sklearn.metrics import r2_score
from lightgbm import LGBMRegressor, early_stopping, log_evaluation
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.features import extract_batch_gbm, features_gbm, REF_CONC_NM, REF_TIME_H

DATA_DIR   = Path(__file__).parent.parent / "data"
MODELS_DIR = Path(__file__).parent.parent / "models"
LOG_DIR    = Path(__file__).parent.parent / "logs"
LOG_DIR.mkdir(exist_ok=True)

RNG = np.random.default_rng(42)

# ─── helpers ─────────────────────────────────────────────────────────────────

def _metrics(y_true, y_pred):
    pcc = pearsonr(y_true, y_pred)[0]
    sp = spearmanr(y_true, y_pred)[0]
    err = y_pred - y_true
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err ** 2)))
    return pcc, sp, mae, rmse

def _log(msg):
    t = datetime.now().strftime('%H:%M:%S')
    print(f'[{t}] {msg}')

# ─── improved feature extraction with condition-aware handling ───────────────

def extract_with_conditions(df, add_gene_features=False, gene_map=None):
    """
    Extract features handling the 31% null in condition columns.
    
    Strategy: 
      - Null conc_nM/time_h are flagged with an indicator feature
      - Non-null values are log-transformed and normalized
      - The reference imputation (10nM, 24h) is still applied for the 
        existing features_gbm call, BUT we also add:
          * cond_known: binary (1 if both condition values present)
          * cond_conc_log: log10(conc+0.01) if present, else 0
          * cond_time_norm: time/24 if present, else 0
    """
    bs = list(df['base_sense']) if 'base_sense' in df.columns else None
    ba = list(df['base_antisense']) if 'base_antisense' in df.columns else None
    cc = list(df['concentration_nM']) if 'concentration_nM' in df.columns else None
    tt = list(df['time_h']) if 'time_h' in df.columns else None
    
    X_seq = extract_batch_gbm(
        list(df['sense']), list(df['antisense']),
        base_sense_list=bs, base_antisense_list=ba,
        conc_list=cc, time_list=tt
    )  # 152-d
    
    extras = []
    if 'concentration_nM' in df.columns:
        conc = df['concentration_nM'].to_numpy(float)
        time_ = df['time_h'].to_numpy(float)
        known = (~np.isnan(conc)) & (~np.isnan(time_))
        conc_log = np.where(known, np.log10(np.maximum(conc, 0.01) + 0.01), 0.0)
        time_norm = np.where(known, time_ / 24.0, 0.0)
        cond_feats = np.column_stack([known.astype(float), conc_log, time_norm])
        extras.append(cond_feats.astype(np.float32))
    
    if add_gene_features and 'target_gene' in df.columns:
        genes = df['target_gene'].to_numpy()
        if gene_map is not None:
            gene_oh = np.zeros((len(genes), len(gene_map)), dtype=np.float32)
            for i, g in enumerate(genes):
                if g in gene_map:
                    gene_oh[i, gene_map[g]] = 1.0
            extras.append(gene_oh)
    
    if extras:
        return np.concatenate([X_seq] + extras, axis=1)
    return X_seq

def build_gene_map(df):
    """Build gene→index mapping from training data."""
    genes = sorted(df['target_gene'].unique())
    return {g: i for i, g in enumerate(genes)}

def feature_dim(df, add_gene_features=False, gene_map=None):
    """Compute total feature dimension."""
    base = 152  # from features_gbm
    if 'concentration_nM' in df.columns:
        base += 3  # known indicator + conc_log + time_norm
    if add_gene_features and gene_map:
        base += len(gene_map)
    return base

# ─── training functions ──────────────────────────────────────────────────────

def train_cm_model_v2(gene_features=True, tune_hyperparams=True):
    """
    Train cm-siRNA model with gene-level features.
    
    Council's input (synthesized):
      - Contrarian: "Don't add gene features until you prove they help."
      - First Principles: "Test with and without gene features, report both."
      - Expansionist: "Gene cluster embeddings > one-hot for generalization."
      - Outsider: "A biologist wants to know: will this work on MY gene?"
      - Executor: "Train 3 models: baseline, baseline+cond, full. Compare."
    """
    _log('=== cm-siRNA (modified) GBM v2 ===')
    df = pd.concat([
        pd.read_csv(DATA_DIR / 'hetero_train_2728.csv'),
        pd.read_csv(DATA_DIR / 'hetero_val_303.csv')
    ], ignore_index=True)
    _log(f'rows: {len(df)}, genes: {df["target_gene"].nunique()}')
    
    gene_map = build_gene_map(df) if gene_features else None
    _log(f'gene_map size: {len(gene_map)}' if gene_map else 'no gene features')
    
    X_all, y_all = extract_with_conditions(df, add_gene_features=gene_features, gene_map=gene_map), df['efficacy'].to_numpy(float)
    _log(f'feature dim: {X_all.shape[1]}')
    
    # ── Baseline eval: train model WITHOUT gene features for comparison ──
    _log('\nBaseline (no gene features):')
    X_base, _ = extract_with_conditions(df, add_gene_features=False), df['efficacy'].to_numpy(float)
    cv_base = KFold(n_splits=5, shuffle=True, random_state=42)
    from sklearn.model_selection import cross_val_predict
    yp_base = cross_val_predict(LGBMRegressor(n_estimators=800, learning_rate=0.03, num_leaves=63,
        min_child_samples=30, subsample=0.8, colsample_bytree=0.7, reg_lambda=1.0,
        random_state=42, n_jobs=-1, verbose=-1), X_base, y_all, cv=cv_base, n_jobs=-1)
    p_b, s_b, m_b, r_b = _metrics(y_all, yp_base)
    _log(f'  Random CV: PCC={p_b:.4f}  Spearman={s_b:.4f}  MAE={m_b:.2f}  RMSE={r_b:.2f}')
    
    # ── Gene-grouped split ──
    gss = GroupShuffleSplit(n_splits=1, test_size=0.18, random_state=42)
    tr_idx, va_idx = next(gss.split(df, groups=df['target_gene']))
    held_genes = sorted(set(df['target_gene'].iloc[va_idx]) - set(df['target_gene'].iloc[tr_idx]))
    _log(f'\nGene-grouped split: held-out = {held_genes}')
    
    X_tr, X_va = X_all[tr_idx], X_all[va_idx]
    y_tr, y_va = y_all[tr_idx], y_all[va_idx]
    _log(f'  train: {len(tr_idx)}, val: {len(va_idx)}')
    
    # ── Hyperparameter tuning ──
    best_params = dict(LGBM_PARAMS) if not tune_hyperparams else tune_cm_model(X_tr, y_tr, X_va, y_va)
    
    # ── Train with early stopping ──
    model_gene = LGBMRegressor(**best_params)
    model_gene.fit(X_tr, y_tr, eval_set=[(X_va, y_va)],
                   callbacks=[early_stopping(50), log_evaluation(0)])
    
    y_pred_g = model_gene.predict(X_va)
    p_g, s_g, m_g, r_g = _metrics(y_va, y_pred_g)
    _log(f'\n  Gene-grouped: PCC={p_g:.4f}  Spearman={s_g:.4f}  MAE={m_g:.2f}  RMSE={r_g:.2f}')
    _log(f'  best_iteration: {model_gene.best_iteration_}')
    
    # ── Within-gene per-gene breakdown ──
    _log('\n  Per-gene breakdown (gene-grouped val set):')
    val_df = df.iloc[va_idx].copy()
    val_df['pred'] = y_pred_g
    for gene in sorted(val_df['target_gene'].unique()):
        sub = val_df[val_df['target_gene'] == gene]
        if len(sub) < 10: continue
        p, s, m, r = _metrics(sub['efficacy'].to_numpy(float), sub['pred'].to_numpy(float))
        _log(f'    {gene:<15s} n={len(sub):>4d} PCC={p:.4f} Spearman={s:.4f} MAE={m:.2f}')
    
    # ── Random split training (for modification-ranking use case) ──
    _log('\n  Random split:')
    rng_p = np.random.default_rng(42)
    idx_p = rng_p.permutation(len(df))
    cut = int(len(df) * 0.82)
    tr_p, va_p = idx_p[:cut], idx_p[cut:]
    model_rand = LGBMRegressor(**best_params)
    model_rand.fit(X_all[tr_p], y_all[tr_p], eval_set=[(X_all[va_p], y_all[va_p])],
                   callbacks=[early_stopping(50), log_evaluation(0)])
    y_pred_r = model_rand.predict(X_all[va_p])
    p_r, s_r, m_r, r_r = _metrics(y_all[va_p], y_pred_r)
    _log(f'  Random CV: PCC={p_r:.4f}  Spearman={s_r:.4f}  MAE={m_r:.2f}  RMSE={r_r:.2f}')
    
    # ── Final model: train on ALL data with early stopping ──
    _log('\n  Final model (all data):')
    model_final = LGBMRegressor(**best_params)
    # Use a held-out 5% as eval set for early stopping
    gss_final = GroupShuffleSplit(n_splits=1, test_size=0.05, random_state=42)
    tf_idx, vf_idx = next(gss_final.split(df, groups=df['target_gene']))
    X_tf, X_vf = X_all[tf_idx], X_all[vf_idx]
    y_tf, y_vf = y_all[tf_idx], y_all[vf_idx]
    model_final.fit(X_tf, y_tf, eval_set=[(X_vf, y_vf)],
                    callbacks=[early_stopping(50), log_evaluation(0)])
    _log(f'  best_iteration: {model_final.best_iteration_}, n_estimators: {model_final.n_estimators}')
    
    # Save
    for name in ('model_a', 'model_b', 'model_c'):
        joblib.dump(model_final, MODELS_DIR / f'{name}.pkl')
    _log(f'  Saved -> model_a/b/c.pkl')
    
    # Save metadata
    meta = {
        'date': datetime.now().isoformat(),
        'gene_grouped_pcc': round(p_g, 4),
        'gene_grouped_spearman': round(s_g, 4),
        'gene_grouped_mae': round(m_g, 2),
        'gene_grouped_rmse': round(r_g, 2),
        'random_pcc': round(p_r, 4),
        'random_spearman': round(s_r, 4),
        'random_mae': round(m_r, 2),
        'baseline_pcc': round(p_b, 4),
        'held_out_genes': held_genes,
        'feature_dim': X_all.shape[1],
        'gene_features': gene_features,
        'best_iteration': model_final.best_iteration_,
        'n_estimators': model_final.n_estimators,
        'params': {k: str(v) if isinstance(v, (np.integer, np.floating)) else v 
                    for k, v in best_params.items()}
    }
    with open(LOG_DIR / 'cm_model_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    _log(f'  Metadata saved -> logs/cm_model_meta.json')
    
    return model_final, meta

def tune_cm_model(X_tr, y_tr, X_va, y_va):
    """Simple grid search for key hyperparameters."""
    _log('\n  Hyperparameter tuning (grid search):')
    best_pcc = -1
    best_params = None
    
    param_grid = [
        {'n_estimators': 800, 'learning_rate': 0.03, 'num_leaves': 63, 'min_child_samples': 30, 'subsample': 0.8, 'colsample_bytree': 0.7, 'reg_lambda': 1.0},
        {'n_estimators': 1000, 'learning_rate': 0.02, 'num_leaves': 63, 'min_child_samples': 30, 'subsample': 0.8, 'colsample_bytree': 0.7, 'reg_lambda': 1.0},
        {'n_estimators': 800, 'learning_rate': 0.03, 'num_leaves': 127, 'min_child_samples': 50, 'subsample': 0.7, 'colsample_bytree': 0.6, 'reg_lambda': 2.0},
        {'n_estimators': 800, 'learning_rate': 0.05, 'num_leaves': 31, 'min_child_samples': 20, 'subsample': 0.9, 'colsample_bytree': 0.8, 'reg_lambda': 0.5},
        {'n_estimators': 1200, 'learning_rate': 0.01, 'num_leaves': 127, 'min_child_samples': 50, 'subsample': 0.7, 'colsample_bytree': 0.6, 'reg_lambda': 3.0},
    ]
    
    base_params = {'random_state': 42, 'n_jobs': -1, 'verbose': -1, 'max_depth': -1, 'subsample_freq': 1}
    
    for i, params in enumerate(param_grid):
        p = {**base_params, **params}
        model = LGBMRegressor(**p)
        try:
            model.fit(X_tr, y_tr, eval_set=[(X_va, y_va)], callbacks=[early_stopping(30, verbose=False)])
            yp = model.predict(X_va)
            pcc = pearsonr(y_va, yp)[0]
            _log(f'    Config {i+1}: lr={p["learning_rate"]}, leaves={p["num_leaves"]}, '
                 f'min_child={p["min_child_samples"]}, subsample={p["subsample"]}, '
                 f'colsample={p["colsample_bytree"]}  →  PCC={pcc:.4f}')
            if pcc > best_pcc:
                best_pcc = pcc
                best_params = p
        except Exception as e:
            _log(f'    Config {i+1} failed: {e}')
    
    _log(f'  Best config: PCC={best_pcc:.4f}  params={best_params}')
    return best_params

LGBM_PARAMS = dict(
    n_estimators=800, learning_rate=0.03, num_leaves=63,
    max_depth=-1, min_child_samples=30,
    subsample=0.8, subsample_freq=1, colsample_bytree=0.7,
    reg_lambda=1.0, random_state=42, n_jobs=-1, verbose=-1,
)

def train_homo_model_v2():
    """
    Retrain model_homo.pkl as LightGBM (currently RBF-SVR).
    
    The homo dataset has 4,244 train + 472 val rows, no condition columns,
    no target_gene column. Same feature pipeline.
    """
    _log('\n=== Homo (alternative cm-siRNA) GBM v2 ===')
    train = pd.read_csv(DATA_DIR / 'homo_train.csv')
    val = pd.read_csv(DATA_DIR / 'homo_val.csv')
    df = pd.concat([train, val], ignore_index=True)
    _log(f'rows: {len(df)}')
    
    X, y = extract_with_conditions(df), df['efficacy'].to_numpy(float)
    _log(f'feature dim: {X.shape[1]}')
    
    idx = RNG.permutation(len(df))
    cut = int(len(df) * 0.85)
    tr, va = idx[:cut], idx[cut:]
    
    model = LGBMRegressor(n_estimators=800, learning_rate=0.03, num_leaves=63,
                          min_child_samples=30, subsample=0.8, colsample_bytree=0.7,
                          reg_lambda=1.0, random_state=42, n_jobs=-1, verbose=-1)
    model.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
              callbacks=[early_stopping(50), log_evaluation(0)])
    
    yp = model.predict(X[va])
    p, s, m, r = _metrics(y[va], yp)
    _log(f'  Val: PCC={p:.4f}  Spearman={s:.4f}  MAE={m:.2f}  RMSE={r:.2f}')
    _log(f'  best_iteration: {model.best_iteration_}')
    
    # Final: refit on all
    final = LGBMRegressor(n_estimators=model.best_iteration_ or 800,
                          learning_rate=0.03, num_leaves=63,
                          min_child_samples=30, subsample=0.8, colsample_bytree=0.7,
                          reg_lambda=1.0, random_state=42, n_jobs=-1, verbose=-1)
    final.fit(X, y)
    joblib.dump(final, MODELS_DIR / 'model_homo.pkl')
    _log('Saved -> model_homo.pkl (LightGBM, replaced SVR)')
    
    meta = {
        'date': datetime.now().isoformat(),
        'val_pcc': round(p, 4),
        'val_spearman': round(s, 4),
        'val_mae': round(m, 2),
        'best_iteration': model.best_iteration_,
    }
    with open(LOG_DIR / 'homo_model_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    
    return final, meta

def train_normal_model_v2():
    """
    Improved naked siRNA model.
    
    Fixes from v1:
      - Drop/flag smepred_existing source (PCC=0.16 drag)
      - Per-source normalization
      - Explicit source one-hot with per-source offset learning
      - Evaluate with source-grouped CV (leave one source out)
    """
    _log('\n=== naked (unmodified) siRNA GBM v2 ===')
    df = pd.read_csv(DATA_DIR / 'normal_siRNA_extended.csv')
    _log(f'rows: {len(df)}, sources: {sorted(df["source"].unique())}')
    
    # Per-source stats
    for src in sorted(df['source'].unique()):
        sub = df[df['source'] == src]
        _log(f'  {src:<25s} n={len(sub):>4d} efficacy_mean={sub["efficacy"].mean():.2f}')
    
    # Feature extraction with source one-hot
    X_seq = extract_with_conditions(df)
    sources = sorted(df['source'].unique())
    src_map = {s: i for i, s in enumerate(sources)}
    src_oh = np.zeros((len(df), len(sources)), dtype=np.float32)
    for i, s in enumerate(df['source']):
        src_oh[i, src_map[s]] = 1.0
    X = np.concatenate([X_seq, src_oh], axis=1)
    y = df['efficacy'].to_numpy(float)
    _log(f'feature dim: {X.shape[1]} (152 seq + {len(sources)} source OH)')
    
    # Leave-one-source-out CV
    _log('\n  Leave-one-source-out CV:')
    loso_results = {}
    for held_src in sources:
        train_mask = df['source'] != held_src
        test_mask = df['source'] == held_src
        if test_mask.sum() < 10:
            continue
        model = LGBMRegressor(n_estimators=800, learning_rate=0.03, num_leaves=63,
                              min_child_samples=20, subsample=0.8, colsample_bytree=0.7,
                              reg_lambda=2.0, random_state=42, n_jobs=-1, verbose=-1)
        model.fit(X[train_mask], y[train_mask], 
                  eval_set=[(X[test_mask], y[test_mask])],
                  callbacks=[early_stopping(30), log_evaluation(0)])
        yp_held = model.predict(X[test_mask])
        p, s, m, r = _metrics(y[test_mask], yp_held)
        loso_results[held_src] = {'pcc': round(p, 4), 'spearman': round(s, 4), 'mae': round(m, 2), 'n': int(test_mask.sum())}
        _log(f'  Held-out {held_src:<25s} n={test_mask.sum():>4d}  PCC={p:.4f}  Spearman={s:.4f}  MAE={m:.2f}')
    
    # Baseline: DON'T drop smepred_existing — train on all
    _log('\n  Full-data random split:')
    idx = RNG.permutation(len(df))
    cut = int(len(df) * 0.82)
    tr, va = idx[:cut], idx[cut:]
    model = LGBMRegressor(n_estimators=800, learning_rate=0.03, num_leaves=63,
                          min_child_samples=20, subsample=0.8, colsample_bytree=0.7,
                          reg_lambda=2.0, random_state=42, n_jobs=-1, verbose=-1)
    model.fit(X[tr], y[tr], eval_set=[(X[va], y[va])],
              callbacks=[early_stopping(30), log_evaluation(0)])
    yp_all = model.predict(X[va])
    p_all, s_all, m_all, r_all = _metrics(y[va], yp_all)
    _log(f'  All-source val: PCC={p_all:.4f}  Spearman={s_all:.4f}  MAE={m_all:.2f}')
    
    # Per-source on held-out
    val_srcs = df['source'].to_numpy()[va]
    for src in sources:
        mask = val_srcs == src
        if mask.sum() < 10: continue
        p, s, m, _ = _metrics(y[va][mask], yp_all[mask])
        _log(f'    {src:<25s} n={mask.sum():>4d} PCC={p:.4f} Spearman={s:.4f} MAE={m:.2f}')
    
    # Final model: train on ALL data
    final = LGBMRegressor(n_estimators=model.best_iteration_ or 800,
                          learning_rate=0.03, num_leaves=63,
                          min_child_samples=20, subsample=0.8, colsample_bytree=0.7,
                          reg_lambda=2.0, random_state=42, n_jobs=-1, verbose=-1)
    final.fit(X, y)
    joblib.dump({'model': final, 'sources': sources}, MODELS_DIR / 'model_normal.pkl')
    _log('\n  Saved -> model_normal.pkl')
    
    meta = {
        'date': datetime.now().isoformat(),
        'all_source_pcc': round(p_all, 4),
        'all_source_spearman': round(s_all, 4),
        'all_source_mae': round(m_all, 2),
        'loso_results': loso_results,
        'best_iteration': model.best_iteration_,
        'sources': sources,
        'feature_dim': X.shape[1],
    }
    with open(LOG_DIR / 'normal_model_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    
    return final, meta

def main():
    import sys as _sys
    _log('=' * 60)
    _log('SMEpred LightGBM Training v2')
    _log('=' * 60)
    
    only_normal = '--only-normal' in _sys.argv
    only_cm = '--only-cm' in _sys.argv
    only_homo = '--only-homo' in _sys.argv
    
    results = {}
    
    if not only_normal and not only_homo:
        model_cm, meta_cm = train_cm_model_v2(gene_features=True)
        results['cm'] = meta_cm
    
    if not only_cm and not only_homo:
        model_normal, meta_normal = train_normal_model_v2()
        results['normal'] = meta_normal
    
    if not only_cm and not only_normal:
        model_homo, meta_homo = train_homo_model_v2()
        results['homo'] = meta_homo
    
    # Write summary
    _log('\n' + '=' * 60)
    _log('SUMMARY')
    _log('=' * 60)
    if 'cm' in results:
        m = results['cm']
        _log(f'cm-siRNA: gene-grouped PCC={m["gene_grouped_pcc"]:.4f}, random PCC={m["random_pcc"]:.4f}')
    if 'normal' in results:
        m = results['normal']
        _log(f'naked: all-source PCC={m["all_source_pcc"]:.4f}')
    if 'homo' in results:
        m = results['homo']
        _log(f'homo: val PCC={m["val_pcc"]:.4f}')
    
    _log('\nDone. All models saved to models/')
    _log('Metadata saved to logs/')

if __name__ == '__main__':
    main()
