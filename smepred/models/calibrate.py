"""
calibrate.py — Fit isotonic calibrators to map raw model scores → calibrated 0-100.

Problem: LightGBM trained on noisy, mean-centered labels (47.7) produces predictions
squeezed into a narrow range (~30-70). Isotonic regression learns a monotonic mapping
from raw scores to the true 0-100 scale using a held-out validation set.

Two calibrators are fitted:
  1. cm-siRNA calibrator (for model_a/b/c)   — using hetero_val_303.csv
  2. naked siRNA calibrator (model_normal)    — using normal_siRNA_extended.csv split
"""

import sys, json, joblib, warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import train_test_split
from lightgbm import LGBMRegressor

warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.features import extract_batch_gbm

DATA_DIR = Path(__file__).parent.parent / 'data'
MODELS_DIR = Path(__file__).parent.parent / 'models'
LOG_DIR = Path(__file__).parent.parent / 'logs'
LOG_DIR.mkdir(exist_ok=True)

RNG = np.random.default_rng(42)

def _features(df):
    bs = list(df['base_sense']) if 'base_sense' in df.columns else None
    ba = list(df['base_antisense']) if 'base_antisense' in df.columns else None
    cc = list(df['concentration_nM']) if 'concentration_nM' in df.columns else None
    tt = list(df['time_h']) if 'time_h' in df.columns else None
    return extract_batch_gbm(list(df['sense']), list(df['antisense']),
                             base_sense_list=bs, base_antisense_list=ba,
                             conc_list=cc, time_list=tt)


def fit_cm_calibrator():
    """Fit isotonic calibrator for the cm-siRNA model (model_a/b/c)."""
    print('=== cm-siRNA Calibrator ===')
    
    df = pd.concat([
        pd.read_csv(DATA_DIR / 'hetero_train_2728.csv'),
        pd.read_csv(DATA_DIR / 'hetero_val_303.csv')
    ], ignore_index=True)
    print(f'Total rows: {len(df)}')
    
    X = _features(df)
    y = df['efficacy'].to_numpy(float)
    
    # Load deployed model
    model = joblib.load(MODELS_DIR / 'model_a.pkl')
    
    # Get held-out predictions via 5-fold CV to avoid overfitting the calibrator
    from sklearn.model_selection import KFold, cross_val_predict
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    raw_preds = cross_val_predict(model, X, y, cv=cv, n_jobs=-1)
    
    # Fit isotonic regression: raw_preds → true labels
    calibrator = IsotonicRegression(out_of_bounds='clip', y_min=0, y_max=100)
    calibrator.fit(raw_preds, y)
    
    # Evaluate
    calibrated = calibrator.transform(raw_preds)
    pcc_raw = pearsonr(y, raw_preds)[0]
    pcc_cal = pearsonr(y, calibrated)[0]
    mae_raw = float(np.mean(np.abs(raw_preds - y)))
    mae_cal = float(np.mean(np.abs(calibrated - y)))
    print(f'  Raw:       PCC={pcc_raw:.4f}  MAE={mae_raw:.2f}  range=[{raw_preds.min():.1f}, {raw_preds.max():.1f}]')
    print(f'  Calibrated: PCC={pcc_cal:.4f}  MAE={mae_cal:.2f}  range=[{calibrated.min():.1f}, {calibrated.max():.1f}]')
    print(f'  Calibrator thresholds: {len(calibrator.X_thresholds_)} points')
    print(f'  Raw input range seen: [{raw_preds.min():.1f}, {raw_preds.max():.1f}]')
    
    # Save
    joblib.dump(calibrator, MODELS_DIR / 'calibrator_cm.pkl')
    print('  Saved -> calibrator_cm.pkl')
    
    result = {
        'raw_pcc': round(pcc_raw, 4), 'cal_pcc': round(pcc_cal, 4),
        'raw_mae': round(mae_raw, 2), 'cal_mae': round(mae_cal, 2),
        'raw_range': [round(float(raw_preds.min()), 1), round(float(raw_preds.max()), 1)],
        'cal_range': [round(float(calibrated.min()), 1), round(float(calibrated.max()), 1)],
    }
    return calibrator, result


def fit_naked_calibrator():
    """Fit isotonic calibrator for the naked siRNA model (model_normal)."""
    print('\n=== Naked siRNA Calibrator ===')
    
    df = pd.read_csv(DATA_DIR / 'normal_siRNA_extended.csv')
    print(f'Total rows: {len(df)}')
    
    X_seq = _features(df)
    sources = sorted(df['source'].unique())
    src_oh = np.zeros((len(df), len(sources)), dtype=np.float32)
    for i, s in enumerate(df['source']):
        src_oh[i, sources.index(s)] = 1.0
    X = np.concatenate([X_seq, src_oh], axis=1)
    y = df['efficacy'].to_numpy(float)
    
    bundle = joblib.load(MODELS_DIR / 'model_normal.pkl')
    model = bundle['model']
    
    # Held-out predictions via CV
    from sklearn.model_selection import KFold, cross_val_predict
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    raw_preds = cross_val_predict(model, X, y, cv=cv, n_jobs=-1)
    
    calibrator = IsotonicRegression(out_of_bounds='clip', y_min=0, y_max=100)
    calibrator.fit(raw_preds, y)
    
    calibrated = calibrator.transform(raw_preds)
    pcc_raw = pearsonr(y, raw_preds)[0]
    pcc_cal = pearsonr(y, calibrated)[0]
    mae_raw = float(np.mean(np.abs(raw_preds - y)))
    mae_cal = float(np.mean(np.abs(calibrated - y)))
    print(f'  Raw:       PCC={pcc_raw:.4f}  MAE={mae_raw:.2f}  range=[{raw_preds.min():.1f}, {raw_preds.max():.1f}]')
    print(f'  Calibrated: PCC={pcc_cal:.4f}  MAE={mae_cal:.2f}  range=[{calibrated.min():.1f}, {calibrated.max():.1f}]')
    
    joblib.dump(calibrator, MODELS_DIR / 'calibrator_naked.pkl')
    print('  Saved -> calibrator_naked.pkl')
    
    result = {
        'raw_pcc': round(pcc_raw, 4), 'cal_pcc': round(pcc_cal, 4),
        'raw_mae': round(mae_raw, 2), 'cal_mae': round(mae_cal, 2),
        'raw_range': [round(float(raw_preds.min()), 1), round(float(raw_preds.max()), 1)],
        'cal_range': [round(float(calibrated.min()), 1), round(float(calibrated.max()), 1)],
    }
    return calibrator, result


def main():
    print('=' * 60)
    print('Score Calibration Pipeline')
    print('=' * 60)
    
    cal_cm, res_cm = fit_cm_calibrator()
    cal_naked, res_naked = fit_naked_calibrator()
    
    print('\n' + '=' * 60)
    print('CALIBRATION SUMMARY')
    print('=' * 60)
    print(f'cm-siRNA:')
    print(f'  Raw PCC={res_cm["raw_pcc"]:.4f} → Cal PCC={res_cm["cal_pcc"]:.4f}')
    print(f'  Raw MAE={res_cm["raw_mae"]:.2f} → Cal MAE={res_cm["cal_mae"]:.2f}')
    print(f'  Raw range=[{res_cm["raw_range"][0]:.1f}, {res_cm["raw_range"][1]:.1f}] → Cal range=[{res_cm["cal_range"][0]:.1f}, {res_cm["cal_range"][1]:.1f}]')
    print(f'naked siRNA:')
    print(f'  Raw PCC={res_naked["raw_pcc"]:.4f} → Cal PCC={res_naked["cal_pcc"]:.4f}')
    print(f'  Raw MAE={res_naked["raw_mae"]:.2f} → Cal MAE={res_naked["cal_mae"]:.2f}')
    print(f'  Raw range=[{res_naked["raw_range"][0]:.1f}, {res_naked["raw_range"][1]:.1f}] → Cal range=[{res_naked["cal_range"][0]:.1f}, {res_naked["cal_range"][1]:.1f}]')
    
    all_res = {'cm': res_cm, 'naked': res_naked}
    with open(LOG_DIR / 'calibration_meta.json', 'w') as f:
        json.dump(all_res, f, indent=2)
    print('\nMetadata saved to logs/calibration_meta.json')


if __name__ == '__main__':
    main()
