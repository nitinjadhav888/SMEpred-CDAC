"""
visualize_lgbm.py — produce validation curve, learning curve, permutation test,
and cross-validation scores for the deployed LightGBM model.

Usage:
    python -m smepred.tools.visualize_lgbm --model models/model_a.pkl --data data/hetero_train_2728.csv --outdir figures

The script is defensive: it checks files and dependencies and writes PNGs into `--outdir`.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

def _ensure_modules():
    missing = []
    try:
        import numpy as np
        import pandas as pd
        import joblib
        import matplotlib
        import matplotlib.pyplot as plt
        from sklearn.model_selection import validation_curve, learning_curve, cross_val_score, GroupShuffleSplit
        from sklearn.model_selection import KFold
        from sklearn.inspection import permutation_importance
    except Exception as e:
        missing.append(str(e))
    if missing:
        print("Missing or failing imports:")
        for m in missing:
            print(" -", m)
        print("Install dependencies: pip install -r requirements.txt")
        sys.exit(1)

def load_data(path: Path):
    import pandas as pd
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")
    return pd.read_csv(path)

def main():
    _ensure_modules()
    import numpy as np
    import pandas as pd
    import joblib
    import matplotlib.pyplot as plt
    from sklearn.model_selection import validation_curve, learning_curve, cross_val_score, GroupShuffleSplit, KFold
    from sklearn.metrics import r2_score

    p = argparse.ArgumentParser()
    p.add_argument('--model', type=Path, default=Path('models/model_a.pkl'))
    p.add_argument('--data', type=Path, default=Path('data/hetero_train_2728.csv'))
    p.add_argument('--outdir', type=Path, default=Path('figures'))
    p.add_argument('--cv', type=int, default=5)
    args = p.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if not args.model.exists():
        print(f'Model file not found: {args.model}')
        sys.exit(1)
    print('Loading model:', args.model)
    model = joblib.load(args.model)

    print('Loading data:', args.data)
    df = load_data(args.data)
    if 'sense' in df.columns and 'antisense' in df.columns and 'efficacy' in df.columns:
        X_s = list(df['sense'])
        X_a = list(df['antisense'])
        y = df['efficacy'].to_numpy(dtype=float)
        # We need to extract features using the project's feature extractor
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from src.features import extract_batch_gbm
        X = extract_batch_gbm(X_s, X_a)
    else:
        # try wide-format numeric features
        df_num = df.select_dtypes(include=[float, int])
        if df_num.shape[1] == 0:
            raise ValueError('No recognizable features in data file.')
        X = df_num.to_numpy(dtype=float)
        y = df['efficacy'].to_numpy(dtype=float)

    print('Data shape:', X.shape)

    # If the saved model knows feature names, wrap X in a DataFrame so LightGBM
    # doesn't warn about missing feature names and cross_val respects the names.
    try:
        fnames = None
        # sklearn wrapper (LGBMRegressor) exposes feature_name_
        fnames = getattr(model, 'feature_name_', None)
        # fallback: if a dict was saved with {'model': ...}
        if fnames is None and isinstance(model, dict) and 'model' in model:
            fnames = getattr(model['model'], 'feature_name_', None)
        if fnames is not None and len(fnames) == X.shape[1]:
            import pandas as _pd
            X = _pd.DataFrame(X, columns=fnames)
            print('Wrapped X as DataFrame with feature names to match model')
    except Exception:
        pass

    # Cross-validation scores
    print('Computing cross-validation R^2 and PCC/Spearman scores...')
    cv = KFold(n_splits=args.cv, shuffle=True, random_state=42)
    try:
        from sklearn.model_selection import cross_val_predict
        from scipy.stats import pearsonr, spearmanr
        y_pred = cross_val_predict(model, X, y, cv=cv, n_jobs=-1)
        from sklearn.metrics import r2_score
        r2 = r2_score(y, y_pred)
        pcc = pearsonr(y, y_pred)[0]
        sp = spearmanr(y, y_pred)[0]
        print(f'  CV (all folds aggregated)  R^2={r2:.4f}  PCC={pcc:.4f}  Spearman={sp:.4f}')
        # also per-fold R2 via cross_val_score for plotting
        scores = cross_val_score(model, X, y, cv=cv, scoring='r2', n_jobs=-1)
        print('  Per-fold R^2 scores:', scores)
        import numpy as np
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(range(1, len(scores) + 1), scores, marker='o')
        plt.title('Cross-validation R^2 (KFold)')
        plt.xlabel('Fold')
        plt.ylabel('R^2')
        plt.grid(True)
        plt.savefig(outdir / 'cv_r2.png', dpi=150)
        plt.close()
    except Exception as e:
        print('cross_val evaluation failed:', e)
        scores = None

    # Learning curve
    print('Computing learning curve...')
    try:
        train_sizes, train_scores, test_scores = learning_curve(
            model, X, y, cv=cv, scoring='r2', n_jobs=-1,
            train_sizes=np.linspace(0.1, 1.0, 5), shuffle=True, random_state=42,
        )
        train_scores_mean = train_scores.mean(axis=1)
        test_scores_mean = test_scores.mean(axis=1)
        plt.figure()
        plt.plot(train_sizes, train_scores_mean, 'o-', label='Train R2')
        plt.plot(train_sizes, test_scores_mean, 'o-', label='Validation R2')
        plt.title('Learning Curve (R^2)')
        plt.xlabel('Training set size')
        plt.ylabel('R^2')
        plt.legend()
        plt.grid(True)
        plt.savefig(outdir / 'learning_curve_r2.png', dpi=150)
        plt.close()
    except Exception as e:
        print('learning_curve failed:', e)

    # Validation curve (vary n_estimators)
    print('Computing validation curve (n_estimators)...')
    try:
        param_name = 'n_estimators'
        param_range = [50, 100, 200, 400, 800]
        train_scores_v, test_scores_v = validation_curve(
            model, X, y, param_name=param_name, param_range=param_range,
            scoring='r2', cv=cv, n_jobs=-1,
        )
        import numpy as np
        train_mean = train_scores_v.mean(axis=1)
        test_mean = test_scores_v.mean(axis=1)
        plt.figure()
        plt.plot(param_range, train_mean, 'o-', label='Train R2')
        plt.plot(param_range, test_mean, 'o-', label='Validation R2')
        plt.title('Validation Curve (n_estimators)')
        plt.xlabel('n_estimators')
        plt.ylabel('R^2')
        plt.xscale('log')
        plt.legend()
        plt.grid(True)
        plt.savefig(outdir / 'validation_curve_n_estimators.png', dpi=150)
        plt.close()
    except Exception as e:
        print('validation_curve failed:', e)

    # Permutation importance (global)
    print('Computing permutation importance (this can be slow)...')
    try:
        from sklearn.inspection import permutation_importance
        perm = permutation_importance(model, X, y, n_repeats=10, random_state=42, n_jobs=-1)
        importances = perm.importances_mean
        idx = importances.argsort()[::-1][:30]
        plt.figure(figsize=(8,6))
        plt.bar(range(len(idx)), importances[idx])
        plt.xticks(range(len(idx)), [str(i) for i in idx], rotation=90)
        plt.title('Permutation importance (top features)')
        plt.tight_layout()
        plt.savefig(outdir / 'permutation_importance.png', dpi=150)
        plt.close()
    except Exception as e:
        print('permutation_importance failed:', e)

    print('Plots saved to', outdir)

if __name__ == '__main__':
    main()
