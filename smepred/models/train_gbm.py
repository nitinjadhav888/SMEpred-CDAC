"""
train_gbm.py — LightGBM training for SMEpred (accuracy rebuild).

WHY THIS REPLACES THE RBF-SVR
  The paper-faithful SVR + MNC pipeline plateaued at ~0.37 PCC on the real HelixZero
  catalog because (a) RBF-SVR scales ~O(n^2) so we could only use a 6k subsample, and
  (b) MNC composition ignores WHERE modifications sit and ignores the assay condition,
  which is the dominant confound in patent data (same sequence, different dose → very
  different inhibition).

WHAT THIS DOES DIFFERENTLY
  • LightGBM (gradient-boosted trees): trains on ALL rows in seconds, captures
    non-linear position × modification × condition interactions.
  • Richer features (src/features.extract_batch_gbm): base+modified MNC, position-aware
    modification density, GC content, and the experimental condition (dose, time).
  • Honest evaluation via a GENE-GROUPED split: whole target genes are held out, so the
    reported accuracy reflects performance on NEW targets, not memorized patent motifs.
  • Predicts efficacy DIRECTLY on the 0–100 inhibition scale (no min-max rescaling).

OUTPUTS (overwrites the SVR pkls so the API/CLI pick them up transparently)
  model_a.pkl / model_b.pkl / model_c.pkl  — the cm-siRNA efficacy GBM (A/B/C share it;
        the GBM unifies what used to be three separate feature recipes)
  model_normal.pkl                         — naked (unmodified) siRNA efficacy GBM
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import joblib
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import GroupShuffleSplit
from lightgbm import LGBMRegressor, early_stopping

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.features import extract_batch_gbm

DATA_DIR   = Path(__file__).parent.parent / "data"
MODELS_DIR = Path(__file__).parent

HETERO_TRAIN = DATA_DIR / "hetero_train_2728_mrna.csv"
HETERO_VAL   = DATA_DIR / "hetero_val_303_mrna.csv"
_HAS_MRNA = HETERO_TRAIN.exists() and HETERO_VAL.exists()
if not _HAS_MRNA:
    HETERO_TRAIN = DATA_DIR / "hetero_train_2728.csv"
    HETERO_VAL   = DATA_DIR / "hetero_val_303.csv"

MRNA_MEAN_FILE = DATA_DIR / "mrna_feature_means.json"
# Prefer the multi-source naked dataset when present (4k rows vs 661)
_EXT = DATA_DIR / "normal_siRNA_extended.csv"
NORMAL_CSV   = _EXT if _EXT.exists() else (DATA_DIR / "normal_siRNA.csv")

LGBM_PARAMS = dict(
    n_estimators=800,
    learning_rate=0.03,
    num_leaves=63,
    max_depth=-1,
    min_child_samples=30,
    subsample=0.8,
    subsample_freq=1,
    colsample_bytree=0.7,
    reg_lambda=1.0,
    random_state=42,
    n_jobs=-1,
    verbose=-1,
)


def _metrics(y_true, y_pred):
    """Return (PCC, Spearman, MAE, RMSE). MAE/RMSE are in % inhibition points."""
    pcc = pearsonr(y_true, y_pred)[0]
    sp = spearmanr(y_true, y_pred)[0]
    err = y_pred - y_true
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err ** 2)))
    return pcc, sp, mae, rmse


_EXCLUDED_COLS = {'sense', 'antisense', 'base_sense', 'base_antisense', 'efficacy',
                  'concentration_nM', 'time_h', 'target_gene', 'source'}
_MRNA_FEAT_NAMES: list[str] = []


def _detect_mrna_cols(df: pd.DataFrame) -> list[str]:
    """Auto-detect extra feature columns (all columns not in the base schema)."""
    return [c for c in df.columns if c not in _EXCLUDED_COLS]


def _mrna_features(df: pd.DataFrame) -> list[np.ndarray]:
    """Extract extra feature arrays per row; nan for missing."""
    cols = _detect_mrna_cols(df)
    if not cols:
        return [None] * len(df)
    arr = df[cols].to_numpy(dtype=np.float32)
    return [arr[i] for i in range(len(arr))]


def _impute_mrna(train_df: pd.DataFrame):
    """Compute mean over training set and save for inference."""
    cols = _detect_mrna_cols(train_df)
    arr = train_df[cols].to_numpy(dtype=np.float32)
    means = np.nanmean(arr, axis=0)
    for i, col in enumerate(cols):
        train_df[col] = train_df[col].fillna(means[i])
    import json
    MRNA_MEAN_FILE.write_text(
        json.dumps({col: float(means[i]) for i, col in enumerate(cols)}, indent=2))
    global _MRNA_FEAT_NAMES
    _MRNA_FEAT_NAMES = cols


def _features(df: pd.DataFrame) -> np.ndarray:
    bs = list(df["base_sense"]) if "base_sense" in df.columns else None
    ba = list(df["base_antisense"]) if "base_antisense" in df.columns else None
    cc = list(df["concentration_nM"]) if "concentration_nM" in df.columns else None
    tt = list(df["time_h"]) if "time_h" in df.columns else None
    mf = _mrna_features(df)
    return extract_batch_gbm(list(df["sense"]), list(df["antisense"]),
                             base_sense_list=bs, base_antisense_list=ba,
                             conc_list=cc, time_list=tt, mrna_feat_list=mf)


def _n_feat_extra() -> int:
    """Number of extra features beyond the 152-d base."""
    try:
        import json
        raw = json.loads(MRNA_MEAN_FILE.read_text())
        return len(raw)
    except Exception:
        return 0


def _gene_grouped_split(df: pd.DataFrame, test_frac=0.18, seed=42):
    """Hold out whole target genes for validation (no gene appears in both sets)."""
    if "target_gene" not in df.columns or df["target_gene"].nunique() < 3:
        # fallback: random split
        val = df.sample(frac=test_frac, random_state=seed)
        return df.drop(val.index).reset_index(drop=True), val.reset_index(drop=True)
    gss = GroupShuffleSplit(n_splits=1, test_size=test_frac, random_state=seed)
    tr_idx, va_idx = next(gss.split(df, groups=df["target_gene"]))
    return df.iloc[tr_idx].reset_index(drop=True), df.iloc[va_idx].reset_index(drop=True)


def train_cm_model():
    print("\n=== cm-siRNA (modified) GBM ===")
    df = pd.concat([pd.read_csv(HETERO_TRAIN), pd.read_csv(HETERO_VAL)], ignore_index=True)
    extra_cols = _detect_mrna_cols(df)
    has_extra = len(extra_cols) > 0
    print(f"  rows: {len(df)}  genes: {df['target_gene'].nunique()}  extra_feats: {len(extra_cols)}")

    if has_extra:
        _impute_mrna(df)
        print(f"  Extra features -> dim {152 + len(extra_cols)} (152 + {len(extra_cols)})")

    X_all, y_all = _features(df), df["efficacy"].to_numpy(float)

    # ── Eval 1: gene-grouped (honest "new target" accuracy → Rank-tab use case) ──
    train_df, val_df = _gene_grouped_split(df)
    held = sorted(set(val_df["target_gene"]) - set(train_df["target_gene"]))
    Xg_tr, yg_tr = _features(train_df), train_df["efficacy"].to_numpy(float)
    Xg_va, yg_va = _features(val_df),   val_df["efficacy"].to_numpy(float)
    mg = LGBMRegressor(**LGBM_PARAMS)
    mg.fit(Xg_tr, yg_tr, eval_set=[(Xg_va, yg_va)],
           callbacks=[early_stopping(50, verbose=False)])
    pcc_g, sp_g, mae_g, rmse_g = _metrics(yg_va, mg.predict(Xg_va))
    print(f"  train {len(train_df)} / val {len(val_df)}  | held-out genes: {held}")
    print(f"  [new-target / gene-grouped]   PCC={pcc_g:.4f}  Spearman={sp_g:.4f}  MAE={mae_g:.2f}  RMSE={rmse_g:.2f}")

    # ── Eval 2: random split (within-target → Single/Multi-Mod ranking use case) ──
    rng = np.random.default_rng(42)
    idx = rng.permutation(len(df)); cut = int(len(df) * 0.82)
    tr, va = idx[:cut], idx[cut:]
    mr = LGBMRegressor(**LGBM_PARAMS)
    mr.fit(X_all[tr], y_all[tr], eval_set=[(X_all[va], y_all[va])],
           callbacks=[early_stopping(50, verbose=False)])
    pcc_r, sp_r, mae_r, rmse_r = _metrics(y_all[va], mr.predict(X_all[va]))
    print(f"  [modification-ranking / random] PCC={pcc_r:.4f}  Spearman={sp_r:.4f}  MAE={mae_r:.2f}  RMSE={rmse_r:.2f}  (SVR baseline 0.37)")

    # ── Final model: refit on ALL data for deployment ──
    model = LGBMRegressor(**LGBM_PARAMS)
    model.fit(X_all, y_all)
    for name in ("model_a", "model_b", "model_c"):
        joblib.dump(model, MODELS_DIR / f"{name}.pkl")
    print(f"  Saved (refit on all {len(df)} rows) -> model_a/b/c.pkl")
    return model


def train_normal_model():
    print("\n=== naked (unmodified) siRNA GBM ===")
    df = pd.read_csv(NORMAL_CSV)
    print(f"  rows: {len(df)}")

    # Source-aware feature stacking: each dataset (Huesken / Mix / Taka / our existing)
    # has a different experimental distribution (cell line, assay). Without telling the
    # model which source a row came from, it fits the population mean and loses signal.
    # We append a small one-hot of the source so the model can learn per-source offsets.
    if "source" in df.columns:
        sources = sorted(df["source"].unique())
        src_idx = {s: i for i, s in enumerate(sources)}
        src_onehot = np.zeros((len(df), len(sources)), dtype=np.float32)
        for i, s in enumerate(df["source"]):
            src_onehot[i, src_idx[s]] = 1.0
        print(f"  sources: {sources}")
    else:
        src_onehot = np.zeros((len(df), 0), dtype=np.float32)

    val_idx = df.sample(frac=0.18, random_state=42).index
    train_mask = ~df.index.isin(val_idx)

    # sequence features
    X_seq_all = extract_batch_gbm(list(df["sense"]), list(df["antisense"]))
    X_all = np.concatenate([X_seq_all, src_onehot], axis=1)
    y_all = df["efficacy"].to_numpy(float)

    X_tr, X_va = X_all[train_mask], X_all[~train_mask]
    y_tr, y_va = y_all[train_mask], y_all[~train_mask]

    # With 4k+ rows we can support more capacity than the 661-row model
    model = LGBMRegressor(**{**LGBM_PARAMS, "n_estimators": 1000, "num_leaves": 63,
                             "learning_rate": 0.02, "min_child_samples": 20})
    model.fit(X_tr, y_tr, eval_set=[(X_va, y_va)],
              callbacks=[early_stopping(80, verbose=False)])
    pcc, sp, mae, rmse = _metrics(y_va, model.predict(X_va))
    print(f"  Held-out (all-source)  PCC={pcc:.4f}  Spearman={sp:.4f}  MAE={mae:.2f}  RMSE={rmse:.2f}")

    # Per-source PCC so we can tell whether the model ranks well WITHIN each dataset
    if "source" in df.columns:
        pred_va = model.predict(X_va)
        src_va = df["source"].to_numpy()[~train_mask]
        for s in sources:
            mask = src_va == s
            if mask.sum() >= 30:
                p, _, _, _ = _metrics(y_va[mask], pred_va[mask])
                print(f"     within {s:<22} (n={int(mask.sum()):>4}):  PCC={p:.4f}")

    # Refit on all rows for deployment
    model_final = LGBMRegressor(**{**LGBM_PARAMS, "n_estimators": model.best_iteration_ or 800,
                                   "num_leaves": 63, "learning_rate": 0.02, "min_child_samples": 20})
    model_final.fit(X_all, y_all)
    joblib.dump({"model": model_final, "sources": list(src_idx.keys()) if "source" in df.columns else []},
                MODELS_DIR / "model_normal.pkl")
    print(f"  Saved (refit on all {len(df)} rows) -> model_normal.pkl")


def main():
    import sys as _sys
    only_normal = "--only-normal" in _sys.argv
    only_cm     = "--only-cm" in _sys.argv
    print("=" * 60)
    print("HelixZero-CMS LightGBM Training")
    print(f"naked dataset: {NORMAL_CSV.name}")
    print("=" * 60)
    if not only_normal:
        train_cm_model()
    if not only_cm:
        train_normal_model()
    print("\nDone. Models saved to models/")


if __name__ == "__main__":
    main()
