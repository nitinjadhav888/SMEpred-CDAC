"""
augment_mrna_features.py — adds mRNA target-site features to training data.

For every row in hetero_train_2728.csv / hetero_val_303.csv where the
target site can be located in the gene's mRNA transcript, computes:
  - target_position_norm: position in transcript (0–1)
  - target_site_gc: GC% of 21-nt target site
  - upstream_20_gc: GC% of 20 nt upstream
  - downstream_20_gc: GC% of 20 nt downstream
  - target_mfe: ViennaRNA MFE of the local region

Outputs augmented CSVs with 5 extra columns.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import pandas as pd
from src.mrna_features import (
    load_mrna_cache, find_target_position, _compute, feature_count,
)

DATA_DIR = Path(__file__).parent.parent
TRAIN_CSV = DATA_DIR / "hetero_train_2728.csv"
VAL_CSV = DATA_DIR / "hetero_val_303.csv"
OUT_TRAIN = DATA_DIR / "hetero_train_2728_mrna.csv"
OUT_VAL = DATA_DIR / "hetero_val_303_mrna.csv"


def augment(df: pd.DataFrame, cache: dict) -> pd.DataFrame:
    """Add mRNA feature columns to a dataframe."""
    n_feat = feature_count()
    feat_names = ["target_position_norm", "target_site_gc",
                  "upstream_20_gc", "downstream_20_gc", "target_mfe"]

    feat_cols = []
    for _, row in df.iterrows():
        gene = row["target_gene"]
        as_seq = row["base_antisense"]
        entry = cache.get(gene)
        if entry is None:
            feat_cols.append([np.nan] * n_feat)
            continue
        pos = find_target_position(as_seq, entry["sequence"])
        if pos is None:
            feat_cols.append([np.nan] * n_feat)
            continue
        feat = _compute(as_seq, entry["sequence"], pos)
        feat_cols.append(feat.tolist())

    feat_df = pd.DataFrame(feat_cols, columns=feat_names, index=df.index)
    return pd.concat([df, feat_df], axis=1)


def main():
    cache = load_mrna_cache()
    print(f"Loaded {len(cache)} mRNA sequences")

    for name, csv_in, csv_out in [("train", TRAIN_CSV, OUT_TRAIN),
                                    ("val", VAL_CSV, OUT_VAL)]:
        df = pd.read_csv(csv_in)
        print(f"\n{name}: {len(df)} rows, {df['target_gene'].nunique()} genes")
        aug = augment(df, cache)
        aug.to_csv(csv_out, index=False)

        # Report coverage
        ok = aug["target_position_norm"].notna().sum()
        pct = 100 * ok / len(aug)
        print(f"  {ok}/{len(aug)} aligned ({pct:.1f}%)")
        if ok > 0:
            print(f"  Features: {aug[['target_position_norm','target_site_gc','upstream_20_gc','downstream_20_gc','target_mfe']].describe().round(3)}")

    print("\nDone. Augmented files saved.")
    print(f"  {OUT_TRAIN}")
    print(f"  {OUT_VAL}")


if __name__ == "__main__":
    main()
