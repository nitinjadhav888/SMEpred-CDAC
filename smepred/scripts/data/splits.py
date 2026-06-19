"""
splits.py — Paper-accurate train/validation split.

How the original SMEpred paper (Dar 2016) split its data:
  1. Arrange all sequences in DECREASING order of their efficacy (activity).
  2. Pick every 10th sequence starting from the 5th (5th, 15th, 25th, ...) for the
     INDEPENDENT VALIDATION set.
  3. The remaining sequences become the TRAINING/TESTING set.

Why this works:
  Sorting by efficacy first, then taking a regular stride, guarantees the validation set
  spans the entire efficacy range (high, medium, low) in the same proportion as the
  training set. This prevents a biased split where, say, all the high-efficacy sequences
  land in training and the model never learns to validate on them.

For Hetero-3031 this yields 303 validation + 2728 training.
For Homo-2110 this yields 210 validation + 1900 training.
The same rule scales to whatever number of clean rows we recover from the real data.
"""

import pandas as pd


def paper_split(df: pd.DataFrame, efficacy_col: str = "efficacy"):
    """
    Split a dataframe into (train, validation) using the paper's deterministic rule.

    Parameters
    ----------
    df           : dataframe with an efficacy column.
    efficacy_col : name of the efficacy column to sort by.

    Returns
    -------
    (train_df, val_df) : two dataframes. Indices are reset.
    """
    # Step 1: sort by descending efficacy (stable so ties keep input order).
    ordered = df.sort_values(efficacy_col, ascending=False, kind="mergesort").reset_index(drop=True)

    # Step 2: every 10th row starting at the 5th (0-based index 4, 14, 24, ...).
    val_mask = (ordered.index % 10) == 4
    val_df = ordered[val_mask].reset_index(drop=True)
    train_df = ordered[~val_mask].reset_index(drop=True)

    return train_df, val_df
