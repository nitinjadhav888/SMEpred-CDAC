"""
Validate SMEpred model predictions against independent CMsiRNAdb database.

Usage:
    cd D:\Helixx\smepred
    python tools\validate_on_cmsirnadb.py

Requires:
    - CMsiRNAdb TSV files in %TEMP%\cmsirnadb_PCSK9.tsv, _PNPLA3.tsv, etc.
    - SMEpred models (model_a.pkl, calibrator_cm.pkl, etc.) in models/
"""

import re
import sys
import warnings
from pathlib import Path
from typing import Dict, Optional, List, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
warnings.filterwarnings('ignore', message='X does not have valid feature names')

from src.features import extract_batch_gbm
from src.predictor import _get_model, _normalize_scores

# ─── mapping from IUPAC-style names to SMEpred single-letter codes ────────

# Ordered alias rules matching modification_codes.json logic
ALIAS_RULES = [
    ("unlocked nucleic acid",               "6"),
    ("vinyl-phosphonate-2",                 "3"),  # 5'-vinyl-P combined with 2'-OMe
    ("vinyl phosphonate",                   "3"),
    ("2-o-methyl",                          "M"),
    ("2-deoxy-2-fluoro",                    "F"),
    ("2-fluoro",                            "F"),
    ("2-deoxy",                             "D"),
    ("locked nucleic acid",                 "L"),
    ("methoxyethyl",                        "E"),
    ("glycol nucleic acid",                 "8"),
    ("threose nucleic acid",               "9"),
    ("arabino",                             "I"),
    ("4-thio",                              "N"),
    ("2-o-benzyl",                          "B"),
    ("pseudouridine",                       "W"),
    ("methylcytidine",                      "V"),
    ("5-methyl cytidine",                   "V"),
    ("inosine",                             "J"),
    ("2-thio uridine",                      "K"),
    ("dihydrouridine",                      "O"),
    ("abasic",                              "Q"),
    ("methylphosphonate",                   "R"),
    ("phosphoramidate",                     "H"),
    ("phosphorothioate",                    "S"),
    ("5-phosphate",                         "1"),
    ("3-phosphate",                         "2"),
]

BASE_LOOKUP = {
    "adenosine": "A",
    "uridine": "U",
    "cytidine": "C",
    "guanosine": "G",
    "thymidine": "T",
}

def iupac_to_symbol(iupac_name: str) -> Optional[str]:
    """Map a CMsiRNAdb IUPAC-style modification name to our single-letter code."""
    name = iupac_name.strip()
    if not name:
        return None

    # Plain base letter (already a symbol)
    if name in ("A", "U", "G", "C", "T"):
        return name

    # Extract the base for fallback
    name_lower = name.lower().replace("'", "").replace("-", " ").replace("_", " ")

    # Try alias rules in order
    for pattern, symbol in ALIAS_RULES:
        p = pattern.lower().replace("'", "").replace("-", " ").replace("_", " ")
        if p in name_lower:
            return symbol

    # Fallback: if it's just a plain nucleoside name (e.g., just "uridine"), return base
    for bname, bcode in BASE_LOOKUP.items():
        if bname in name_lower:
            return bcode

    return None


def parse_mod_string(mod_str: str) -> Dict[int, str]:
    """
    Parse a CMsiRNAdb Modification_Types column.
    Format: "1*2'-O-Methyluridine || 2*A || 3*G || ..."
    Returns: {1: 'M', 2: 'A', 3: 'G', ...}
    """
    result = {}
    if not mod_str or pd.isna(mod_str):
        return result
    parts = str(mod_str).split("||")
    for part in parts:
        part = part.strip()
        if not part or "*" not in part:
            continue
        pos_str, mod_name = part.split("*", 1)
        try:
            pos = int(pos_str.strip())
        except ValueError:
            continue
        # Some entries have trailing base letters like " || 21*T" which are already symbols
        mod_name_clean = mod_name.strip()
        code = iupac_to_symbol(mod_name_clean)
        if code is not None:
            result[pos] = code
    return result


def build_modified_string(base_seq: str, mod_map: Dict[int, str]) -> str:
    """
    Apply position-specific modifications to a base sequence.
    Only positions 1..len(base_seq) are considered.
    """
    if not base_seq or pd.isna(base_seq):
        return ""
    chars = []
    for i, base_char in enumerate(base_seq):
        pos = i + 1  # 1-indexed
        if pos in mod_map:
            # If the mapped symbol is a lowercase base letter, use it directly
            chars.append(mod_map[pos])
        else:
            chars.append(base_char)
    return "".join(chars)


def load_cmsirnadb(tsv_paths: List[str]) -> pd.DataFrame:
    """Load and concatenate CMsiRNAdb TSV files."""
    dfs = []
    for path in tsv_paths:
        p = Path(path)
        if p.exists():
            df = pd.read_csv(p, sep="\t", encoding="utf-8", low_memory=False)
            df["source_gene"] = p.stem.replace("patent_dataset_", "")
            dfs.append(df)
            print(f"  Loaded {p.stem}: {len(df)} entries")
        else:
            print(f"  File not found: {p}")
    if not dfs:
        raise FileNotFoundError("No CMsiRNAdb TSV files found")
    return pd.concat(dfs, ignore_index=True)


def prepare_validation_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert CMsiRNAdb rows to SMEpred-compatible format.
    Returns DataFrame with columns:
        sense_mod, antisense_mod, base_sense, base_antisense, inhibition
    """
    rows = []
    total = len(df)
    skipped = 0

    for idx, row in df.iterrows():
        if idx % 500 == 0:
            print(f"  Processing row {idx}/{total}")

        base_sense = str(row.get("Sense_seqence", ""))
        base_anti = str(row.get("Antisense_seqence", ""))
        mod_types_sense = row.get("Modification_Types_Sense_strand", "")
        mod_types_anti = row.get("Modification_Types_Antisense_strand", "")

        # Skip entries with no modification info and no sequence
        if (pd.isna(mod_types_sense) or not str(mod_types_sense).strip()) and \
           (pd.isna(mod_types_anti) or not str(mod_types_anti).strip()):
            skipped += 1
            continue

        inhibition = row.get("Inhibition", None)
        if pd.isna(inhibition) or inhibition == "":
            skipped += 1
            continue
        inhibition = float(inhibition)

        # Skip if sequences are missing or too short
        if pd.isna(base_sense) or pd.isna(base_anti) or len(base_sense) < 19 or len(base_anti) < 19:
            skipped += 1
            continue

        # Parse modification types
        sense_mods = parse_mod_string(str(mod_types_sense))
        anti_mods = parse_mod_string(str(mod_types_anti))

        # Build modified strings
        sense_mod = build_modified_string(base_sense.strip(), sense_mods)
        anti_mod = build_modified_string(base_anti.strip(), anti_mods)

        # Validate lengths
        if len(sense_mod) != len(base_sense.strip()) or len(anti_mod) != len(base_anti.strip()):
            skipped += 1
            continue

        # Filter: only 0-100 inhibition range (our model's training range)
        if inhibition < 0 or inhibition > 100:
            skipped += 1
            continue

        # Uppercase to ensure canonical symbols
        sense_mod = sense_mod.upper()
        anti_mod = anti_mod.upper()

        rows.append({
            "sense_mod": sense_mod,
            "antisense_mod": anti_mod,
            "base_sense": base_sense.strip().upper(),
            "base_antisense": base_anti.strip().upper(),
            "inhibition": inhibition,
        })

    result = pd.DataFrame(rows)
    print(f"\n  Prepared {len(result)} entries for validation ({skipped} skipped)")
    return result


def run_validation(data: pd.DataFrame, batch_size: int = 512) -> pd.DataFrame:
    """
    Run SMEpred model predictions on prepared validation data.
    """
    model = _get_model("A")
    n = len(data)
    all_preds = []

    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch = data.iloc[start:end]

        s_list = batch["sense_mod"].tolist()
        a_list = batch["antisense_mod"].tolist()
        bs_list = batch["base_sense"].tolist()
        ba_list = batch["base_antisense"].tolist()

        X = extract_batch_gbm(s_list, a_list,
                              base_sense_list=bs_list,
                              base_antisense_list=ba_list)
        raw = model.predict(X)
        cal = _normalize_scores(raw, calibrator_key="cm")
        all_preds.extend(cal.tolist())

        if (start // batch_size) % 5 == 0:
            print(f"  Predicted batch {start//batch_size + 1}/{(n-1)//batch_size + 1}")

    data["predicted_calibrated"] = np.array(all_preds)
    return data


def print_metrics(data: pd.DataFrame):
    """Compute and print validation metrics."""
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import mean_absolute_error

    y_true = data["inhibition"].values
    y_pred = data["predicted_calibrated"].values

    # Filter valid
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_t, y_p = y_true[mask], y_pred[mask]

    pcc, _ = pearsonr(y_t, y_p)
    spr, _ = spearmanr(y_t, y_p)
    mae = mean_absolute_error(y_t, y_p)

    print(f"\n{'='*50}")
    print(f"VALIDATION RESULTS")
    print(f"{'='*50}")
    print(f"  N entries tested : {len(y_t)}")
    print(f"  PCC             : {pcc:.4f}")
    print(f"  Spearman ρ      : {spr:.4f}")
    print(f"  MAE             : {mae:.2f} pts")
    print(f"  RMSE            : {np.sqrt(np.mean((y_t - y_p)**2)):.2f} pts")
    print(f"  Mean true       : {y_t.mean():.2f}")
    print(f"  Mean predicted  : {y_p.mean():.2f}")
    print(f"  Std true        : {y_t.std():.2f}")
    print(f"  Std predicted   : {y_p.std():.2f}")
    print(f"{'='*50}")

    # Per-decile buckets
    print(f"\n  Performance by true-value decile:")
    for i in range(10):
        lo, hi = i * 10, (i + 1) * 10
        mask_d = (y_t >= lo) & (y_t < hi)
        if mask_d.sum() > 5:
            d_mae = mean_absolute_error(y_t[mask_d], y_p[mask_d])
            d_bias = (y_p[mask_d] - y_t[mask_d]).mean()
            print(f"    [{lo:3d}-{hi:3d})  n={mask_d.sum():5d}  MAE={d_mae:5.2f}  bias={d_bias:+6.2f}")

    return {"pcc": pcc, "spearman": spr, "mae": mae}


if __name__ == "__main__":
    import os

    tmp = os.environ.get("TEMP", r"C:\Users\Nilesh\AppData\Local\Temp")
    tsv_files = [
        os.path.join(tmp, "cmsirnadb_PCSK9.tsv"),
        os.path.join(tmp, "cmsirnadb_PNPLA3.tsv"),
    ]

    print("Step 1: Loading CMsiRNAdb data...")
    df = load_cmsirnadb(tsv_files)
    print(f"  Total entries loaded: {len(df)}")

    print("\nStep 2: Preparing validation data...")
    vdata = prepare_validation_data(df)
    print(f"  Validation entries: {len(vdata)}")

    if len(vdata) < 10:
        print("Not enough validation data. Exiting.")
        sys.exit(1)

    print("\nStep 3: Running predictions...")
    vdata = run_validation(vdata)

    print("\nStep 4: Computing metrics...")
    metrics = print_metrics(vdata)

    # Save results
    out_path = Path(__file__).parent.parent / "logs" / "cmsirnadb_validation.csv"
    out_path.parent.mkdir(exist_ok=True)
    vdata.to_csv(out_path, index=False)
    print(f"\nResults saved to {out_path}")
