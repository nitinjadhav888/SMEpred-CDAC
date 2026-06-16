"""
parse_sirnamod.py — AUXILIARY parser for the siRNAmod export (File 1).

FILE: C:\\Helixx\\Data-(1)-csv.csv.xls  (plain CSV despite the extension)
COLUMNS: siRNAmodDB, PMID, sense seq, sense modification name(s),
         antisense seq, antisense modification name(s), inhibition %.

LIMITATION
  This export lists modification NAMES per strand but gives NO positions. Without
  positions we cannot reconstruct the 35-symbol positioned strands that Models A/B/C
  consume, so modified rows here are NOT directly compatible with the File-2 feature
  space. Its clean, valuable contribution is the UNMODIFIED rows (no modification on
  either strand), which only need base sequences -> we use them to augment the normal
  (unmodified) siRNA training set produced from File 2.

WHAT THIS SCRIPT DOES
  1. Parse and clean sequences + efficacy (clip to [0,100]).
  2. Select rows whose sense AND antisense modification fields indicate NO modification
     (value "0", empty, or a plain canonical base) -> genuine unmodified siRNAs.
  3. Merge those into data/normal_siRNA.csv (dedup), growing the normal ranker dataset.
  Modified rows are counted and reported but skipped (positionless).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.collect.clean_utils import clean_sequence, parse_efficacy, valid_length, CANONICAL

SIRNAMOD = Path(r"C:\Helixx\Data-(1)-csv.csv.xls")
DATA_DIR = Path(__file__).parent.parent
NORMAL_CSV = DATA_DIR / "normal_siRNA.csv"

_NO_MOD_TOKENS = {"", "0", "nan", "none"}


def _is_unmodified(mod_field) -> bool:
    """True if a modification field denotes no modification."""
    s = str(mod_field).strip()
    if s.lower() in _NO_MOD_TOKENS:
        return True
    # a single canonical base symbol also means effectively unmodified
    return s.upper() in CANONICAL


def main():
    if not SIRNAMOD.exists():
        sys.exit(f"ERROR: siRNAmod export not found at {SIRNAMOD}")

    df = pd.read_csv(SIRNAMOD)
    print(f"Parsing {SIRNAMOD.name} ... raw rows: {len(df)}")

    sense_col = "Sequence of sense strand"
    anti_col = "Sequence of antisense strand"
    sense_mod = "All Modification (sense strand)"
    anti_mod = "All Modification name (antisense strand)"
    eff_col = "Biological inhibition percentage"

    rows = []
    n_modified = n_badseq = n_badeff = 0
    for _, r in df.iterrows():
        eff = parse_efficacy(r[eff_col])
        if eff is None:
            n_badeff += 1
            continue
        s = clean_sequence(r[sense_col])
        a = clean_sequence(r[anti_col])
        if not (valid_length(s) and valid_length(a)):
            n_badseq += 1
            continue
        if not (_is_unmodified(r[sense_mod]) and _is_unmodified(r[anti_mod])):
            n_modified += 1
            continue
        rows.append({"sense": s, "antisense": a, "efficacy": eff})

    aux = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    print(f"  unmodified rows usable : {len(aux)}")
    print(f"  skipped (modified)     : {n_modified}")
    print(f"  skipped (bad sequence) : {n_badseq}")
    print(f"  skipped (bad efficacy) : {n_badeff}")

    # merge into the normal set from File 2
    if NORMAL_CSV.exists():
        existing = pd.read_csv(NORMAL_CSV)
        before = len(existing)
        combined = pd.concat([existing, aux], ignore_index=True)
        combined = combined.drop_duplicates(subset=["sense", "antisense", "efficacy"]).reset_index(drop=True)
        combined.to_csv(NORMAL_CSV, index=False)
        print(f"\n  normal_siRNA.csv: {before} -> {len(combined)} "
              f"(+{len(combined) - before} from siRNAmod)")
    else:
        aux.to_csv(NORMAL_CSV, index=False)
        print(f"\n  normal_siRNA.csv created with {len(aux)} rows")

    print("Done.")


if __name__ == "__main__":
    main()
