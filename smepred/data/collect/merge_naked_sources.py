"""
merge_naked_sources.py — fold published siRNA datasets into our
naked-siRNA training set.

WHY
  Our LightGBM naked-siRNA model (the Rank tab) was trained on only 661 rows from
  HelixZero + siRNAmod. OligoFormer ships three additional published datasets in their
  data/ folder:
    Hu.csv   (Huesken et al. 2005)        — 2,361 rows, the gold-standard set
    Mix.csv  (Reynolds/Vickers/Ui-Tei/…)  —   472 rows
    Taka.csv (Takayuki 2007)              —   702 rows
  Total: +3,535 published rows, all with experimentally measured % inhibition.

FORMAT CONVERSION
  OligoFormer stores:
    siRNA      — 19-nt antisense (guide) strand
    label      — efficacy in [0, 1] (1 = fully silencing)
  SMEpred expects:
    sense      — 21-nt sense strand
    antisense  — 21-nt antisense strand
    efficacy   — % inhibition in [0, 100]

  To convert we:
    1. Pad each 19-nt antisense with the standard "UU" 3' overhang → 21 nt
    2. Derive the 21-nt sense as reverse-complement of the antisense
    3. Scale the label × 100 → percent
"""
from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).parent.parent
OLIGO_DATA = DATA_DIR / "oligoformer"
NORMAL_CSV = DATA_DIR / "normal_siRNA.csv"
OUT_PATH   = DATA_DIR / "normal_siRNA_extended.csv"


def _rev_comp(seq: str) -> str:
    t = str.maketrans("AUCG", "UAGC")
    return seq.upper().replace("T", "U").translate(t)[::-1]


def _convert(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Each row: 19-nt antisense + 0-1 label → 21-nt sense/antisense + 0-100 efficacy."""
    rows = []
    for _, row in df.iterrows():
        anti19 = str(row["siRNA"]).upper().replace("T", "U")
        if len(anti19) != 19 or any(c not in "AUGC" for c in anti19):
            continue
        antisense = anti19 + "UU"                # 21-nt with 3' overhang
        sense = _rev_comp(anti19) + "UU"         # 21-nt sense, standard duplex form
        efficacy = max(0.0, min(100.0, float(row["label"]) * 100.0))
        rows.append({
            "sense": sense, "antisense": antisense,
            "efficacy": round(efficacy, 2), "source": source,
        })
    return pd.DataFrame(rows)


def main():
    parts = []
    if NORMAL_CSV.exists():
        df = pd.read_csv(NORMAL_CSV)
        df["source"] = "smepred_existing"
        parts.append(df[["sense", "antisense", "efficacy", "source"]])
        print(f"  existing normal_siRNA.csv : {len(df)} rows")

    for name in ["Hu", "Mix", "Taka"]:
        p = OLIGO_DATA / f"{name}.csv"
        if not p.exists():
            print(f"  ! {p.name} missing — skipped")
            continue
        df = pd.read_csv(p)
        conv = _convert(df, source=f"oligoformer_{name}")
        parts.append(conv)
        print(f"  {p.name:>12}: {len(df):>5} raw -> {len(conv):>5} converted")

    combined = pd.concat(parts, ignore_index=True)
    before = len(combined)
    combined = combined.drop_duplicates(subset=["sense", "antisense"]).reset_index(drop=True)
    print(f"\n  combined        : {before} rows  ->  {len(combined)} after dedup")
    print(f"  efficacy stats  : min={combined.efficacy.min():.1f}  "
          f"mean={combined.efficacy.mean():.1f}  max={combined.efficacy.max():.1f}")
    print(f"  by source       : {combined['source'].value_counts().to_dict()}")

    combined.to_csv(OUT_PATH, index=False)
    print(f"\n  -> wrote {OUT_PATH}  ({len(combined)} rows, 4 columns)")


if __name__ == "__main__":
    main()
