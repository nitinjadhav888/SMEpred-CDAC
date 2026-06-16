"""
parse_helix_catalog.py — PRIMARY parser for the HelixZero 43k catalog.

WHY A CUSTOM PARSER
  The file `HelixZero_Biological_Catalog_43k.csv.xls` is a CSV in name only. Several
  columns (`Modification_locations_*`, `position_*`) contain UNQUOTED comma-lists like
  `1,2,3,4,5,...`, so the number of comma-separated fields varies wildly per row
  (only ~3500 of 43k rows have the nominal 26 fields). pandas.read_csv cannot parse it.

  Fortunately the per-position modification data is fully recoverable from one robust
  signal: each strand is annotated as a stream of `position*name` tokens separated by
  ` || `, e.g.  `1*2'-O-Methylcytidine || 2*2'-O-Methyladenosine || ... || 25*...`.
  We verified across 857k tokens that NO modification name contains a comma, so a simple
  token regex extracts every (position, name) pair cleanly.

HOW WE RECONSTRUCT EACH siRNA
  1. Efficacy: the row ID ends in the inhibition value (e.g. `...-48h-88.00` -> 88.0),
     which we verified matches the `Inhibition` column. We read it from the ID.
  2. Token stream: extract all (position, name) pairs in order.
  3. Strand split: positions run 1..N for the ANTISENSE strand, then reset to 1..M for
     the SENSE strand (antisense block appears first in the file). We detect the reset
     (position not strictly increasing) to split the stream into the two strands.
  4. Per token we derive TWO things:
       - the canonical BASE  (from the name's base suffix: "...uridine" -> U, etc.)
       - the 35-alphabet SYMBOL (via the alias map; canonical positions keep their base)
     This yields both the unmodified base strand and the modified-symbol strand without
     needing the malformed sequence columns at all.

OUTPUTS (written to data/)
  - hetero_train_2728.csv / hetero_val_303.csv : all cm-siRNAs (paper split)
  - normal_siRNA.csv                           : rows with NO modification (both strands
                                                 fully canonical) for the normal ranker
  - homo_train.csv / homo_val.csv              : chemically homogeneous subset
                                                 (a single modification class used)
  - data/collect/unmapped_report.txt           : modification names that matched no alias
"""

import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.collect.clean_utils import (
    parse_efficacy, map_modification, valid_length, write_unmapped_report,
)
from data.collect.splits import paper_split

# ─── paths ────────────────────────────────────────────────────────────────────

CATALOG = Path(r"C:\Helixx\HelixZero_Biological_Catalog_43k.csv.xls")
DATA_DIR = Path(__file__).parent.parent
COLLECT_DIR = Path(__file__).parent

# ─── regexes ──────────────────────────────────────────────────────────────────

# position*name, name ends at ' || ' or a comma (column boundary). Names have no commas.
_TOKEN_RE = re.compile(r"(\d+)\*([^|,]+?)(?=\s*\|\||,)")
# trailing inhibition value in the ID. NOTE: negatives appear as a double dash, e.g.
# '...-24h--8.87' (separator dash + negative number), so the value itself may start
# with '-'. Capturing the sign is essential — otherwise '-8.87' reads as '+8.87'.
_ID_EFF_RE = re.compile(r"-(-?\d+(?:\.\d+)?)\s*$")
# experimental condition embedded in the ID: '...-<conc>n-<time>h-...'
# e.g. '...-100n-48h-88.00' -> concentration '100n', time '48h'.
_ID_COND_RE = re.compile(r"-([0-9.]+n)-([0-9.]+h)-")

# base-suffix -> canonical RNA base
_BASE_SUFFIX = [
    ("adenosine", "A"), ("adenine", "A"),
    ("uridine", "U"), ("uracil", "U"),
    ("thymidine", "U"), ("thymine", "U"),   # DNA T -> RNA U
    ("cytidine", "C"), ("cytosine", "C"),
    ("guanosine", "G"), ("guanine", "G"),
]


def _base_from_name(name: str):
    """Derive the underlying canonical RNA base from a position name."""
    s = name.strip()
    up = s.upper()
    if up in ("A", "U", "G", "C", "T"):
        return "U" if up == "T" else up
    low = s.lower()
    for suffix, base in _BASE_SUFFIX:
        if suffix in low:
            return base
    return None   # e.g. abasic / unknown


def _parse_line(line: str, unmapped: Counter):
    """
    Parse one raw line into a record dict, or return None if unusable.
    """
    # efficacy from the ID (field before the first comma)
    id_field = line.split(",", 1)[0]
    m = _ID_EFF_RE.search(id_field)
    if not m:
        return None
    efficacy = parse_efficacy(m.group(1))   # negatives are clipped to 0 here
    if efficacy is None:
        return None

    # experimental condition (concentration, time) for the dose-controlled subset
    cm = _ID_COND_RE.search(id_field)
    condition = f"{cm.group(1)}_{cm.group(2)}" if cm else "unknown"
    # numeric condition values (model features): '100n'->100.0 nM, '48h'->48.0 h.
    # NaN when not parseable so the trainer can fill with a reference condition.
    if cm:
        concentration_nM = float(cm.group(1).rstrip("n"))
        time_h = float(cm.group(2).rstrip("h"))
    else:
        concentration_nM = float("nan")
        time_h = float("nan")

    # target gene: column index 4. The leading columns (ID, patent_ID, auth_status,
    # accession, Target_Gene) are clean simple values, so a plain comma split is safe.
    parts = line.split(",")
    target_gene = parts[4].strip() if len(parts) > 4 else ""

    # token stream
    tokens = [(int(p), n.strip()) for p, n in _TOKEN_RE.findall(line)]
    if not tokens:
        return None

    # split into strands at each position reset (pos not strictly increasing)
    runs = []
    cur = []
    prev = 0
    for pos, name in tokens:
        if pos <= prev and cur:
            runs.append(cur)
            cur = []
        cur.append((pos, name))
        prev = pos
    if cur:
        runs.append(cur)

    if len(runs) < 2:
        return None
    antisense_run, sense_run = runs[0], runs[1]   # antisense block is first in the file

    def build(run):
        # positions must be contiguous 1..N
        positions = [p for p, _ in run]
        n = len(run)
        if sorted(positions) != list(range(1, n + 1)):
            return None, None, False
        base = [""] * n
        symb = [""] * n
        modified = False
        for pos, name in run:
            b = _base_from_name(name)
            sym = map_modification(name)
            if sym is None:                # unmapped modification
                unmapped[name] += 1
                sym = ""                   # treat as canonical for this position
            if b is None:
                b = "A"                    # placeholder for abasic/unknown (rare)
            base[pos - 1] = b
            if sym == "":
                symb[pos - 1] = b          # canonical -> keep base
            else:
                symb[pos - 1] = sym
                modified = True
        return "".join(base), "".join(symb), modified

    base_as, mod_as, mod_as_flag = build(antisense_run)
    base_ss, mod_ss, mod_ss_flag = build(sense_run)
    if base_as is None or base_ss is None:
        return None

    # length filter on the canonical base strands
    if not (valid_length(base_as) and valid_length(base_ss)):
        return None

    is_modified = mod_as_flag or mod_ss_flag
    # distinct non-canonical symbols used (for homogeneous-subset detection)
    distinct_mods = {c for c in (mod_as + mod_ss) if c not in ("A", "U", "G", "C", "T")}

    return {
        "sense": mod_ss,             # modified-symbol sense strand (model input)
        "antisense": mod_as,         # modified-symbol antisense strand (model input)
        "base_sense": base_ss,
        "base_antisense": base_as,
        "efficacy": efficacy,
        "is_modified": is_modified,
        "n_distinct_mods": len(distinct_mods),
        "condition": condition,
        "concentration_nM": concentration_nM,
        "time_h": time_h,
        "target_gene": target_gene,
    }


def main():
    if not CATALOG.exists():
        sys.exit(f"ERROR: catalog not found at {CATALOG}")

    print(f"Parsing {CATALOG.name} ...")
    rows = []
    unmapped = Counter()
    total = kept = 0
    with CATALOG.open(encoding="utf-8") as f:
        next(f)  # skip header
        for line in f:
            total += 1
            rec = _parse_line(line, unmapped)
            if rec is not None:
                rows.append(rec)
                kept += 1

    df = pd.DataFrame(rows)
    print(f"  raw rows           : {total}")
    print(f"  parsed/kept        : {kept}")
    print(f"  dropped (unusable) : {total - kept}")

    # dedup on the modelled fields
    before = len(df)
    df = df.drop_duplicates(subset=["sense", "antisense", "efficacy"]).reset_index(drop=True)
    print(f"  duplicates removed : {before - len(df)}")
    print(f"  unique rows        : {len(df)}")
    print(f"  modified rows      : {int(df['is_modified'].sum())}")
    print(f"  unmodified rows    : {int((~df['is_modified']).sum())}")
    print(f"  efficacy: min={df['efficacy'].min():.1f} max={df['efficacy'].max():.1f} "
          f"mean={df['efficacy'].mean():.1f}  zeros(non-functional)={int((df['efficacy']==0).sum())}")
    print(f"  top conditions     : "
          + ", ".join(f"{c}={n}" for c, n in df['condition'].value_counts().head(5).items()))

    # unmapped report
    write_unmapped_report(unmapped, COLLECT_DIR / "unmapped_report.txt",
                          title="HelixZero catalog - unmapped modification names")
    print(f"  unmapped mod names : {len(unmapped)} distinct "
          f"(report -> {COLLECT_DIR / 'unmapped_report.txt'})")

    # ── Hetero (all cm-siRNAs) ──
    # Include base (unmodified) strands: combined base+modified composition predicts
    # markedly better, because near-fully-modified strands otherwise collapse to M/F
    # and lose all underlying sequence signal.
    hetero = df[["sense", "antisense", "base_sense", "base_antisense", "efficacy",
                 "concentration_nM", "time_h", "target_gene"]]
    h_train, h_val = paper_split(hetero)
    h_train.to_csv(DATA_DIR / "hetero_train_2728.csv", index=False)
    h_val.to_csv(DATA_DIR / "hetero_val_303.csv", index=False)
    print(f"\n  Hetero -> train {len(h_train)} / val {len(h_val)}")

    # ── Normal (unmodified) ──
    normal = df[~df["is_modified"]][["base_sense", "base_antisense", "efficacy"]].copy()
    normal.columns = ["sense", "antisense", "efficacy"]
    normal = normal.drop_duplicates().reset_index(drop=True)
    normal.to_csv(DATA_DIR / "normal_siRNA.csv", index=False)
    print(f"  Normal (unmodified) -> {len(normal)} rows")

    # ── Homo (dose-controlled: the single most common experimental condition) ──
    # The paper's "homogeneous" set held experimental conditions constant. Inhibition
    # depends heavily on dose/time, so fixing (concentration, time) removes that
    # confound and gives a cleaner efficacy signal than the mixed-condition Hetero set.
    known = df[df["condition"] != "unknown"]
    if len(known):
        top_cond = known["condition"].value_counts().index[0]
        homo = df[df["condition"] == top_cond][["sense", "antisense", "base_sense", "base_antisense", "efficacy"]]
        if len(homo) >= 20:
            homo_train, homo_val = paper_split(homo)
            homo_train.to_csv(DATA_DIR / "homo_train.csv", index=False)
            homo_val.to_csv(DATA_DIR / "homo_val.csv", index=False)
            print(f"  Homo (condition {top_cond}) -> train {len(homo_train)} / val {len(homo_val)}")
        else:
            print(f"  Homo subset too small ({len(homo)}) - skipped")
    else:
        print("  No parseable conditions - Homo skipped")

    print("\nDone.")


if __name__ == "__main__":
    main()
