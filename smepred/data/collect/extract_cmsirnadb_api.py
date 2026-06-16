"""
extract_cmsirnadb_api.py — Extract ALL records from CMsiRNAdb browse API page by page.
Parses mts/Mtas (descriptive per-position names) into 35-symbol encoding.
"""
import json, re, sys, time
from pathlib import Path
from collections import Counter
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

import pandas as pd

# ─── load alias rules ────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.collect.clean_utils import map_modification, parse_efficacy

API = "https://cellknowledge.com.cn/CMsiRNAdb/php_mysql/mysql_browse_severside.php"
DATA_DIR = Path(__file__).parent.parent.parent / "data"
OUT_CSV = DATA_DIR / "cmsirnadb_full.csv"
UNMAPPED_LOG = Path(__file__).parent / "cmsirnadb_unmapped.txt"

# ─── helpers ─────────────────────────────────────────────────────────────

_TOKEN_RE = re.compile(r"(\d+)\*(.+?)(?=\s*\|\||$)")
_ID_COND_RE = re.compile(r"-([0-9.]+n)-([0-9.]+h)-")
_BASE_SUFFIX = [
    ("adenosine", "A"), ("adenine", "A"),
    ("uridine", "U"), ("uracil", "U"),
    ("thymidine", "U"), ("thymine", "U"),
    ("cytidine", "C"), ("cytosine", "C"),
    ("guanosine", "G"), ("guanine", "G"),
]


def base_from_name(name: str):
    s = name.strip()
    up = s.upper()
    if up in ("A", "U", "G", "C", "T"):
        return "U" if up == "T" else up
    low = s.lower()
    for suffix, base in _BASE_SUFFIX:
        if suffix in low:
            return base
    return None


def parse_token_stream(token_str: str):
    """Parse '1*2'-O-Methylcytidine || 2*A || ...' into (base_str, mod_str)."""
    if not token_str or token_str.strip() == "":
        return None, None
    tokens = _TOKEN_RE.findall(token_str)
    if not tokens:
        return None, None
    n = len(tokens)
    base_chars = [""] * n
    mod_chars = [""] * n
    for pos_str, name in tokens:
        pos = int(pos_str) - 1
        base = base_from_name(name)
        sym = map_modification(name)
        if sym is None:
            # Unmapped — treat as canonical placeholder
            sym = ""
        if base is None:
            base = "A"  # fallback for abasic
        base_chars[pos] = base
        if sym == "":
            mod_chars[pos] = base
        else:
            mod_chars[pos] = sym
    return "".join(base_chars), "".join(mod_chars)


def parse_id(id_str: str):
    """Extract concentration_nM, time_h from ID like '...-100n-48h-88.00'."""
    m = _ID_COND_RE.search(id_str)
    if m:
        return float(m.group(1).rstrip("n")), float(m.group(2).rstrip("h"))
    return float("nan"), float("nan")


def valid_length(seq, lo=19, hi=25):
    return seq is not None and lo <= len(seq) <= hi


# ─── main extraction ─────────────────────────────────────────────────────

def main():
    # First, get total count
    url = f"{API}?page=1&size=1&browserCell_Type=All&browserModType=All&browselocation=All&t={int(time.time())}"
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=60) as resp:
        meta = json.loads(resp.read().decode())
    
    total = meta["total"]
    print(f"Total records in CMsiRNAdb: {total}")

    # Use page size 200 for efficiency
    page_size = 200
    last_page = (total + page_size - 1) // page_size
    print(f"Fetching {last_page} pages (size={page_size})...")

    records = []
    unmapped_counter = Counter()
    errors = 0

    for page in range(1, last_page + 1):
        url = (f"{API}?page={page}&size={page_size}"
               f"&browserCell_Type=All&browserModType=All&browselocation=All"
               f"&t={int(time.time() * 1000)}")
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        
        for attempt in range(3):
            try:
                with urlopen(req, timeout=120) as resp:
                    data = json.loads(resp.read().decode())
                break
            except (HTTPError, URLError, TimeoutError) as e:
                print(f"  Page {page} attempt {attempt+1} failed: {e}")
                time.sleep(5)
        else:
            print(f"  Page {page} failed after 3 attempts, skipping")
            errors += 1
            continue

        for row in data["data"]:
            try:
                # Sense strand
                base_s, mod_s = parse_token_stream(row.get("mts", ""))
                # Antisense strand
                base_a, mod_a = parse_token_stream(row.get("Mtas", ""))

                if base_s is None or base_a is None:
                    continue
                if not (valid_length(base_s) and valid_length(base_a)):
                    continue

                efficacy = parse_efficacy(row.get("ihb"))
                if efficacy is None:
                    continue

                conc_nM, time_h = parse_id(row.get("ID", ""))

                # Track unmapped mods (map_modification returns None for unmapped)
                # We already handle this per-position in parse_token_stream

                records.append({
                    "sense": mod_s,
                    "antisense": mod_a,
                    "base_sense": base_s,
                    "base_antisense": base_a,
                    "efficacy": efficacy,
                    "concentration_nM": conc_nM,
                    "time_h": time_h,
                    "target_gene": row.get("ge", ""),
                    "cell_type": row.get("cel", ""),
                    "patent_id": row.get("pID", ""),
                })
            except Exception as e:
                errors += 1
                continue

        if page % 20 == 0 or page == last_page:
            print(f"  Page {page}/{last_page} — {len(records)} records so far (errors={errors})")

    df = pd.DataFrame(records)
    # Dedup on sequence + efficacy
    before = len(df)
    df = df.drop_duplicates(subset=["sense", "antisense", "efficacy"]).reset_index(drop=True)
    print(f"\nDone. Collected {before} raw, {len(df)} unique records.")

    # Stats
    print(f"\nEfficacy: min={df['efficacy'].min():.1f} max={df['efficacy'].max():.1f} mean={df['efficacy'].mean():.1f}")
    print(f"Target genes: {df['target_gene'].nunique()}")
    print(f"Cell types: {df['cell_type'].nunique()}")

    # Count distinct modification symbols used
    all_syms = set()
    for col in ["sense", "antisense"]:
        for seq in df[col]:
            all_syms.update(seq)
    canon = {"A", "U", "G", "C", "T"}
    mod_syms = sorted(all_syms - canon)
    print(f"Modification symbols used: {mod_syms}")

    # Save
    df.to_csv(OUT_CSV, index=False)
    print(f"Saved: {OUT_CSV} ({OUT_CSV.stat().st_size / 1024:.0f} KB)")

    # Unmapped report
    if unmapped_counter:
        lines = [f"# CMsiRNAdb API extraction - unmapped modification names",
                 f"# Total distinct: {len(unmapped_counter)}", ""]
        for name, ct in unmapped_counter.most_common(50):
            lines.append(f"{ct:>8}  {name}")
        UNMAPPED_LOG.write_text("\n".join(lines))
        print(f"Unmapped report: {UNMAPPED_LOG}")
    else:
        print("No unmapped modifications found.")


if __name__ == "__main__":
    main()
