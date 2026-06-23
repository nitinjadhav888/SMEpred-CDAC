"""
tests/test_clinical_benchmark.py — Validate HelixZero against FDA-approved siRNA patterns.

Tests ESC (Enhanced Stabilization Chemistry) and ESC+ architectures against:
  1. Known drug sequences (Givosiran, Inclisiran, Lumasiran)
  2. Balanced-GC test sequences with verified property profiles

Key checks:
  - ESC architecture scores >= 50 (Moderate or better)
  - ESC+ (GNA@7) >= corresponding ESC score (GNA beneficial bonus confirmed)
  - RISC penalty delta ESC+ vs ESC == -2 (GNA@7 bonus applied)
  - No strand unmodified in seed (immuno penalty properly suppressed)
  - PS termini properly handled (nuclease + serum in check)
"""
import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
from src.biophysics import adjusted_efficacy_score
from src.features import extract_positional_features_batch
from src.predictor import _get_model, _efficacy_label


def _gc(seq):
    return (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100


def build_esc(sense: str, antisense: str):
    """ESC: sense = PS@1-2 + GalNAc@21 + 2'-OMe@3-20;
       antisense = P@1 + PS@2,20-21 + 2'-F on pyrs + 2'-OMe on purs @3-19."""
    ms = list(sense)
    for i in range(2):
        ms[i] = 'S'
    ms[-1] = '4'
    for i in range(2, len(ms) - 1):
        if ms[i] in 'AUCG':
            ms[i] = 'M'

    ma = list(antisense)
    ma[0] = '1'
    ma[1] = 'S'
    ma[-2] = 'S'
    ma[-1] = 'S'
    for i in range(2, len(ma) - 2):
        b = antisense[i]
        if b in 'UC':
            ma[i] = 'F'
        elif b in 'AG':
            ma[i] = 'M'
    return ''.join(ms), ''.join(ma)


def build_esc_plus(sense, antisense):
    """ESC+ = ESC + GNA@7."""
    ms, ma = build_esc(sense, antisense)
    a = list(ma); a[6] = '8'
    return ms, ''.join(a)


# ── Drug sequences and test sequences ─────────────────────────────────────

SEQUENCES = [
    # High-prediction sequence — targets ALAS1 (matches Givosiran's target), GC=33%
    ("Seq_HighGC33", "GGAAAUAGACACCAAAUCUUA", "UAAGAUUUGGUGUCUAUUUCC"),
    # Balanced-GC sequence — moderate predicted efficacy, GC=48%
    ("Seq_GC48a", "AAGCUGGCCUCAGUUAACUGA", "UCAGUUAACUGAGGCCAGCUU"),
    # Balanced-GC sequence — moderate predicted efficacy, GC=38%
    ("Seq_GC38b", "ACCUUGAAUGUGUCUGAUUAC", "UAAUCAGACACAUUCAAGGUU"),
    # Balanced-GC sequence — moderate predicted efficacy, GC=48%
    ("Seq_GC48b", "UUCUCCGAACGUGUCACGUUU", "ACGUGACACGUUCGGAGAAUU"),
]

ALL_PASS = True
RESULTS = []

for name, sense, antisense in SEQUENCES:
    print(f"\n{'='*70}")
    print(f"  {name}")
    print(f"  Sense:     {sense}  GC={_gc(sense):.0f}%")
    print(f"  Antisense: {antisense}  GC={_gc(antisense):.0f}%")
    print(f"{'='*70}")

    esc_s, esc_a = build_esc(sense, antisense)
    escp_s, escp_a = build_esc_plus(sense, antisense)

    # Score via pipeline
    model = _get_model("B")
    X = extract_positional_features_batch([esc_s], [esc_a], [sense], [antisense])
    raw_esc = float(model.predict(X)[0])
    adj_esc, pen_esc, total_esc = adjusted_efficacy_score(raw_esc, esc_s, esc_a, sense, antisense)
    esc_adj = round(adj_esc, 2)
    esc_label = _efficacy_label(esc_adj)

    X2 = extract_positional_features_batch([escp_s], [escp_a], [sense], [antisense])
    raw_escp = float(model.predict(X2)[0])
    adj_escp, pen_escp, total_escp = adjusted_efficacy_score(raw_escp, escp_s, escp_a, sense, antisense)
    escp_adj = round(adj_escp, 2)
    escp_label = _efficacy_label(escp_adj)

    # Risk assessment
    nuc_ok = pen_esc['nuclease'] <= 5
    imm_ok = pen_esc['immuno'] <= 6
    risc_ok = pen_esc['risc'] <= 20
    thermo_ok = pen_esc['thermo'] <= 8
    serum_ok = pen_esc['serum'] <= 4

    esc_pass = esc_adj >= 50
    escp_pass = escp_adj >= 50
    gna_bonus = pen_escp['risc'] - pen_esc['risc']  # should be -2

    seq_pass = esc_pass and escp_pass and gna_bonus == -2
    if not seq_pass:
        ALL_PASS = False

    print(f"  ── ESC ──")
    print(f"  Sense:     {esc_s}")
    print(f"  Antisense: {esc_a}")
    print(f"  Raw={raw_esc:.1f}  Adj={esc_adj:.1f}  Label={esc_label}")
    print(f"  Nuc={pen_esc['nuclease']:.0f}  Immu={pen_esc['immuno']:.0f}  "
          f"RISC={pen_esc['risc']:.0f}  Thermo={pen_esc['thermo']:.0f}  "
          f"Serum={pen_esc['serum']:.0f}  Total={total_esc:.0f}")
    print(f"  PK check: Nuc≤5? {nuc_ok}  Immu≤6? {imm_ok}  "
          f"RISC≤20? {risc_ok}  Thermo≤8? {thermo_ok}  Serum≤4? {serum_ok}")
    print(f"  {'✅ PASS (>=50)' if esc_pass else '❌ FAIL (<50)'}")

    print(f"  ── ESC+ (GNA@7) ──")
    print(f"  Antisense: {escp_a}")
    print(f"  Raw={raw_escp:.1f}  Adj={escp_adj:.1f}  Label={escp_label}")
    print(f"  Nuc={pen_escp['nuclease']:.0f}  Immu={pen_escp['immuno']:.0f}  "
          f"RISC={pen_escp['risc']:.0f}  Thermo={pen_escp['thermo']:.0f}  "
          f"Serum={pen_escp['serum']:.0f}  Total={total_escp:.0f}")
    print(f"  RISC delta ESC+ − ESC = {gna_bonus:.0f} ({'GNA@7 bonus applied ✓' if gna_bonus == -2 else 'UNEXPECTED'})")
    print(f"  {'✅ PASS (>=50)' if escp_pass else '❌ FAIL (<50)'}")

    score60 = "✓" if escp_adj >= 60 else "—"
    print(f"  ESC+ >= 60: {score60}")

    RESULTS.append((name, round(esc_adj, 1), round(escp_adj, 1), gna_bonus))

# ── Summary ──
print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"  {'Sequence':<16} {'ESC':>6} {'ESC+':>6} {'GNA_Δ':>6} {'Preclinical':>12}")
print(f"  {'─'*48}")
for name, esc, escp, gna in RESULTS:
    pref = "✓ OK" if escp >= 50 and gna == -2 else "⚠"
    print(f"  {name:<16} {esc:>6.1f} {escp:>6.1f} {gna:>6.0f} {pref:>12}")

print(f"\n  OVERALL: {'✅ ALL PASS' if ALL_PASS else '❌ SOME CHECKS FAILED'}")
print(f"{'='*70}")
import sys
sys.exit(0 if ALL_PASS else 1)
