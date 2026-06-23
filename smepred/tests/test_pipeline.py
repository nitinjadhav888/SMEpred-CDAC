"""
tests/test_pipeline.py — Unit tests for each pipeline module.

These tests verify correctness without needing trained model files.
They test: parsing, siRNA generation, feature extraction, and modification engine.
"""

import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.parser import load_sequence, _normalize, _parse_fasta
from src.sirna_generator import generate_candidates, _reverse_complement
from src.features import extract_positional_features_batch, extract_batch_v4
from src.modification_engine import single_mod_scan, multimod_gen, _apply_mod


# ─── parser tests ─────────────────────────────────────────────────────────────

def test_normalize_dna_to_rna():
    seq = load_sequence("ATGCATGCATGCATGCATGCA")  # DNA T → RNA U
    assert "T" not in seq
    assert "U" in seq

def test_normalize_already_rna():
    seq = load_sequence("AUGCAUGCAUGCAUGCAUGCA")
    assert seq == "AUGCAUGCAUGCAUGCAUGCA"

def test_parse_fasta_inline():
    fasta = ">gene1\nAUGCAUGCAUGCAUGCAUGCA"
    seq = load_sequence(fasta)
    assert seq == "AUGCAUGCAUGCAUGCAUGCA"

def test_sequence_too_short():
    try:
        load_sequence("AUGCAU")
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "too short" in str(e).lower()


# ─── siRNA generation tests ───────────────────────────────────────────────────

def test_reverse_complement():
    sense = "AUGCAUGCAUGCAUGCAUGCA"
    rc = _reverse_complement(sense)
    assert len(rc) == len(sense)
    # first base of sense is A → last base of RC must be U
    assert rc[-1] == "U"

def test_candidate_count():
    seq = "A" * 100  # 100-nt sequence
    candidates = generate_candidates(seq)
    assert len(candidates) == 80  # 100 - 21 + 1

def test_candidate_length():
    seq = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUG"  # 31-nt
    candidates = generate_candidates(seq)
    for c in candidates:
        assert len(c.sense) == 21
        assert len(c.antisense) == 21


# ─── feature tests ────────────────────────────────────────────────────────────

def test_positional_features_shape():
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    X = extract_positional_features_batch([s], [a])
    assert X.shape == (1, 1467), f"Expected (1, 1467), got {X.shape}"

def test_positional_features_multi_batch():
    s1 = "GCAGCACGACUUCUUCAAGUU"
    a1 = "CUUGAAGAAGUCGUGCUGCUU"
    s2 = "AUGCAUGCAUGCAUGCAUGCA"
    a2 = "UGCAUGCAUGCAUGCAUGCAU"
    X = extract_positional_features_batch([s1, s2], [a1, a2])
    assert X.shape == (2, 1467)

def test_positional_features_detects_mod():
    base_s = "GCAGCACGACUUCUUCAAGUU"
    base_a = "CUUGAAGAAGUCGUGCUGCUU"
    mod_s  = "FCAGCACGACUUCUUCAAGUU"  # F at position 1
    X_base = extract_positional_features_batch([base_s], [base_a], [base_s], [base_a])
    X_mod  = extract_positional_features_batch([mod_s], [base_a], [base_s], [base_a])
    # The two vectors should differ (mod introduces flags)
    assert not np.allclose(X_base, X_mod), "modified and unmodified should differ"

def test_v4_features_shape():
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGTCGUGCUGCUU"
    X = extract_batch_v4([s], [a])
    assert X.shape == (1, 214), f"Expected (1, 214), got {X.shape}"

def test_v4_features_multi_batch():
    s1 = "GCAGCACGACUUCUUCAAGUU"
    a1 = "CUUGAAGAAGUCGUGCUGCUU"
    s2 = "AUGCAUGCAUGCAUGCAUGCA"
    a2 = "UGCAUGCAUGCAUGCAUGCAU"
    X = extract_batch_v4([s1, s2], [a1, a2])
    assert X.shape == (2, 214)

def test_v4_features_onehot_sums_to_one():
    s = "AUGCAUGCAUGCAUGCAUGCA"
    a = "UGCAUGCAUGCAUGCAUGCAU"
    X = extract_batch_v4([s], [a])
    # First 84 features = position one-hot (21 pos × 4 bases)
    onehot = X[0, :84].reshape(21, 4)
    row_sums = onehot.sum(axis=1)
    assert np.allclose(row_sums, 1.0), "each position row should sum to 1"


# ─── modification engine tests ────────────────────────────────────────────────

def test_apply_mod_position():
    seq = "AUGCAUGCAUGCAUGCAUGCA"
    modified = _apply_mod(seq, 3, "F")  # position 3 (1-based) = index 2
    assert modified[2] == "F"
    assert modified[:2] == seq[:2]
    assert modified[3:] == seq[3:]

def test_single_mod_scan_count():
    sense = "GCAGCACGACUUCUUCAAGUU"
    antisense = "CUUGAAGAAGUCGUGCUGCUU"
    variants = single_mod_scan(sense, antisense)
    # 30 mods × 21 positions × 2 strands = 1260
    assert len(variants) == 1260, f"Expected 1260, got {len(variants)}"

def test_single_mod_all_positions_covered():
    sense = "GCAGCACGACUUCUUCAAGUU"
    antisense = "CUUGAAGAAGUCGUGCUGCUU"
    variants = single_mod_scan(sense, antisense)
    positions = {v.mod_position for v in variants}
    assert positions == set(range(1, 22)), "All 21 positions should be covered"

def test_multimod_gen_basic():
    sense = "GCAGCACGACUUCUUCAAGUU"
    antisense = "CUUGAAGAAGUCGUGCUGCUU"
    result = multimod_gen(
        sense, antisense,
        sense_mods="F", sense_positions="2,5",
    )
    assert result.sense[1] == "F"   # position 2 (0-indexed: 1)
    assert result.sense[4] == "F"   # position 5 (0-indexed: 4)

def test_multimod_gen_multi_type():
    sense = "GCAGCACGACUUCUUCAAGUU"
    antisense = "CUUGAAGAAGUCGUGCUGCUU"
    result = multimod_gen(
        sense, antisense,
        sense_mods="F,,M", sense_positions="2,,5",
    )
    assert result.sense[1] == "F"   # position 2
    assert result.sense[4] == "M"   # position 5


# ─── biophysics tests ────────────────────────────────────────────────────────

def test_biophysics_nuclease_penalty_ps():
    from src.biophysics import nuclease_penalty, adjusted_efficacy_score
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = nuclease_penalty(s, a, s, a)
    assert p0 >= 0
    # PS at sense 5' terminus reduces penalty
    mod_s = "S" + s[1:]
    p_ps = nuclease_penalty(mod_s, a, s, a)
    assert p_ps <= p0, f"PS at terminus should reduce nuclease penalty: {p_ps} > {p0}"
    # Adjusted score should be lower than raw
    adj, _, _ = adjusted_efficacy_score(100, s, a, s, a)
    assert 0 <= adj <= 100
    assert adj < 100, "Penalties should reduce raw score"

def test_biophysics_immuno_uridine_penalty():
    from src.biophysics import immuno_penalty
    s = "U" * 21
    a = "U" * 21
    # All unmodified U → high penalty
    p0 = immuno_penalty(s, a, s, a)
    # Modify all U with M → lower penalty
    mod_s = "M" * 21
    mod_a = "M" * 21
    p_mod = immuno_penalty(mod_s, mod_a, s, a)
    assert p_mod <= p0, "Modifying all U should reduce immuno penalty"

def test_biophysics_risc_5p_reduces_penalty():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = risc_penalty(s, a, s, a)
    mod_a = "1" + a[1:]  # 5'-P on antisense
    p_5p = risc_penalty(s, mod_a, s, a)
    assert p_5p <= p0, "5'-P should reduce RISC penalty"

def test_biophysics_risc_moe_penalty():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = risc_penalty(s, a, s, a)
    # MOE at position 6 (index 5) in antisense seed
    mod_a = a[:5] + "E" + a[6:]
    p_moe = risc_penalty(s, mod_a, s, a)
    assert p_moe > p0, "MOE in seed should increase RISC penalty"

def test_biophysics_risc_gna_disruptive():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = risc_penalty(s, a, s, a)
    # GNA at position 3 (index 2) — disruptive zone (pos 2-5)
    mod_a = a[:2] + "8" + a[3:]
    p_gna = risc_penalty(s, mod_a, s, a)
    assert p_gna > p0, "GNA at pos 3 should increase RISC penalty"

def test_biophysics_risc_gna_beneficial():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    # GNA at position 7 (index 6) — beneficial zone (pos 6-8, Schlegel 2022 ESC+)
    gna_a = a[:6] + "8" + a[7:]
    p_gna = risc_penalty(s, gna_a, s, a)
    # Same position with 2'-OMe instead — neutral mod, no bonus
    ome_a = a[:6] + "M" + a[7:]
    p_ome = risc_penalty(s, ome_a, s, a)
    # GNA should be lower than 2'-OMe at same position
    assert p_gna < p_ome, f"GNA at pos 7 should have lower penalty than 2'-OMe: {p_gna} >= {p_ome}"

def test_biophysics_risc_una_exempt():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    # UNA at position 7 (index 6) — should be exempt from seed penalty (Bramsen 2010)
    una_a = a[:6] + "6" + a[7:]
    p_una = risc_penalty(s, una_a, s, a)
    # 2'-OMe at same position — NOT exempt, should pay seed mod penalty
    ome_a = a[:6] + "M" + a[7:]
    p_ome = risc_penalty(s, ome_a, s, a)
    assert p_una < p_ome, f"UNA at pos 7 should have lower penalty than 2'-OMe: {p_una} >= {p_ome}"

def test_biophysics_risc_ena_penalty():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = risc_penalty(s, a, s, a)
    # ENA at position 5 (index 4) — seed zone (pos 2-8)
    mod_a = a[:4] + "Y" + a[5:]
    p_ena = risc_penalty(s, mod_a, s, a)
    assert p_ena > p0, "ENA in seed should increase RISC penalty"

def test_biophysics_risc_tna_penalty():
    from src.biophysics import risc_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = risc_penalty(s, a, s, a)
    # TNA at position 4 (index 3) — seed body (pos 2-6)
    mod_a = a[:3] + "9" + a[4:]
    p_tna = risc_penalty(s, mod_a, s, a)
    assert p_tna > p0, "TNA in seed body should increase RISC penalty"

def test_biophysics_risc_missing_2f_penalty():
    from src.biophysics import risc_penalty
    # Use sequence with only 6 U/C and test 0 vs 4 vs 6 2'-F
    s = "AAAAAAAAAAAAAAAAAAAAA"
    a_base = "AAAAAAAUGCAUGCAUGCAU"  # ~6 U/C total
    p_no_f = risc_penalty(s, a_base, s, a_base)
    # 4 F on pyrimidines → ~67% coverage (~4/6) → no penalty, but total_mods=4
    mod_a = list(a_base)
    count = 0
    for i, b in enumerate(mod_a):
        if b in 'UC' and count < 4:
            mod_a[i] = 'F'
            count += 1
    mod_a = ''.join(mod_a)
    p_with_f = risc_penalty(s, mod_a, s, a_base)
    # 2 F on pyrimidines → ~33% coverage → partial penalty, total_mods=2
    mod_a_partial = list(a_base)
    count = 0
    for i, b in enumerate(mod_a_partial):
        if b in 'UC' and count < 2:
            mod_a_partial[i] = 'F'
            count += 1
    mod_a_partial = ''.join(mod_a_partial)
    p_partial = risc_penalty(s, mod_a_partial, s, a_base)
    assert p_with_f < p_partial, "Fuller 2'-F coverage should have lower penalty than partial"
    assert p_partial < p_no_f, "Partial 2'-F coverage should have lower penalty than none"

def test_biophysics_risc_exotic_penalty():
    from src.biophysics import risc_penalty
    s = "AAAAAAAAAAAAAAAAAAAAA"
    a = "AAAAAAAAAAAAAAAAAAAAA"
    # Benzyl (B) in guide strand → exotic penalty
    mod_b = list(a); mod_b[5] = "B"; mod_b = "".join(mod_b)
    p_b = risc_penalty(s, mod_b, s, a)
    # 2'-OMe (M) in guide strand → no exotic penalty
    mod_m = list(a); mod_m[5] = "M"; mod_m = "".join(mod_m)
    p_m = risc_penalty(s, mod_m, s, a)
    assert p_b > p_m, "Benzyl should have higher penalty than 2'-OMe"
    # Multiple exotic mods
    mod_bb = list(a); mod_bb[5] = "B"; mod_bb[10] = "J"; mod_bb = "".join(mod_bb)
    p_bb = risc_penalty(s, mod_bb, s, a)
    assert p_bb > p_b, "Multiple exotic mods should increase penalty"

def test_biophysics_thermo_low_gc_penalty():
    from src.biophysics import thermo_penalty
    # Ideal GC = 43% (9/21) → within 35-50% sweet spot
    ideal = "AUGCAUGCAAUGCAUGCAUGC"
    ideal_a = "GCAUGCAUUGCAUGCAUGCAU"
    p_ideal = thermo_penalty(ideal, ideal_a, ideal, ideal_a)
    # 0% GC → extreme, penalty should be higher
    low_gc = "AU" * 10 + "AA"
    low_gc_a = "UA" * 10 + "UU"
    p_low = thermo_penalty(low_gc, low_gc_a, low_gc, low_gc_a)
    assert p_low > p_ideal, f"Low GC should have higher penalty than ideal: {p_low} <= {p_ideal}"

def test_biophysics_serum_ps_reduces_penalty():
    from src.biophysics import serum_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    p0 = serum_penalty(s, a, s, a)
    mod_a = "S" + a[1:20] + "S"
    p_ps = serum_penalty(s, mod_a, s, a)
    assert p_ps <= p0, "PS at AS termini should reduce serum penalty"

def test_biophysics_adjusted_score_range():
    from src.biophysics import adjusted_efficacy_score, nuclease_penalty, immuno_penalty, risc_penalty, thermo_penalty, serum_penalty
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    for fn in [nuclease_penalty, immuno_penalty, risc_penalty, thermo_penalty, serum_penalty]:
        p = fn(s, a, s, a)
        assert 0 <= p <= 60, f"{fn.__name__} returned {p} outside expected range"
    adj, penalties, total = adjusted_efficacy_score(80, s, a, s, a)
    assert 0 <= adj <= 100
    assert set(penalties.keys()) == {"nuclease", "immuno", "risc", "thermo", "serum"}
    assert total >= 0


# ─── run tests ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    tests = [
        test_normalize_dna_to_rna,
        test_normalize_already_rna,
        test_parse_fasta_inline,
        test_sequence_too_short,
        test_reverse_complement,
        test_candidate_count,
        test_candidate_length,
        test_positional_features_shape,
        test_positional_features_multi_batch,
        test_positional_features_detects_mod,
        test_v4_features_shape,
        test_v4_features_multi_batch,
        test_v4_features_onehot_sums_to_one,
        test_apply_mod_position,
        test_single_mod_scan_count,
        test_single_mod_all_positions_covered,
        test_multimod_gen_basic,
        test_multimod_gen_multi_type,
        # Biophysics
        test_biophysics_nuclease_penalty_ps,
        test_biophysics_immuno_uridine_penalty,
        test_biophysics_risc_5p_reduces_penalty,
        test_biophysics_risc_moe_penalty,
        test_biophysics_risc_gna_disruptive,
        test_biophysics_risc_gna_beneficial,
        test_biophysics_risc_una_exempt,
        test_biophysics_risc_ena_penalty,
        test_biophysics_risc_tna_penalty,
        test_biophysics_risc_missing_2f_penalty,
        test_biophysics_risc_exotic_penalty,
        test_biophysics_thermo_low_gc_penalty,
        test_biophysics_serum_ps_reduces_penalty,
        test_biophysics_adjusted_score_range,
    ]

    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {t.__name__}: {e}")
            failed += 1

    print(f"\n{passed} passed, {failed} failed out of {len(tests)} tests.")
    sys.exit(0 if failed == 0 else 1)
