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
from src.features import mnc_vector, mnc_full, bin_vector, features_model_a, features_model_b, features_model_c
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

def test_mnc_sums_to_one():
    seq = "AUGCAUGCAUGCAUGCAUGCA"
    vec = mnc_vector(seq)
    total = vec.sum()
    assert abs(total - 1.0) < 1e-5, f"MNC should sum to 1, got {total}"

def test_mnc_full_shape():
    sense = "AUGCAUGCAUGCAUGCAUGCA"
    antisense = "UGCAUUGCAUGCAUGCAUGCU"
    vec = mnc_full(sense, antisense)
    assert vec.shape == (70,), f"Expected (70,), got {vec.shape}"

def test_model_a_feature_shape():
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    f = features_model_a(s, a)
    assert f.shape == (140,)   # 70 base + 70 modified

def test_model_b_feature_shape():
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    f = features_model_b(s, a)
    expected = 140 + 35 * 13  # 140 + 455 = 595
    assert f.shape == (expected,), f"Expected ({expected},), got {f.shape}"

def test_model_c_feature_shape():
    s = "GCAGCACGACUUCUUCAAGUU"
    a = "CUUGAAGAAGUCGUGCUGCUU"
    f = features_model_c(s, a)
    expected = 140 + 35 * 8   # 140 + 280 = 420
    assert f.shape == (expected,), f"Expected ({expected},), got {f.shape}"

def test_model_a_base_differs_from_modified():
    # when base != modified, the two MNC halves should differ
    base_s, base_a = "GCAGCACGACUUCUUCAAGUU", "CUUGAAGAAGUCGUGCUGCUU"
    mod_s = "MCAGCACGACUUCUUCAAGUU"   # one 2'-OMe at position 1
    f = features_model_a(mod_s, base_a, base_sense=base_s, base_antisense=base_a)
    assert f.shape == (140,)
    assert not (f[:70] == f[70:]).all(), "base and modified halves should differ"

def test_bin_vector_binary():
    seq = "AUGCAUGCAUGCAUGCAUGCA"
    vec = bin_vector(seq)
    unique = set(vec.tolist())
    assert unique <= {0.0, 1.0}, "BIN vector should be binary"


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
        test_mnc_sums_to_one,
        test_mnc_full_shape,
        test_model_a_feature_shape,
        test_model_b_feature_shape,
        test_model_c_feature_shape,
        test_model_a_base_differs_from_modified,
        test_bin_vector_binary,
        test_apply_mod_position,
        test_single_mod_scan_count,
        test_single_mod_all_positions_covered,
        test_multimod_gen_basic,
        test_multimod_gen_multi_type,
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
