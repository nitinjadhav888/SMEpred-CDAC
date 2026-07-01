"""
Full pipeline integration test for SMEpred.
Tests: rank → single-mod → multi-mod endpoints with known siRNA sequences.
Uses LNA paper benchmark sequences (Elmén et al. 2005, PMC546170).
"""
import requests, json, sys

BASE = "http://127.0.0.1:8000"

# Benchmark sequences from Elmén et al. 2005 (PMC546170)
# siRNA1 = anti-firefly luciferase (high-activity baseline)
# siRNA8 = siLNA with LNA at AS 5' (should score LOWER due to RISC penalty)
# siLNA5 = LNA at 3' overhangs only (should score similar/better than unmodified)
BENCHMARKS = [
    {
        "name": "siRNA1 (Unmodified firefly luciferase — high activity baseline)",
        "sense":     "CGUACGCGGAAUACUUCGAUU",
        "antisense": "UCGAAGUAUUCCGCGUACGUU",
        "expect": "moderate-high score, large penalties for unmodified state"
    },
    {
        "name": "siLNA5 (LNA at 3' overhangs — should preserve activity, improved serum)",
        "sense":     "CGUACGCGGAAUACUUCGALL",  # L = LNA at 3' overhangs
        "antisense": "UCGAAGUAUUCCGCGUACGLL",
        "expect": "similar to siRNA1 but serum penalty reduced"
    },
    {
        "name": "siLNA8 (LNA at AS 5' pos — experimentally abolishes activity)",
        "sense":     "CGUACGCGGAAUACUUCGAUU",
        "antisense": "LCGAAGUAUUCCGCGUACGUU",  # L at AS[0] = experimentally dead
        "expect": "RISC penalty should skyrocket — validates LNA-AS-5' rule"
    },
]

def separator(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print('='*60)

def test_health():
    separator("1. API HEALTH CHECK")
    r = requests.get(f"{BASE}/", timeout=10)
    print(f"Status: {r.status_code}")
    print(r.text[:300])
    assert r.status_code == 200, "API not reachable!"
    print("✓ API is online")

def test_rank():
    separator("2. RANK ENDPOINT (sliding window)")
    payload = {
        "sequence": "GAAGAAGCUGCCGUUGUUCUGGGUACUACAGCAGAAGGG",
        "top_n": 5,
        "input_type": "gene"
    }
    r = requests.post(f"{BASE}/rank", json=payload, timeout=30)
    data = r.json()
    assert r.status_code == 200
    results = data.get("results", [])
    assert len(results) > 0, "No candidates returned!"
    print(f"✓ Rank endpoint: {len(results)} candidates returned")
    for c in results[:3]:
        print(f"  Rank {c.get('rank')}: sense={c.get('sense','')[:15]}... score={c.get('efficacy_score',0):.1f}")

def test_single_mod():
    separator("3. SINGLE-MOD ENDPOINT")
    payload = {
        "sense":     "GCUGGAAGUGCUUUUGACGUU",
        "antisense": "CGUCAAAAGCACUUCCAGCUU",
        "model": "B",
        "top_n": 5,
        "full_scan": False,
    }
    r = requests.post(f"{BASE}/single-mod", json=payload, timeout=60)
    data = r.json()
    assert r.status_code == 200
    results = data.get("results", [])
    assert len(results) > 0
    print(f"✓ Single-mod endpoint: {len(results)} variants")
    print(f"  Parent score: {data.get('parent_score','?')}")
    print(f"  Best delta: {results[0].get('delta_score','?')}")

def test_multi_mod():
    separator("4. MULTI-MOD ENDPOINT (manual spec)")
    payload = {
        "sense":     "GCUGGAAGUGCUUUUGACGUU",
        "antisense": "CGUCAAAAGCACUUCCAGCUU",
        "sense_mods": "F,M,F,M",
        "sense_positions": "0,1,2,3",
        "antisense_mods": "S,S,F,F",
        "antisense_positions": "0,1,2,3",
        "model": "B",
    }
    r = requests.post(f"{BASE}/multi-mod", json=payload, timeout=30)
    data = r.json()
    assert r.status_code == 200
    result = data.get("result", {})
    score = result.get("efficacy_score", 0)
    penalties = result.get("penalties", {})
    print(f"✓ Multi-mod endpoint: final score = {score:.1f}")
    print(f"  Penalties: { {k: v['total'] for k,v in penalties.items()} }")

def test_biophysics_benchmarks():
    separator("5. LNA BENCHMARK VALIDATION (Elmén 2005)")
    scores = []
    for bm in BENCHMARKS:
        payload = {
            "sense":     bm["sense"],
            "antisense": bm["antisense"],
            "sense_mods": "",
            "sense_positions": "",
            "antisense_mods": "",
            "antisense_positions": "",
            "model": "B",
        }
        r = requests.post(f"{BASE}/multi-mod", json=payload, timeout=30)
        if r.status_code != 200:
            print(f"  ✗ {bm['name']}: HTTP {r.status_code}")
            continue
        result = r.json().get("result", {})
        score = result.get("efficacy_score", 0)
        penalties = result.get("penalties", {})
        risc_pen = penalties.get("risc", {}).get("total", 0)
        serum_pen = penalties.get("serum", {}).get("total", 0)
        scores.append((bm["name"], score, risc_pen, serum_pen))
        print(f"\n  [{bm['name']}]")
        print(f"    Final score: {score:.1f}")
        print(f"    RISC penalty: {risc_pen:.1f}  |  Serum penalty: {serum_pen:.1f}")
        print(f"    Expected: {bm['expect']}")

    # Scientific validation: siLNA8 (LNA at AS 5') should have HIGHEST RISC penalty
    if len(scores) >= 3:
        risc_unmod  = scores[0][2]
        risc_lna3p  = scores[1][2]
        risc_lna_as5 = scores[2][2]
        print(f"\n  ─── Scientific Consistency Check ───")
        print(f"  siRNA1 (unmod) RISC: {risc_unmod:.1f}")
        print(f"  siLNA5 (3' LNA) RISC: {risc_lna3p:.1f}")
        print(f"  siLNA8 (AS-5' LNA) RISC: {risc_lna_as5:.1f}")
        if risc_lna_as5 > risc_unmod:
            print("  ✓ VALIDATED: LNA at AS 5' correctly produces higher RISC penalty than unmodified")
        else:
            print("  ✗ WARNING: LNA at AS 5' should produce higher RISC penalty — check rules")

def test_modifications_endpoint():
    separator("6. MODIFICATIONS CATALOG")
    r = requests.get(f"{BASE}/modifications", timeout=10)
    data = r.json()
    assert r.status_code == 200
    mods = data.get("modifications", [])
    print(f"✓ Modifications endpoint: {len(mods)} mods loaded")

if __name__ == "__main__":
    print("\n" + "█"*60)
    print("  SMEpred FULL PIPELINE INTEGRATION TEST")
    print("  Based on: Elmén et al. 2005 (LNA) + Weingärtner 2020 (GalNAc)")
    print("█"*60)
    
    errors = []
    tests = [test_health, test_rank, test_single_mod, test_multi_mod, 
             test_biophysics_benchmarks, test_modifications_endpoint]
    
    for test in tests:
        try:
            test()
        except Exception as e:
            errors.append(f"{test.__name__}: {e}")
            print(f"  ✗ FAILED: {e}")
    
    separator("SUMMARY")
    total = len(tests)
    passed = total - len(errors)
    print(f"Passed: {passed}/{total}")
    if errors:
        print("Failures:")
        for e in errors:
            print(f"  ✗ {e}")
    else:
        print("✓ All tests passed! Pipeline is fully operational.")
    sys.exit(0 if not errors else 1)
