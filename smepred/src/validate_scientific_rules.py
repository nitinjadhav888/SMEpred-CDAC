import os
from offtarget import get_offtarget_engine

def print_report(name, report):
    print(f"\n{'='*60}")
    print(f"TEST: {name}")
    print(f"Overall Safety Score: {report['overallSafetyScore']}")
    print(f"Status: {report['status']}")
    print("Risk Factors:")
    for r in report['riskFactors']:
        print(f"  - {r}")
    print("Safety Notes:")
    for s in report['safetyNotes']:
        print(f"  - {s}")
    print(f"{'='*60}")

def run_validations():
    engine = get_offtarget_engine()
    
    # 1. 15-mer Transcriptome Match
    # Let's find a 15-mer from the transcriptome and make it the antisense
    # The transcriptome might be huge, let's just make a dummy string if we can't find one
    # Wait, the engine has `self.sequence`. Let's just pick the first 15 chars from it!
    toxic_15mer = engine.sequence[:15] if engine.sequence else "ACGTACGTACGTACG"
    antisense_15 = toxic_15mer + "AAAAAA" # 21 nt
    sense_15     = "UUUUUU" + toxic_15mer[::-1] # dummy
    report_1 = engine.validate_safety(sense_15, antisense_15, antisense_15, sense_15)
    print_report("1. Transcriptome-Wide Off-Target Screening (15-mer toxicity)", report_1)
    
    # 2. Thermodynamic Asymmetry
    # Make Antisense 5' end stronger (G/C rich) and 3' end weaker (A/U rich)
    as_thermo = "GGCGCGAAAAAUUUUUUUAAA"
    ss_thermo = "UUUAAAAAAAUUUUUCGCGCC"
    report_2 = engine.validate_safety(ss_thermo, as_thermo, as_thermo, ss_thermo)
    print_report("2. True Thermodynamic Asymmetry Profiling (ΔG)", report_2)
    
    # 3. TLR7/TLR8 Motif Masking (Unmasked vs Masked)
    # UGGC motif in antisense
    as_tlr = "AAUUGGCUAAAAAAAAAAAAA"
    ss_tlr = "UUUUUUUUUUUUUAGCCAAUU"
    
    # UNMASKED (no M)
    mod_as_unmasked = as_tlr
    report_3a = engine.validate_safety(ss_tlr, as_tlr, mod_as_unmasked, ss_tlr)
    print_report("3a. Toll-Like Receptor Motif (Unmasked UGGC)", report_3a)
    
    # MASKED (M at position of U in UGGC)
    # as_tlr: A A U U G G C U
    # idx:    0 1 2 3 4 5 6
    # motif UGGC starts at index 2
    mod_as_masked = "AAMUGGCUAAAAAAAAAAAAA" # M at pos 2
    report_3b = engine.validate_safety(ss_tlr, as_tlr, mod_as_masked, ss_tlr)
    print_report("3b. Toll-Like Receptor Motif (Masked UGGC)", report_3b)
    
    # 4. PK & Delivery Conjugate Validation (GalNAc)
    # No GalNAc
    report_4a = engine.validate_safety(ss_tlr, as_tlr, mod_as_masked, ss_tlr)
    print_report("4a. PK Validation (Missing GalNAc 4 token)", report_4a)
    
    # With GalNAc (4 token in sense)
    mod_ss_galnac = "4UUUUUUUUUUUUAGCCAAUU"
    report_4b = engine.validate_safety(ss_tlr, as_tlr, mod_as_masked, mod_ss_galnac)
    print_report("4b. PK Validation (Present GalNAc 4 token)", report_4b)

if __name__ == "__main__":
    run_validations()
