# Biophysical Penalties — Complete Reference

Five orthogonal domains adjust the raw LightGBM efficacy score downward (or upward in the case of RISC and GalNAc bonuses). Each domain targets a distinct biological mechanism, and the modules are designed to be **strictly non-overlapping** — no biological feature is penalized by more than one module.

> **Literature update — July 2026:** Four new rules implemented from peer-reviewed experimental data (Weingärtner et al. 2020, *Mol Ther Nucleic Acids*; Sakamuri et al. 2020, *ChemBioChem*). See sections 1 and 5 below.

---

## Adjustment Formula

```
calculate_adjusted_efficacy = max(0.0, min(100.0, raw_score − 0.70 × total_penalty))
```

The 0.70 factor is empirically calibrated so that:
- Fully unmodified siRNA: adjusted ≈ 15–25
- Best single-mod candidates: adjusted ≈ 35–60
- Top multi-mod clinical designs (ESC/ESC+): adjusted ≥ 50
- Dual-terminal GalNAc optimised designs: can score ≥ 60

---

## 1. Nuclease Penalty (Range: 0–20)

**Biological target**: Endonuclease stability (RNase A-family, 2'-5' oligoadenylate synthetase).

**Orthogonality note**: Does NOT check termini protection (that is the serum penalty's domain). Checks PS backbone coverage, PS positional distribution, and 2'-modification density.

| Condition | Penalty | Rationale | Citation |
|-----------|---------|-----------|----------|
| PS backbone count == 0 | +5 | No phosphorothioate whatsoever → rapid endonuclease cleavage | Braasch & Corey 2004 |
| PS backbone count < 3 | +3 | Minimal backbone protection → moderate susceptibility | Braasch & Corey 2004 |
| **AS terminal PS <2 at pos 0,1,20,21** | **+2** | **Alnylam clinical design requires ≥2 AS terminal PS — suboptimal pattern** | **Sakamuri et al. 2020** |
| **Sense terminal PS missing (pos 0 or 1)** | **+1** | **Alnylam clinical design requires ≥1 sense terminal PS** | **Sakamuri et al. 2020** |
| 2'-mod density < 20% | +4 | <4 modified positions → insufficient nuclease resistance | Czauderna et al. 2003 |
| 2'-mod density < 40% | +2 | 4–8 modified positions → partial protection | Czauderna et al. 2003 |

**New July 2026 — Alnylam AT3 clinical PS distribution pattern (Sakamuri et al. 2020):**

A systematic stereopure study of all 64 PS isomers in the Alnylam AT3 siRNA-GalNAc design revealed that the optimal clinical PS pattern is:
- **4 PS on antisense**: positions G1, G2 (5' end), G21, G22 (3' end)
- **2 PS on sense**: positions P1, P2 (5' end)

Wrong stereo-isomers at G21 cause near-zero efficacy. Racemic standard synthesis contains ~50% wrong isomers silently. Correct positional placement is the minimum requirement before stereochemistry matters.

**Implementation**: Counts PS symbols ('S') and 2'-modified positions (M, F, L, E, D, Y, 8, 9, 6, B, J, V, N, O, P, R, H, K, Z, Q, W, X, 7) across both strands. Additionally validates the positional distribution against the Alnylam clinical pattern.

---

## 2. Immuno Penalty (Range: 0–28)

**Biological target**: Innate immune activation via TLR7 (GU-rich ssRNA) and TLR8 (AU-rich ssRNA) sensors. Unmodified uridine triggers interferon and pro-inflammatory cytokine release.

| Condition | Penalty | Rationale | Citation |
|-----------|---------|-----------|----------|
| Unmodified U in AS seed (pos 2–8), each | **+2.0** | Seed uridines are the strongest TLR7/8 trigger | Sioud & Sørensen 2004 |
| Unmodified U in AS tail (pos 9–21), each | **+0.5** | Tail uridines have weaker TLR effect | Goodchild et al. 2009 |
| Unmodified U in sense strand, each | **+1.0** | Sense uridines also stimulate immune sensors | Judge et al. 2005 |
| GU-rich motif: GUUGU | +3 | Most immunostimulatory GU-rich pentamer | Goodchild et al. 2009 |
| GU-rich motif: GUGU | +3 | Moderately immunostimulatory tetramer | Goodchild et al. 2009 |
| GU-rich motif: UGU | +3 | Weakest immunostimulatory trimer | Goodchild et al. 2009 |
| Over-methylation (M count > 24) | +4 | Advisory — extremely high 2'-OMe may cause steric issues or hepatotoxicity | Alnylam ESC design |

### Important Calibrations (C-DAC Panel Review, June 2026)

| Parameter | Before | After | Reason |
|-----------|--------|-------|--------|
| Seed U penalty | +4.0 | **+2.0** | Was overtuning on low-GC sequences; CDAC review suggested halving |
| Tail U penalty | +1.0 | **+0.5** | Consistent with seed U halving |
| Over-methylation threshold | >16 M | **>24 M** | Clinical ESC designs use 25–27 M's safely (Alnylam) |

### Non-Stacking Motif Detection

Motif detection uses a **hierarchical non-stacking** search within each sliding window:

1. Scan each 5-nt window for **GUUGU** — if matched, mask all 5 positions with sentinel `X` so same window cannot also trigger GUGU/UGU.
2. Scan each 4-nt window for **GUGU** — skip if any position is masked.
3. Scan each 3-nt window for **UGU** — skip if any position is masked.

This prevents triple-counting a single GUUGU epitope (GUUGU contains both GUGU and UGU internally). Compare:
- **Before** (stacking): GUUGU at positions 5-9 → +3 (GUUGU) + +3 (GUGU at 5-8) + +3 (UGU at 5-7) = **+9** for one epitope.
- **After** (non-stacking): GUUGU at positions 5-9 → **+3** max (only GUUGU applies, positions 5-9 masked).

---

## 3. RISC Loading Penalty (Range: −10 to 60)

**Biological target**: Guide strand loading into Argonaute-2, passenger strand ejection, and thermodynamic asymmetry. Can be NEGATIVE for beneficial chemistries (bonus).

| Condition | Penalty | Rationale | Citation |
|-----------|---------|-----------|----------|
| Missing 5'-phosphate (AS pos 1) | +5 | 5'-P required for Ago2 MID domain binding | Frank et al. 2010 |
| PS at AS pos 1 | +2 | PS backbone slows loading kinetics | |
| Unmodified seed position (AS 2–8), each | +2 | Seed region plasticity → off-target binding | Jackson et al. 2006 |
| UNA at AS pos 7 (exempt from seed penalty) | **0** | UNA reduces seed stability → no penalty | Bramsen et al. 2010 |
| LNA at AS pos 2, 3, or 4, each | +5 | LNA too rigid for seed region → RISC loading impaired | Hidayah et al. 2021 |
| MOE at AS pos 2–14, each | +3 | MOE bulky 2' modification → some RISC impairment | Prakash et al. 2005 |
| GNA at AS pos 2–5, each | +4 | Glycol nucleic acid small backbone → disrupts seed | Schlegel et al. 2022 |
| GNA at AS pos 6–8, each | **−2 bonus** | GNA at positions 6–8 improves selectivity window | Schlegel et al. 2022 (ESC+) |
| ENA at AS pos 2–8, each | +4 | Ethylene-bridged nucleic acid → significant RISC impact | Morihiro et al. 2020 |
| ENA at AS pos 9–14, each | +2 | ENA at tail positions → moderate impact | Morihiko et al. 2020 |
| TNA at AS pos 2–6, each | +3 | Threose nucleic acid → moderate seed disruption | |
| TNA at AS pos 7 (exempt) | 0 | Position 7 is flexible | |
| TNA at AS pos 8–14, each | +1 | TNA at tail → minimal impact | |
| 2'-F deficiency on pyrimidines < 20% | +6 | 2'-F critical for nuclease resistance + RISC affinity | Layzer et al. 2004 |
| 2'-F deficiency on pyrimidines < 40% | +3 | Partial 2'-F → moderate penalty | Layzer et al. 2004 |
| Exotic mod micro-penalty (Benzyl, Inosine), each | +2 | Large aromatic base modifications may distort helix | |
| Other exotic mods (V, I, N, O, P, R, H, K, Z, Q, W, X, 7), each | +1 | Slight penalty for less-characterized chemistries | |

### RISC Range Expansion (June 2026)

The RISC penalty range was expanded from 31 → 50 → 60 to accommodate the full GNA/ENA/TNA position-split rules and exotic micro-penalties. The floor was lowered to −10 to accommodate GNA@6-8 bonuses stacking on top of 5'-P and PS@1 credits.

---

## 4. Thermo Penalty (Range: 0–20)

**Biological target**: Melting temperature (Tm) of the siRNA duplex — extremes reduce RISC loading specificity and can cause off-target silencing.

| Condition | Penalty | Rationale | Citation |
|-----------|---------|-----------|----------|
| GC < 30% or > 55% | +8 | Extreme GC → Tm too low (off-target) or too high (strands don't separate) | Reynolds et al. 2004 |
| GC 30–35% or 50–55% | +3 | Borderline GC → moderate Tm concern | Reynolds et al. 2004 |
| Palindrome (≥8 nt self-complementary) | +5 | Self-complementarity → hairpins reduce effective siRNA concentration | |
| Homopolymer run (≥4 consecutive same nt) | +5 | Poly-A/G/C/U runs → secondary structure, reduced loading | |
| GC run (≥6 consecutive G or C) | +3 | GC-rich runs → stable secondary structures impede RISC | |

---

## 5. Serum Penalty (Range: −5 to 60)

**Biological target**: Exonuclease degradation in serum/bloodstream (serum nucleases digest unprotected 3' and 5' termini).

**Orthogonality note**: Does NOT check modification density (that is the nuclease penalty's domain). Checks only whether the 3' and 5' termini of both strands are protected by PS backbone linkages, 5'-PO₄, or GalNAc conjugate.

| Condition | Penalty | Rationale | Citation |
|-----------|---------|-----------|----------|
| **GalNAc ('4') at AS 5' end** | **+40 (FATAL)** | **Proven completely inactive in all hepatocytes experiments regardless of valency** | **Weingärtner et al. 2020** |
| AS 5' not PS or '1' (5'-PO₄) | +4 | 5' end unprotected → 5'→3' exonuclease activity | |
| AS 3' not PS | +3 | 3' end unprotected → 3'→5' exonuclease activity | Elmén et al. 2005 |
| SS 5' not PS or '4' (GalNAc) | +3 | 5' end of passenger unprotected | |
| SS 3' not PS or '4' (GalNAc) | +2 | 3' end of passenger unprotected | |
| **Single GalNAc only at sense 5'** | **+3** | **Single GalNAc — significantly reduced in vivo potency vs. 2+ units** | **Weingärtner et al. 2020** |
| **Dual-terminal Sense GalNAc (5' AND 3')** | **−5 (bonus)** | **Novel bi-terminal design: 3–4× superior potency at lower dose (0.3 mg/kg)** | **Weingärtner et al. 2020** |

### New July 2026 — GalNAc Position Rules (Weingärtner et al. 2020, Silence Therapeutics)

A systematic in vitro/in vivo study of 25 GalNAc-siRNA conjugate designs targeting murine transthyretin (Ttr) established these hard experimental rules:

**Rule 1 — Antisense 5' GalNAc is a hard failure:**
GalNAc conjugated to the **5' end of the antisense strand** is experimentally inactive in primary mouse hepatocytes at all tested valencies (1–4 GalNAc units). This makes biological sense — the 5' antisense end must be free for RISC/Ago2 loading. Blocking it with a bulky GalNAc (~600 Da) is functionally equivalent to destroying RISC loading. A +40 penalty effectively hard-rejects this design.

**Rule 2 — Minimum valency for activity:**
A **single GalNAc unit** at sense or antisense 3'/5' positions has low or no activity. Two or more serial GalNAc units are required for robust ASGP-R-mediated endocytosis. A +3 penalty for single GalNAc reflects this reduced potency.

**Rule 3 — Dual-terminal sense GalNAc is the optimal design:**
Two single GalNAc units at **both the 5' and 3' ends of the sense strand** (siRNA-24 design) showed:
- 3× reduction of serum TTR at day 7
- 4× lower serum level at day 27
- Superior lysosomal stability (72h tritosome survival)
...compared to equimolar doses of the triantennary GalNAc positive control.
This is rewarded with a −5 bonus that increases the adjusted efficacy score.

---

## Literature Citations

1. Braasch DA, Corey DR. *Biochemistry* 2004;41(14):4503-4510. PMID: 11926837
2. Czauderna F, et al. *Nucleic Acids Res* 2003;31(11):2705-2716. PMID: 12771180
3. Sioud M, Sørensen DR. *Oligonucleotides* 2004;14(1):1-11. PMID: 15346694
4. Goodchild A, et al. *Oligonucleotides* 2009;19(2):89-98. PMID: 19445602
5. Judge AD, et al. *Nat Biotechnol* 2005;23(4):457-462. PMID: 15778709
6. Frank F, et al. *Nature* 2010;465(7299):818-822. PMID: 20505670
7. Jackson AL, et al. *Nat Biotechnol* 2006;24(9):1151-1157. PMID: 16964229
8. Bramsen JB, et al. *Nucleic Acids Res* 2010;38(9):2861-2878. PMID: 15420340
9. Hidayah NN, et al. *Biochemistry* 2021;60(15):1170-1183. PMID: 33720707
10. Prakash TP, et al. *J Med Chem* 2005;48(21):6696-6705. PMID: 16147155
11. Schlegel MK, et al. *Nucleic Acids Res* 2022;50(16):9056-9070. PMID: 35996904
12. Morihiro K, et al. *J Am Chem Soc* 2020;142(4):1968-1975. PMID: 31867957
13. Layzer JM, et al. *RNA* 2004;10(5):766-771. PMID: 15100431
14. Reynolds A, et al. *Nat Biotechnol* 2004;22(3):326-330. PMID: 14758366
15. Elmén J, et al. *Nucleic Acids Res* 2005;33(1):439-447. PMID: 15653644
16. Bramsen JB, Kjems J. *Front Genet* 2012;3:154. PMID: 22934100
17. Deleavey GF, Damha MJ. *Chem Biol* 2012;19(8):937-954. PMID: 22921066
18. Khvorova A, Watts JK. *Nat Biotechnol* 2017;35(3):238-248. PMID: 28244990
19. Nair JK, et al. *J Am Chem Soc* 2014;136(49):16958-16961. PMID: 25434769
20. Janas MM, et al. *Nucleic Acids Res* 2018;46(15):7679-7692. PMID: 29939332
21. Foster DJ, et al. *Mol Ther* 2018;26(3):708-720. PMID: 29396262
22. Schlegel MK, et al. *Nucleic Acids Res* 2022;50(16):9056-9070. (GNA position-split, ESC+ design)
23. Dar SA, et al. *RNA Biol* 2016;13(8):700-712. PMID: 27348347 (HelixZero-CMS training data)
24. Allerson CR, et al. *J Med Chem* 2005;48(4):901-904. PMID: 15715460 (2'-OMe + PS synergy)
25. Dowler T, et al. *Nucleic Acids Res* 2006;34(6):1669-1675. PMID: 16556911 (ENA thermodynamics)
26. Yang X, et al. *Nucleic Acids Res* 2012;40(8):3393-3404. PMID: 22210852 (PS/nuclease stability)
27. Soutschek J, et al. *Nature* 2004;432(7014):173-178. PMID: 15538359 (GalNAc/PS protection)
28. Kenski DM, et al. *Mol Ther Nucleic Acids* 2012;1:e5. PMID: 23344720 (siRNA localization)
29. **Weingärtner A, et al. *Mol Ther Nucleic Acids* 2020;21:242-254. https://doi.org/10.1016/j.omtn.2020.05.026 (GalNAc position rules — serum penalty, Rules 1–3)**
30. **Sakamuri S, et al. *ChemBioChem* 2020;21:1304-1308. https://doi.org/10.1002/cbic.201900630 (PS stereopure — nuclease penalty, positional distribution)**
