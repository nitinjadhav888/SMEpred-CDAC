# HelixZero-CMS

**Chemical Modification Space Prediction for siRNA Therapeutics**

HelixZero-CMS predicts the silencing efficacy of chemically modified siRNA candidates. Given a 21-nt siRNA duplex, it scores 1,302 single-modification variants (31 symbols × 21 positions × 2 strands) and performs beam search across the multi-modification space (up to 14 simultaneous modifications) — all in ~20 seconds per sequence.

Live Web UI at `http://localhost:8000` (see [Installation](#installation)).

---

## What It Solves

| Problem | Magnitude |
|---------|-----------|
| Single-mod enumeration | 31 mods × 21 pos × 2 strands = **1,302** variants per siRNA |
| Multi-mod space (14 max) | Σ C(1,302, k) for k=1..14 ≈ **10⁶⁸** combinatorial candidates |
| Beam search traversal | ~35k candidates evaluated in **~20 seconds** |
| Biophysics integration | **5 orthogonal penalty domains** (28+ literature citations) |
| Clinical validation | **ESC/ESC+ designs** score ≥50, PK bounds satisfied |

---

## Installation

```bash
pip install -e .
uvicorn api.main:app --reload --port 8000
# → http://localhost:8000
```

**Requires**: Python 3.10+, LightGBM 4.x, scikit-learn 1.6+, NumPy, Pandas, FastAPI, Uvicorn.

---

## Model Accuracy

| Model | Features | Rows | PCC | Spearman | Trees | Leaves |
|-------|----------|------|-----|----------|-------|--------|
| **Naked V4** (unmodified screening) | 214-d sequence | 83,535 | 0.55 | — | — | — |
| **HelixZero v4** (Model B) | **1,467-d** position-aware | 83,535 | **0.822** | **0.823** | 1,115 | 127 |

---

## 6-Step Workflow

```
┌─────┐  ┌──────┐  ┌────────┐  ┌──────────┐  ┌───────────┐  ┌───────┐
│Input│→ │Rank  │→ │Select  │→ │Single-Mod│→ │Multi-Mod  │→ │Export │
│Seq  │  │Tab   │  │Best Hit│  │Scan      │  │Beam Search│  │JSON   │
└─────┘  └──────┘  └────────┘  └──────────┘  └───────────┘  └───────┘
         (gene or      1,302       beam width=30
          dsirna       variants    max_mods=14
          mode)                    ~20 seconds
```

---

## Web UI (4 Tabs)

### 1. 🔬 Rank siRNAs
- **Input**: mRNA/gene sequence or FASTA file.
- **Input type**: `Gene/Transcript` (sliding window) or `DsiRNA (27-mer)` (Dicer cleavage).
- **Output**: All 21-mer candidates sorted by Naked Model score, with seed toxicity and functional checks. Click `Multi-Mod` on any row to jump directly to a beam search.

### 2. 🔧 Single-Mod Scan
- **Input**: 21-nt sense + antisense strands.
- **Output**: All 1,302 single-mod variants ranked by Model B efficacy with biophysical penalty breakdown and seed toxicity.
- Shows **dual baselines**: Naked Model (from Rank tab) vs Model B (recalibrated for chemical space).

### 3. 🔄 Multi-Mod (Beam Search)
- **Input**: 21-nt sense + antisense strands. Optional manual modification strings.
- **Output**: Top multi-mod candidates scored by Model B with delta vs parent.
- Uses plateau-based early stopping (stops when best score improves <0.5 over 3 rounds).
- Expandable rows show color-coded sequence heatmaps with penalty details.

### 4. 📖 Modifications
- Complete legend of all 31 modification symbols with names and chemical types. Tap to copy any symbol.

---

## Biophysical Penalties

Five orthogonal domains adjust the raw efficacy score — strictly non-overlapping to prevent double-counting.

```
adjusted = max(0, min(100, raw − 0.70 × total_penalty))
```

| Domain | Range | What It Checks |
|--------|-------|----------------|
| **Nuclease** | 0–16 | PS backbone coverage, 2'-mod density (endo-nuclease). No termini. |
| **Immunogenicity** | 0–28 | Unmodified U in seed (+2 each), tail (+0.5), sense (+1); GU-rich motifs GUUGU/GUGU/UGU (non-stacking); over-methylation (M>24). |
| **RISC Loading** | −10 to 60 | 5'-P, seed mods (UNA@7 exempt), LNA/MOE/GNA/ENA/TNA position rules, GNA@6-8 bonus (−2), 2'-F deficiency, exotic micro-penalties. |
| **Thermo** | 0–20 | GC extremes, palindrome, homopolymer, GC runs. |
| **Serum** | 0–17 | Termini protection (PS, 5'-PO₄, GalNAc). No density checks. |

**Key calibrations** (C-DAC panel review, June 2026):
- Seed U penalty: +4 → **+2.0**
- Tail U penalty: +1 → **+0.5**
- Over-methylation threshold: >16 → **>24**
- Nuclease/serum orthogonality enforced (no cross-module double-counting)
- Motif detection: non-stacking hierarchical (GUUGU→GUGU→UGU)

---

## API Reference

### POST `/rank`

```json
// Request
{ "sequence": "AUGCAUGCAUG...", "top_n": 20, "input_type": "gene" }
// Response
{ "total_candidates": 42, "input_type": "gene", "results": [...] }
```

### POST `/single-mod`

```json
// Request
{ "sense": "GGAAAUAGACACCAAAUCUUA", "antisense": "UAAGAUUUGGUGUCUAUUUCC", "full_scan": false }
// Response
{ "parent_score": 29.88, "naked_baseline": 27.09, "model_b_baseline": 29.88,
  "total_variants": 1302, "model": "B", "results": [...] }
```

### POST `/multi-mod`

```json
// Request
{ "sense": "GGAAAUAGACACCAAAUCUUA", "antisense": "UAAGAUUUGGUGUCUAUUUCC",
  "sense_mods": "F,,M", "sense_positions": "2,5,,10,12", "antisense_mods": "", "antisense_positions": "" }
```

### POST `/multi-mod-scan`

```json
// Request
{ "sense": "...", "antisense": "...", "max_mods": 14, "beam_width": 30, "full_scan": true }
// Response
{ "parent_score": 29.88, "naked_baseline": 27.09, "model_b_baseline": 29.88,
  "total_variants": 352, "results": [...] }
```

### GET `/modifications`

Returns all 31 modification symbols with `{"canonical": [...], "modifications": [...]}`.

---

## CLI Reference

```bash
python cli/run.py rank --seq "AUGCAUG..." --top-n 10
python cli/run.py single-mod --sense "..." --antisense "..." --top-n 20
python cli/run.py multi-mod --sense "..." --antisense "..." --sense-mods "M,F" --sense-pos "2,5"
python cli/run.py multi-mod-scan --sense "..." --antisense "..." --max-mods 14 --beam-width 30
```

---

## Project Structure

```
smepred/
├── app.html                  # Single-file web UI
├── api/main.py               # FastAPI server (9 endpoints)
├── cli/run.py                # Command-line interface
├── src/
│   ├── predictor.py          # Unified prediction (2 models)
│   ├── sirna_generator.py    # Candidate generation (gene + DsiRNA)
│   ├── features.py           # Feature extraction (214-d + 1,467-d)
│   ├── modification_engine.py # 31 mod symbols, beam search
│   ├── biophysics.py         # 5-domain orthgonal penalties
│   ├── filters.py            # Toxicity + functional checks
│   └── parser.py             # Sequence I/O
├── models/                   # LightGBM .pkl files
├── data/                     # Toxicity tables, mod definitions
├── tests/                    # 32 unit tests + clinical benchmark
├── docs/                     # Full architecture + validation docs
└── scripts/                  # Training + paper generation
```

---

## Evaluation & Validation

| Test | Status | Coverage |
|------|--------|----------|
| Pipeline unit tests | **32/32 PASS** | Features, biophysics all 5 domains, modification engine |
| Clinical benchmark | **4/4 PASS** | ESC/ESC+ designs ≥50, PK bounds, GNA@7 −2 delta |
| Dual baseline verification | **Verified** | Both Naked (27.09) and Model B (29.88) reported |
| Biophysics orthogonality | **Verified** | Nuclease ≠ serum, no cross-module double-counting |

---

## Citation

If you use HelixZero-CMS in your research, please cite:

```
Nitin Jadhav. "HelixZero-CMS: Chemical Modification Space Prediction for
siRNA Therapeutics." CDAC-Pune HPC-M&BA, 2026.
```

---

## License

Research use only. Not approved for clinical decision-making.
