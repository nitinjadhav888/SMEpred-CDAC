# Limitations, Improvements, and Validation Notes

## Resolved Issues

| Issue | Resolution | Status |
|-------|-----------|--------|
| Documentation lag | All docs updated to match current codebase. Dual baselines (Naked/Model B), DsiRNA mode, recalibrated penalties, non-stacking motif detection, and orthogonality fixes documented. | ✅ |
| Dataset-dependent quality | Model B achieves PCC=0.822 (83,535 rows). Clinical benchmark validates ESC/ESC+ designs (4/4 pass). | ✅ |
| Off-target/toxicity | Janas (2018) 6-mer seed lookup + rescue mod detection. Seed U penalty recalibrated from +4 → +2.0 to match clinical data. | ✅ |
| Multi-mod search space | Beam search with pairing pool cap (3× beam width) completes in ~20s. Plateau-based early stopping prevents unnecessary rounds. | ✅ |
| Over-engineering bias | No artificial over-mod penalty. Beam search finds natural optimal mod count. Exotic micro-penalties are subtle (+1–2). | ✅ |
| PDF page numbering | Two-pass approach (count → build with "X of Y") fixes reportlab duplicate page bug. | ✅ |
| Score breakdown overflow | JS-driven viewport-aware popup positioning (right → left → below, clamped). | ✅ |
| Cross-module double-counting | Nuclease = PS count + density only. Serum = termini only. No overlapping checks. | ✅ |
| Immuno over-penalty | Non-stacking motif search prevents triple-counting GUUGU. Seed U penalty halved (4→2). | ✅ |

## Remaining Improvements Needed

| Priority | Improvement | Description | Status |
|----------|-------------|-------------|--------|
| P0 | Transcriptome off-target | Full Bowtie alignment against human transcriptome to detect full-sequence homology with off-target genes. Current seed-based check only catches miRNA-like toxicity. | ⏳ Phase 2 |
| P1 | Systems toxicology | TOXsiRNA-style prediction integrating chemical structure, cellular pathways, and clinical safety data. | ⏳ Phase 2 |
| P1 | Co-modification encoding | Current 1-symbol-per-position scheme cannot represent PS backbone + 2'-O-mod at the same nucleotide (e.g., 2'-OMe + PS). Would require a multi-hot encoding or paired symbols. | ⏳ Phase 2 |
| P2 | Delivery integration | LNP/GalNAc/Delivery vehicle simulation — currently assumes GalNAc conjugate (symbol '4') is applied. No liposomal formulation modeling. | ⏳ Phase 3 |
| P3 | Prospective validation | Wet-lab validation of top-ranked designs against actual RISC loading and silencing data. | ⏳ External |
| P3 | Model versioning | Automated experiment tracking for training runs (hyperparameters, data splits, metrics). | ⏳ Phase 4 |

## Validation Against Published Tools

| Tool | Metric | Our Value | Published | Notes |
|------|--------|-----------|-----------|-------|
| si-Fi | PCC | **0.822** | 0.719 | +0.103 improvement |
| OligoFormer | Spearman | **0.823** | 0.797 | +0.026 improvement |
| TOXsiRNA | Seed toxicity | 4,097 seeds | ~2,500 seeds | 1.6× larger lookup table |
| Mandelli SVR | Feature dim | **1,467** | 214 | 6.8× richer representation |
| CMS Therapeutics | ESC validation | **62.0 adj** | Clinical | Score aligns with clinically approved range |

## Encoding Limitation Detail

The current single-symbol-per-position encoding (`A`, `U`, `G`, `C`, `F`, `M`, `S`, `L`, `E`, `D`, `1`, `2`, `3`, `5`, `6`, `8`, `9`, `Y`, `4`, `B`, `J`, `V`, `I`, `N`, `O`, `P`, `R`, `H`, `K`, `Z`, `Q`, `W`, `X`, `7`) cannot represent a nucleotide with **both** a backbone modification (PS) and a sugar modification (2'-OMe) simultaneously. In reality, PS + 2'-OMe at the same position is a common clinical combination. The model sees either `S` (PS-only) or `M` (2'-OMe-only) but not both.

This is a known area for improvement in a future multi-hot encoding layer.
