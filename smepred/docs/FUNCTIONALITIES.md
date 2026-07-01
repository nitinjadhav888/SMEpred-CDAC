# HelixZero-CMS: Comprehensive Functionalities Documentation

This document explicitly outlines **every single feature, module, and prediction** that HelixZero-CMS is capable of performing. It serves as the master catalog of the platform's computational and biological capabilities.

---

## 1. Sequence Ingestion & Normalization

The system is designed to handle raw biological data and prepare it strictly for RNAi simulation.

*   **Format Flexibility**: Accepts raw sequence strings, inline FASTA formats (`>header\nsequence`), or direct FASTA file uploads.
*   **DNA-to-RNA Normalization**: Automatically detects DNA sequences and converts Thymine (T) to Uracil (U).
*   **Strict Character Validation**: Strips whitespace and explicitly rejects any non-canonical bases (anything outside of A, U, G, C) to ensure the machine learning pipeline is never fed invalid states.

---

## 2. In Silico siRNA Generation

Given a target mRNA sequence, the system can systematically slice it into therapeutic candidate molecules.

*   **Sliding Window Scanning (Gene Mode)**: Scans an entire mRNA sequence using a 21-nucleotide sliding window (step size = 1) to generate all possible 21-mer siRNA candidates (`N - 20` candidates).
*   **Dicer-Substrate Mode (27-mer Mode)**: If the user inputs a long 25–30nt DsiRNA, the system bypasses the sliding window and accurately simulates human Dicer cleavage, extracting exactly the biologically active 21-mer mature siRNA duplex.
*   **Reverse Complementation**: For every sense (passenger) strand generated, the system computes the exact Watson-Crick reverse-complement antisense (guide) strand.

---

## 3. LightGBM Machine Learning Predictions

The core statistical prediction engine utilizes two separate, highly specialized LightGBM models to predict in vitro silencing efficacy.

*   **Naked Model (Initial Screening)**:
    *   **Input Features**: 214-dimensional sequence features (Sense one-hot encoding, Trinucleotide composition, GC%).
    *   **Function**: Rapidly screens thousands of unmodified candidates to find the most naturally potent sequences.
*   **Model B (Chemical Awareness)**:
    *   **Input Features**: 1,467-dimensional position-aware features. It encodes exactly *what* chemical modification is at *which* position across all 42 nucleotides (21 sense + 21 antisense). It also encodes the parent sequence simultaneously to learn the "delta" effect of the modification.
    *   **Function**: Predicts the statistical efficacy of highly modified sequences.
*   **Platt Calibration**: Both models have their regression outputs passed through an isotonic calibrator (Platt scaling) to transform arbitrary values into a standardized `0–100` biological relevance scale.

---

## 4. Multi-Domain Biophysical Penalty System

Because statistical ML models often fail to account for rigid physiological realities (like serum nucleases or immune receptors), HelixZero-CMS applies a **deterministic 5-domain penalty engine** on top of the ML score. 

### Domain 1: Nuclease Resistance (Endonuclease)
Penalizes designs that will be rapidly degraded by blood enzymes.
*   Checks for Phosphorothioate (PS) backbone coverage.
*   Calculates 2'-modification density (e.g., 2'-OMe, 2'-F). Applies massive penalties if <20% of the duplex is protected.

### Domain 2: Immunogenicity & Toxicity (TLR7/8)
Penalizes designs that trigger the innate immune system (Interferon response).
*   Identifies exposed, unmodified Uridines (highly toxic) and penalizes based on position (Seed vs. Tail vs. Sense strand).
*   Performs a non-stacking hierarchical search for toxic GU-rich motifs (`GUUGU`, `GUGU`, `UGU`).
*   Flags over-methylation (clinical ESC standard: warning if >24 2'-OMe modifications are used, as this can stunt activity).

### Domain 3: RISC Loading Stereochemistry (Ago2)
Enforces the structural geometry required by Argonaute 2 (Ago2) to successfully load the guide strand and cleave the target.
*   **5'-Phosphate Requirement**: Enforces the absolute necessity of a 5'-PO4 or analog on the antisense strand.
*   **Locked Nucleic Acid (LNA) Rules**: Massively penalizes LNA at the antisense 5' position or in the catalytic cleft (pos 10, 11, 13) due to structural rigidity, but permits it at 3' overhangs.
*   **Positional Bonuses**: Understands that GNA (Glycol Nucleic Acid) is toxic in the early seed (pos 2–5) but provides a **therapeutic bonus** if placed at positions 6–8 (the ESC+ clinical design).
*   Penalizes Bulky modifications (MOE, ENA) if placed in the central catalytic cleft.

### Domain 4: Thermodynamic Profiling
Penalizes sequences with extreme melting temperatures that reduce specificity.
*   **GC Bounds**: Heavily penalizes GC% < 30% or > 65%.
*   **Structural Flaws**: Detects and penalizes homopolymeric runs (≥4 same nt), GC runs (≥6 consecutive G/C), and Palindromes (hairpin formation risk).

### Domain 5: Serum Stability (Exonuclease)
Ensures the termini (ends) of the RNA duplex are protected from immediate degradation in the bloodstream.
*   Checks that all four termini (Sense 5'/3', Antisense 5'/3') are capped with either PS linkages, 5'-Phosphates, or GalNAc conjugates.
*   **GalNAc Positional Rules**: Issues a fatal penalty if GalNAc is placed on the Antisense 5' (blocks RISC), but rewards a massive bonus for dual-terminal Sense (5' + 3') GalNAc placement.

---

## 5. Seed Toxicity & Off-Target Safety Module

HelixZero-CMS predicts whether the siRNA will accidentally silence vital survival genes (miRNA-like off-target effects).

*   **Hexamer Lookup Engine**: Extracts the seed region (Antisense positions 2–7) and cross-references it against a built-in database of 4,097 experimentally validated seed hexamers.
*   **Cell Viability Prediction**: Returns an exact predicted cell viability percentage.
*   **Risk Categorization**: Labels the sequence as **Safe** (≥75%), **Caution** (55-74%), or **Toxic** (<55%).
*   **Seed Rescue AI**: Detects if the user has placed a specific thermodynamic destabilizer (like 2'-F, 2'-OMe, or LNA) in the seed region. If a rescue modification is present on a previously "Toxic" seed, the label is dynamically upgraded to **Mitigated**.

---

## 6. Combinatorial Modification Engine

The platform provides a vast chemical workspace for in silico drug design.

*   **Chemical Catalog**: Supports 31 distinct chemical modifications (Canonical bases, 2'-Sugars, Backbone linkages, Conjugates, and Exotic probe chemistries).
*   **Single-Mod Scanning**: Automatically generates and scores all 1,302 possible single-modification permutations for any given duplex (31 mods × 21 positions × 2 strands) in seconds.
*   **Multi-Mod Beam Search**: An autonomous AI search algorithm that intelligently combines chemical modifications. 
    *   It initializes a diverse "beam" of starting chemicals.
    *   Iteratively pairs modifications together, preventing steric clashes.
    *   Utilizes a plateau-detection algorithm to stop searching once the biological score stops improving.
*   **Custom Manual Design**: Users can manually inject exact combinations of chemicals (e.g., "M,F,,S") at exact positions and instantly receive the adjusted biological score.

---

## 7. Interactive UI & Transparency

The frontend is designed to make complex biophysics highly interpretable.

*   **Penalty Breakdowns**: Hovering over any score produces a detailed breakdown showing exactly which biological rules (Nuclease, Immuno, RISC, Thermo, Serum) subtracted points from the design.
*   **Visual Heatmaps**: Expanding a row visually highlights the sequence, color-coding the exact locations of chemical modifications.
*   **Dual-Baseline Transparency**: When testing modifications, the UI displays both the Naked Baseline and the recalibrated Chemical Baseline to ensure users understand exactly how much value the chemistry added.
*   **Confidence Badges**: Automatically tags chemicals as "Standard" (FDA approved), "Exotic" (Preclinical), or "Rare" to guide user design choices.
