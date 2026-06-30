# Certificate of Biological Safety
**Human RNAi Therapeutics Off-Target Scan**

## Sequence Details
* **Sense Strand:** `GGCAGAGACAAUAAAACAUUC`
* **Antisense Strand:** `GAAUGUUUUAUUGUCUCUGCC`
* **Antisense Modifications:** `None`

## Validation Results
* **Status:** **WARNING_SEED**
* **Overall Safety Score:** **40.0%**

### Risk Factors Identified
* WARNING: Unmasked TLR7/8 motif (UGU) found in Antisense strand. High risk of innate immune activation (Interferon response).
* WARNING: Unmasked TLR7/8 motif (UGU) found in Antisense strand. High risk of innate immune activation (Interferon response).
* WARNING: Unmasked TLR7/8 motif (UUG) found in Antisense strand. High risk of innate immune activation (Interferon response).
* WARNING: Missing GalNAc ('4') delivery conjugate. Predicted hepatic uptake is 0%. In vivo pharmacokinetic profile will fail.

### Safety Notes & Mitigations
* Note: Antisense 5' end is not A or U. This is sub-optimal for Ago2 MID-domain anchoring.

---
*Validated checks include: Thermodynamic Asymmetry end-loading preference, Slicer-mediated 15-mer exclusion, and Seed-region mismatch mitigation against the GRCh38 Human Transcriptome.*