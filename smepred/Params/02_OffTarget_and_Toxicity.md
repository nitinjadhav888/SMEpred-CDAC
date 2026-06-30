# Off-Target Safety: The "Modification-Aware" Engine

A common point of confusion is whether the Off-Target Safety Score is calculated on the "naked" sequence or the "modified" sequence. 

**The answer is BOTH.** The Off-Target engine is explicitly **Modification-Aware**. It calculates the severe biological risks of the bare genetic code, and then calculates exactly how much of that risk you successfully neutralized by applying chemical armor. 

If you run a naked sequence, it might score `40% (Toxic)`. But if you take that exact same sequence and apply the correct chemical modifications, the score will jump to `95% (Cleared)` because the chemistry mitigated the biological risks.

Here is the exact step-by-step breakdown of how the engine calculates the Off-Target score (starting from 100%):

### 1. Thermodynamic Asymmetry (Penalty: up to -40%)
* **The Check:** The drug has two strands (Sense and Antisense/Guide). Ago2 (the cellular machinery) will pick whichever strand has the *weaker* 5' end and throw the other one away. We strictly want it to pick the Antisense strand.
* **The Math:** The model calculates the exact thermodynamic energy of both ends. If the Sense strand is weaker, the cell will accidentally load the wrong strand, turning your entire drug into a massive off-target bomb.
* **The Penalty:** `-40%` if it favors the wrong strand, and `-5%` if the Antisense strand doesn't start with an optimal `A` or `U`.

### 2. The 15-mer Slicer Check (Penalty: Instant 0%)
* **The Check:** The model scans a massive database containing the entire human transcriptome (every known human gene). It checks if a 15-nucleotide chunk of your drug perfectly matches any unintended gene.
* **The Penalty:** If a 15-mer matches an unintended gene, the drug will act like molecular scissors and slice that healthy gene in half. This is catastrophic, so the model drops the score instantly to **0% (TOXIC)**.

### 3. Seed Region Match & Mitigation (Penalty: up to -30%)
* **The Check:** The model takes your 6-mer seed region (positions 2-7) and counts exactly how many times it appears across the entire human genome.
* **The Modification Awareness (The Rescue):** 
  * If the naked seed matches 5 times, the engine applies a massive penalty (up to `-30%`) because it will cause miRNA-like toxicity.
  * **However**, the engine checks your applied chemistry. If it sees you placed a `2'-OMe` at position 2, or a `GNA` at position 7, it knows you have chemically blocked that toxicity. The penalty is drastically reduced from `-30%` down to a maximum of `-5%`. 

### 4. TLR (Toll-Like Receptor) Immune Masking (Penalty: -15% per motif)
* **The Check:** The human immune system searches for specific viral RNA patterns (like `GUUGU` or `UGGC`). If your drug contains these, it triggers a cytokine storm (severe inflammation).
* **The Modification Awareness:** The model checks if those specific dangerous letters have a `2'-OMe` modification placed exactly on top of them. If they are modified, the immune system is "blind" to them. If they are left naked, the model deducts **-15%** for every single unmasked alarm trigger it finds.

### 5. Delivery Conjugation (Penalty: -10%)
* **The Check:** If the drug is designed to target the liver (hepatic delivery), it physically cannot enter liver cells without a specific "key" attached to it called `GalNAc`.
* **The Modification Awareness:** The engine checks your modification string. If it does not find the GalNAc conjugate (modification symbol `4`), it deducts **-10%** because the drug will just float uselessly in the bloodstream.
