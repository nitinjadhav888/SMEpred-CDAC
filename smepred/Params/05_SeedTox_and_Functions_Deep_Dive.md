# Deep Dive: Seed Toxicity and Functional Structural Filters

When you run an siRNA candidate through the HelixZero-CMS ranking pipeline, two critical safety columns are calculated for every single variant: **SEED TOX** and **Functions**. These columns serve as the absolute gatekeepers for clinical viability. 

This document explicitly breaks down the biological science, exact code logic, and clinical impact of both these systems.

---

## 1. Seed Toxicity (SEED TOX)

### The Biological Concept
The **Seed Region** of an siRNA (specifically positions 2 through 7 on the active Guide/Antisense strand) is responsible for target recognition. 

In human biology, natural microRNAs use this exact 6-nucleotide seed to bind to the 3' UTRs (untranslated regions) of hundreds of different genes simultaneously. If an artificial siRNA drug happens to have a "sticky" seed region that accidentally mimics a human microRNA, it will unintentionally silence hundreds of essential, healthy genes. This causes massive, often fatal, cytotoxicity (cell death).

### How It Is Calculated in the Code (`src/filters.py`)
1. **Extraction:** The model isolates exactly positions 2 through 7 (a 6-mer) from the guide strand.
2. **Empirical Lookup:** It cross-references this 6-mer against an enormous, pre-loaded database (`cell_viability.tsv`). This database contains data from high-throughput empirical screening (Janas et al., *Molecular Cell*, 2018), covering over 4,000 specific seed sequences that have been physically tested on human cells to measure exactly how many cells survive.
3. **Scoring & Labeling:**
   - **Safe (Green):** $\ge$ 70% cell viability. 
   - **Caution (Yellow):** 50% - 69% cell viability. The drug might cause side effects.
   - **Toxic (Red):** < 50% cell viability. The drug is heavily lethal to cells and will fail clinical trials.

### The Mitigation Engine (Seed Rescue)
If a sequence has a inherently toxic seed (e.g., cell viability = 12%), it is not immediately discarded if we are using the **Chemical Modification Engine**. 
- **The Science:** Biological literature (Jackson et al., *RNA*, 2006) proved that inserting specific steric blockades—like bulky chemical modifications—into the seed region physically disrupts the drug's ability to act like a microRNA without stopping its ability to cut its primary target.
- **The Code Logic:** The model scans positions 2-7 for rescuing chemistry: **2'-OMe (M), 2'-Fluoro (F), LNA (L), or 2'-MOE (E)**. 
- **The Impact:** If the code detects these modifications inside a toxic seed, it overrides the "Toxic" label and changes it to **Mitigated (Safe)**, displaying a tooltip like `Seed off-target rescue: 2'-OMe @ pos 2`.

---

## 2. Functional Structural Filters (Functions)

Even if a drug is completely non-toxic and perfectly designed, it can still fail if it structurally knots itself up. The **Functions** column ensures the RNA molecule physically behaves properly inside the human body.

The code strictly evaluates **both the Sense and Antisense strands** against classical RNAi biophysical rules (Reynolds / Ui-Tei design paradigms).

### A. GC Content (30% - 65%)
- **The Concept:** G (Guanine) and C (Cytosine) bind to each other using 3 hydrogen bonds, while A and U use only 2. 
- **The Logic:** If the GC content is > 65%, the two strands of the drug are super-glued together. When the human Argonaute (Ago2) machine tries to unzip them, the engine stalls and the drug fails. If GC is < 30%, the strands fall apart randomly in the bloodstream.
- **Warning Output:** `GC XX% out of optimal 30-65% range`.

### B. Internal Palindromes (Hairpins)
- **The Concept:** RNA is a sticky molecule. A palindrome means the sequence reads the same forwards as its complement reads backwards (e.g., `ACGT...ACGT`). 
- **The Logic:** If a sequence is palindromic, it will bend and stick to *itself*, forming a "hairpin" structure instead of binding to the disease gene. The drug becomes a useless knot.
- **Warning Output:** `✗ Internal palindrome detected (forms stable hairpins)`.

### C. 5-Base Homopolymer Runs
- **The Concept:** A homopolymer run is a sequence of the exact same letter repeating 5 or more times (e.g., `AAAAA` or `GGGGG`).
- **The Logic:** Stretches of repeating nucleotides cause the RNA polymerase and cellular unzipping enzymes to "slip" or get stuck. 
- **Warning Output:** `5-base homopolymer run detected (prevents unwinding)`.

### D. 6-Base Contiguous GC Runs
- **The Concept:** Six consecutive Gs and Cs (e.g., `GGCGCG`).
- **The Logic:** Even if the overall GC percentage is perfect (e.g., 50%), having all the Gs and Cs clumped tightly in one spot creates a microscopic "weld" that the Ago2 enzyme cannot break through.
- **Warning Output:** `6-base contiguous GC run detected`.

## Summary for the Panel
The **Functions** and **SEED TOX** columns represent the absolute minimum biological viability filters. Before the AI even attempts to score how effectively the drug silences a disease, these filters prove that the drug won't knot itself into a useless hairpin, won't get jammed inside the human cell's machinery, and won't inadvertently poison the patient's healthy cells.
