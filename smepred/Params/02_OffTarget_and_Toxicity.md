# Off-Target Safety and Structural Toxicity

This document explains the safety checks the model performs to ensure the drug doesn't accidentally harm the patient.

## 1. Off-Target Safety Scan
- **What it is:** The system scans the entire human transcriptome (all RNA in the human body) to see if our siRNA drug accidentally binds to the wrong gene.
- **Why it matters:** If an siRNA designed to cure a liver disease accidentally matches a gene responsible for heart function, it will silence the heart gene, causing severe side effects.
- **How we calculate it:** We align the sequence against the human transcriptome database. A "Cleared (90%)" status means the sequence is highly unique and will not mistakenly destroy essential non-target genes.
- **The Science:** We use strict thermodynamic rules (SantaLucia 2004) and ensure that the binding energy strictly discriminates between the target and random genes. We also specifically test the *Guide Strand* (the active half of the drug) since that is what the cell uses to hunt for targets.

## 2. Seed Toxicity
- **What it is:** The "Seed Region" consists of positions 2 through 8 on the siRNA drug. This region is hyper-sensitive.
- **Why it matters:** In human biology, microRNAs (natural regulators) use positions 2-8 to bind to hundreds of different genes at once. If our drug's seed region is too "sticky", it will act like a rogue microRNA and shut down hundreds of healthy genes, causing massive cellular toxicity.
- **The Fix:** The model checks if we have applied specific chemical modifications (like GNA at position 7) to intentionally weaken the seed region. If properly modified, the Seed Toxicity is neutralized.

## 3. Internal Palindromes & Hairpins (Structural Warnings)
- **What it is:** A warning (e.g., `✗ Internal palindrome detected`) means the sequence is symmetrical and can fold back and bind to itself.
- **Why it matters:** RNA is sticky. If the left side of the sequence perfectly matches the right side, the drug will fold into a "hairpin" knot. 
- **The Impact:** When the drug knots up, the human cellular machinery (RISC) cannot open it. The drug becomes entirely useless. The model automatically flags and rejects these sequences.

## 4. GC Content (Thermodynamic Balance)
- **What it is:** The percentage of G (Guanine) and C (Cytosine) bases in the sequence. 
- **How we calculate it:** We require a GC content strictly between 30% and 65%.
- **Why it matters:** G-C bonds are very strong (3 hydrogen bonds). If the GC content is too high (>65%), the drug is too tightly bound and the cell can't unzip it to use it. If it's too low (<30%), the drug falls apart before it can work.
