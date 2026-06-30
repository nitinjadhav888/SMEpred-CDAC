# Biophysical Penalties Explained

When the model evaluates a chemically modified drug, it outputs a list of penalties (e.g., `Nuclease: -4.0`, `🛡️ Immuno: -5.5`, etc.). These penalties deduct points from the AI's raw score to calculate the true clinical efficacy.

## 1. 🩸 Nuclease Penalty
- **What it is:** Nucleases are enzymes in the human blood that act like molecular scissors. Their job is to chop up and destroy free-floating RNA (like viruses).
- **Why we penalize:** If our drug is injected without armor, it will be destroyed in seconds before it ever reaches the target organ. 
- **How to fix it:** The model requires the "Termini" (the exposed ends of the drug, specifically the 5' and 3' ends) to be armored with Phosphorothioate (PS) linkages. If these are missing, the drug receives a massive Nuclease penalty.

## 2. 🛡️ Immuno Penalty
- **What it is:** The human immune system has sensors (like TLR7 and TLR8) that detect foreign RNA. 
- **Why we penalize:** Unmodified Uridines (U bases) and specific motifs (like GUUGU) trigger the immune system. If injected, the patient's body will launch a massive inflammatory response (cytokine storm), which can be fatal.
- **How to fix it:** The model penalizes high numbers of unmodified Uridines. To clear this penalty, the sequence must be heavily modified with 2'-OMe or 2'-F chemistries, effectively cloaking the drug from the immune system.

## 3. 🎯 RISC Penalty (Argonaute Accommodation)
- **What it is:** RISC (RNA-induced silencing complex) and its core protein, Argonaute 2 (Ago2), is the microscopic machine inside the cell that actually uses our drug to destroy the disease gene.
- **Why we penalize:** Ago2 is extremely picky about the physical shape of the drug. If we put bulky chemical modifications (like LNA, ENA, or MOE) in the wrong places, the drug physically won't fit into the machine. It's like putting a square peg in a round hole.
- **How we calculate it:** The model strictly enforces positional rules. For example, position 1 must be A or U. Bulky modifications in the central catalytic cleft (positions 10-14) are heavily penalized because they jam the machine's scissors.

## 4. 🌡️ Thermo Penalty
- **What it is:** Thermodynamics refers to the binding energy required to unzip the two strands of the drug.
- **Why we penalize:** The drug is administered as a double helix (two strands bound together). Ago2 must rip them apart and throw away the "Passenger" strand to use the "Guide" strand. If the 5' end of the Guide strand is bound too tightly, Ago2 will accidentally throw away the Guide strand and keep the Passenger strand. The drug will fail entirely.
- **How we calculate it:** Based on nearest-neighbor thermodynamics (SantaLucia 2004), the model calculates the exact energy at both ends of the drug. It enforces strict "Thermodynamic Asymmetry"—the correct end must be weaker so it opens first.

## 5. 💉 Serum Penalty
- **What it is:** Checks for overall stability in blood serum during transit to the target organ.
- **Why we penalize:** Similar to nuclease degradation, but focused on overall backbone stability rather than just the tips. A naked RNA backbone rapidly hydrolyzes in blood.

## Summary
The physics engine ensures that the AI only recommends drugs that are fully armored against blood enzymes (Nuclease/Serum), cloaked from the immune system (Immuno), thermodynamically balanced (Thermo), and perfectly shaped to fit the cellular machinery (RISC).
