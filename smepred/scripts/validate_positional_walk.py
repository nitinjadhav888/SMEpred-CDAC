import sys
from pathlib import Path
import json

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.features import extract_positional_features_batch
from src.predictor import _get_model
from src.biophysics import calculate_adjusted_efficacy

def run_walk():
    # Standard naked sequence
    sense = "GGAAAUAGACACCAAAUCUUA"
    antisense = "UAAGAUUUGGUGUCUAUUUCC"

    model_b = _get_model("B")
    results_out = []

    # Positional Walk on Antisense Strand (Positions 1 through 21)
    # We will place a single 2'-OMe ('M') at each position while keeping the rest unmodified
    for pos in range(1, 22):
        # Create unmodified strings
        anti_mods = list("U" * 21) # Placeholder for unmodified, wait, the mod engine uses "" or standard chars
        # Actually, if it's unmodified we just use a blank string, but the feature extractor 
        # expects the actual nucleotide if it's unmodified, or we just pass the naked sequence?
        # In our features, if it's unmodified we pass the nucleotide. 
        # But wait, the exact syntax for un-modified in our extractor is just the nucleotide itself, 
        # OR we just use the `predict_modified` function which takes care of this!
        pass

    # Let's use the exact feature generation logic
    from src.modification_engine import multimod_gen

    for pos in range(1, 22):
        # Create a single modification at 'pos'
        res = multimod_gen(
            sense, antisense, 
            sense_mods="", sense_positions="", 
            antisense_mods="M", antisense_positions=str(pos)
        )
        
        # res.mod_antisense now has the modification
        X = extract_positional_features_batch([res.sense], [res.antisense], [sense], [antisense])
        raw_score = float(model_b.predict(X)[0])
        adj_score, penalties, total_pen = calculate_adjusted_efficacy(raw_score, res.sense, res.antisense, sense, antisense)

        results_out.append({
            "Position": pos,
            "Adjusted_Score": round(adj_score, 1),
            "Raw_Score": round(raw_score, 1),
            "RISC_Penalty": penalties["risc"]
        })

    print(json.dumps(results_out, indent=2))

if __name__ == "__main__":
    run_walk()
