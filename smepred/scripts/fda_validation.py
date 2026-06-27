import sys
from pathlib import Path
import json

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.features import extract_positional_features_batch
from src.predictor import _get_model
from src.biophysics import adjusted_efficacy_score

def run_validation():
    sense = "GGAAAUAGACACCAAAUCUUA"
    antisense = "UAAGAUUUGGUGUCUAUUUCC"

    tests = [
        {
            "name": "Positive Control: ESC Chemistry (FDA-Approved Pattern)",
            "desc": "Alternating 2'-F/2'-OMe with terminal Phosphorothioates and 5'-Phosphate.",
            "sense_mods": "SSMMMMMMMMMMMMMMMMMM4",
            "anti_mods": "1SMMMFFFMMFMFFFMFFFSS"
        },
        {
            "name": "Positive Control: ESC+ Chemistry (GNA@7)",
            "desc": "Same as ESC, but introduces Glycol Nucleic Acid (GNA) at position 7 to rescue seed toxicity.",
            "sense_mods": "SSMMMMMMMMMMMMMMMMMM4",
            "anti_mods": "1SMMMF8FMMFMFFFMFFFSS"
        },
        {
            "name": "Negative Control: Steric Blockade (LNA in Seed)",
            "desc": "Places bulky LNA modifications in the seed region (pos 2-5) which physically blocks Ago2 loading.",
            "sense_mods": "SSMMMMMMMMMMMMMMMMMM4",
            "anti_mods": "1SLLLLFFMMFMFFFMFFFSS"
        },
        {
            "name": "Negative Control: Unprotected Termini",
            "desc": "Removes terminal Phosphorothioates and conjugates. Will be rapidly degraded by blood exonucleases.",
            "sense_mods": "MMMMMMMMMMMMMMMMMMMMM",
            "anti_mods": "MMMMMMMMMMMMMMMMMMMMM"
        },
        {
            "name": "Negative Control: Over-Methylated",
            "desc": "Completely methylated without alternating Fluoro. High immune penalty.",
            "sense_mods": "MMMMMMMMMMMMMMMMMMMMM",
            "anti_mods": "1MMMMMMMMMMMMMMMMMMMM"
        }
    ]

    model_b = _get_model("B")
    results_out = []

    for t in tests:
        mod_sense = t["sense_mods"]
        mod_anti = t["anti_mods"]

        X = extract_positional_features_batch([mod_sense], [mod_anti], [sense], [antisense])
        raw_score = float(model_b.predict(X)[0])
        adj_score, penalties, total_pen = adjusted_efficacy_score(raw_score, mod_sense, mod_anti, sense, antisense)

        results_out.append({
            "Test": t["name"],
            "Description": t["desc"],
            "Raw_Score": round(raw_score, 1),
            "Adjusted_Score": round(adj_score, 1),
            "Total_Penalty": round(total_pen, 1),
            "Penalties": penalties,
            "Label": "High" if adj_score >= 60 else "Low"
        })

    print(json.dumps(results_out, indent=2))

if __name__ == "__main__":
    run_validation()
