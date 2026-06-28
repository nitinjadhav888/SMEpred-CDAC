import os
import sys

# Add parent dir to path so we can import src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.predictor import predict_modified

def run_positional_walk():
    # Use a solid validated siRNA parent sequence
    sense = "GACUUCUAGMAAGGCMAAMUU"
    antisense = "UUAUMABCCSUFCFAGAAGESS"
    
    print("\n" + "="*60)
    print("CLINICAL VALIDATION: POSITIONAL WALK (2'-OMe on Antisense)")
    print("="*60)
    print("Testing the impact of a single 2'-OMe ('M') modification walking from position 1 to 21 on the Antisense strand.\n")
    print(f"{'Position':<10} | {'Efficacy Score':<15} | {'Penalty Details'}")
    print("-" * 60)
    
    # We will use the 'multimod' API directly by passing a single mod at each position
    for pos in range(1, 22):
        # We will use the multimod function which accepts lists
        res = predict_modified(
            sense=sense,
            antisense=antisense,
            mode="multimod",
            antisense_mods="M",
            antisense_positions=str(pos)
        )
        
        # Get the top (only) result
        if res["results"]:
            cand = res["results"][0]
            # Print position, score, and penalties
            p = getattr(cand, 'penalties', {})
            pen_str = ", ".join([f"{k}: -{v:.1f}" for k, v in p.items() if v > 0])
            score = round(cand.efficacy_score, 1)
            print(f"{pos:<10} | {score:<15} | {pen_str}")
        else:
            print(f"{pos:<10} | {'FAILED':<15} |")

if __name__ == "__main__":
    run_positional_walk()
