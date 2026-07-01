"""
cli/run.py — Command-line interface for HelixZero-CMS.

Three sub-commands:
  rank        : rank siRNA candidates by modification potential (default)
  single-mod  : scan all 1260 single-modification variants for one siRNA
  multi-mod   : predict efficacy of a custom multi-modification design

Usage examples:

  # Rank siRNAs from a FASTA file (show top 10 by modification potential)
  python cli/run.py rank --input gene.fasta --top 10

  # Rank from inline sequence (default: modification potential)
  python cli/run.py rank --sequence AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUG

  # Rank by naked score instead
  python cli/run.py rank --sequence AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUG --naked

  # Single-mod scan on a chosen siRNA (all 1260 variants, model A)
  python cli/run.py single-mod \\
      --sense GCAGCACGACUUCUUCAAGUU \\
      --antisense CUUGAAGAAGUCGUGCUGCUU \\
      --model A --top 20

  # Multi-mod: apply F at positions 2,5 on sense and M at positions 10,12 on antisense
  python cli/run.py multi-mod \\
      --sense GCAGCACGACUUCUUCAAGUU \\
      --antisense CUUGAAGAAGUCGUGCUGCUU \\
      --sense-mods F --sense-positions 2,5 \\
      --antisense-mods M --antisense-positions 10,12 \\
      --model A
"""

import sys
import json
from pathlib import Path

import click
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

from src.predictor import rank_sirnas, rank_by_naked_score, predict_modified


# ─── output helpers ───────────────────────────────────────────────────────────

def _print_table(rows: list, fields: list):
    if not rows:
        click.echo("No results.")
        return
    widths = {f: max(len(f), max(len(str(r.get(f, ""))) for r in rows)) for f in fields}
    header = "  ".join(f.upper().ljust(widths[f]) for f in fields)
    sep    = "  ".join("-" * widths[f] for f in fields)
    click.echo(header)
    click.echo(sep)
    for r in rows:
        click.echo("  ".join(str(r.get(f, "")).ljust(widths[f]) for f in fields))


# ─── rank command ─────────────────────────────────────────────────────────────

@click.group()
def cli():
    """HelixZero-CMS — siRNA efficacy prediction pipeline."""


@cli.command("rank")
@click.option("--input",    "-i", "fasta_path", default=None, help="FASTA file path")
@click.option("--sequence", "-s", default=None,  help="Inline mRNA/gene sequence")
@click.option("--top",      "-n", default=20,    show_default=True, help="Number of top candidates to show")
@click.option("--pool",      default=50,           show_default=True, help="Number of candidates to evaluate with mini-scan")
@click.option("--naked",     is_flag=True,         help="Use naked-model ranking instead of modification potential")
@click.option("--output",    "-o", default=None,  help="Save results to CSV file")
def rank_cmd(fasta_path, sequence, top, pool, naked, output):
    """
    Rank siRNA candidates by *modification potential* (default) or naked score (--naked).

    The naked model (PCC=0.55) correlates poorly with modification potential
    (Spearman ρ≈0.44). By default, this command scores all candidates with a quick
    naked filter, then runs a 40-variant mini-scan (top 4 mod types, antisense
    pos 1-10) on the top POOL candidates, and reranks by best modified score.

    Scores range from 0 (no silencing) to 100 (maximum silencing).
    Very High = ≥80, High = 70-79, Moderate = 55-69, Low = <55.
    """
    if not fasta_path and not sequence:
        raise click.UsageError("Provide --input or --sequence.")
    source = fasta_path or sequence

    click.echo(f"\nGenerating and ranking siRNA candidates ...\n")
    results = rank_by_naked_score(source, top_n=top)

    fields = ["rank", "position", "sense", "antisense", "efficacy_score", "efficacy_label"]
    rows = [r.to_dict() for r in results]
    _print_table(rows, fields)

    if output:
        import pandas as pd
        pd.DataFrame(rows).to_csv(output, index=False)
        click.echo(f"\nSaved to {output}")

    click.echo(f"\nShowing top {len(results)} of all candidates.")


# ─── single-mod command ───────────────────────────────────────────────────────

@cli.command("single-mod")
@click.option("--sense",     required=True, help="21-nt sense strand sequence")
@click.option("--antisense", required=True, help="21-nt antisense strand sequence")
@click.option("--model",     default="B", show_default=True, type=click.Choice(["B"]), help="Model: B (default, unified HelixZero model)")
@click.option("--top",  "-n", default=20, show_default=True, help="Top N results to display")
@click.option("--output",    default=None, help="Save results to CSV file")
def single_mod_cmd(sense, antisense, model, top, output):
    """
    Scan all 1260 single-modification variants of a siRNA and rank by predicted efficacy.

    For each of 30 modification types × 21 positions × 2 strands, predicts efficacy.
    The delta_score column shows improvement (positive) or reduction (negative)
    compared to the unmodified parent siRNA.
    """
    click.echo(f"\nGenerating 1260 cm-siRNA variants and predicting efficacy (Model {model}) ...\n")
    out = predict_modified(sense, antisense, mode="scan", model_key=model)
    results = out["results"][:top]

    fields = ["rank","mod_symbol","mod_strand","mod_position","efficacy_score","delta_score","efficacy_label"]
    rows = [r.to_dict() for r in results]
    _print_table(rows, fields)

    # Show top variant's biophysics penalties
    if results:
        r = results[0]
        p = getattr(r, 'biophysics', None)
        if p:
            click.echo(f"\nBiophysical penalties (subtracted from raw):")
            for k, v in p.items():
                click.echo(f"  {k}: -{v:.1f}")

    if output:
        import pandas as pd
        out_all = predict_modified(sense, antisense, mode="scan", model_key=model)
        pd.DataFrame([r.to_dict() for r in out_all["results"]]).to_csv(output, index=False)
        click.echo(f"\nSaved all results to {output}")


# ─── multi-mod command ────────────────────────────────────────────────────────

@cli.command("multi-mod")
@click.option("--sense",               required=True, help="21-nt sense strand")
@click.option("--antisense",           required=True, help="21-nt antisense strand")
@click.option("--sense-mods",          default="",    help="Mod symbols for sense strand (e.g. F,,M)")
@click.option("--sense-positions",     default="",    help="Positions for sense mods (e.g. 2,5,,10,12)")
@click.option("--antisense-mods",      default="",    help="Mod symbols for antisense strand")
@click.option("--antisense-positions", default="",    help="Positions for antisense mods")
@click.option("--model",               default="B",   show_default=True, type=click.Choice(["B"]))
def multi_mod_cmd(sense, antisense, sense_mods, sense_positions, antisense_mods, antisense_positions, model):
    """
    Predict efficacy of one custom multi-modification cm-siRNA design.

    Separate multiple modification types with ',,' in both --sense-mods and --sense-positions.
    Example: --sense-mods F,,M --sense-positions 2,5,,10,12
      means: apply F at positions 2 and 5, apply M at positions 10 and 12 on the sense strand.
    """
    click.echo(f"\nApplying custom modifications and predicting efficacy (Model {model}) ...\n")
    out = predict_modified(
        sense, antisense,
        mode="multimod",
        model_key=model,
        sense_mods=sense_mods,
        sense_positions=sense_positions,
        antisense_mods=antisense_mods,
        antisense_positions=antisense_positions,
    )
    results = out["results"]
    if results:
        r = results[0]
        click.echo(f"Modified sense    : {r.sense}")
        click.echo(f"Modified antisense: {r.antisense}")
        click.echo(f"Efficacy score    : {r.efficacy_score:.2f} / 100 (biophysically adjusted)")
        click.echo(f"Delta vs parent   : {r.delta_score:+.2f}")
        click.echo(f"Efficacy label    : {r.efficacy_label}")
        p = getattr(r, 'biophysics', None)
        if p:
            click.echo(f"\nBiophysical penalties (raw score reduction):")
            for k, v in p.items():
                click.echo(f"  {k}: -{v:.1f}")
    else:
        click.echo("No result generated.")


if __name__ == "__main__":
    cli()
