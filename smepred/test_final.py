"""Verify final beam search behavior."""
import sys
import warnings; warnings.filterwarnings('ignore')
from collections import Counter
from src.modification_engine import multi_mod_scan

sense = 'GGAAAUAGACACCAAAUCUUA'
as_ = 'UAAGAUUUGGUGUCUAUUUCC'

results = multi_mod_scan(sense, as_, max_mods=14, beam_width=30, model_key='B', full_scan=True)
mod_counts = Counter()
for r in results[:20]:
    n = r.mod_positions.count(',') + 1 if r.mod_positions else 1
    mod_counts[n] += 1
last_round = max((r.mod_positions.count(',') + 1 if r.mod_positions else 1) for r in results)
print(f'Total: {len(results)}, Last round: {last_round}')
print(f'Top-20 mod distribution: {dict(sorted(mod_counts.items()))}')
for i, r in enumerate(results[:3]):
    n = r.mod_positions.count(',') + 1 if r.mod_positions else 1
    syms = r.mod_symbol.replace('+', '')
    print(f'  #{i+1}: adj={r.efficacy_score:.1f} mods={n:2d} syms={syms}')
print('---')
# Test exotic micro-penalty differentiation
sense2 = 'AAGCUGGCCUCAGUUAACUGA'
as2 = 'UCAGUUAACUGAGGCCAGCUU'
r2 = multi_mod_scan(sense2, as2, max_mods=14, beam_width=30, model_key='B', full_scan=True)
print(f'Seq2: Total={len(r2)}, Last round={max((r.mod_positions.count(",") + 1 if r.mod_positions else 1) for r in r2)}')
print(f'Seq2 top-1: adj={r2[0].efficacy_score:.1f} mods={(r2[0].mod_positions.count(",") + 1 if r2[0].mod_positions else 1)}')
import os; os.remove(__file__)
