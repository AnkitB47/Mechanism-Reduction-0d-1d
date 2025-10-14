from __future__ import annotations

from pathlib import Path
from typing import List
import json
import random

from .ga_bridge import parse_budget
from .mechanism_io import copy_mechanism
from .evaluator import evaluate


def run_reduction(mech_path: Path, cases_root: Path, cases: List[str], budget: str, constraints_path: Path, outdir: Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    bud = parse_budget(budget)

    # Demo: copy original mechanism as reduced_mech.yaml (no-op pruning) to prove wiring
    red_path = outdir / 'reduced_mech.yaml'
    copy_mechanism(mech_path, red_path)

    # Evaluate reduced vs baseline CSVs (baseline taken from cases_root)
    metrics = evaluate(red_path, cases_root, cases, outdir / 'validate')

    # Scorecard
    scorecard = {
        'mechanism_in': str(mech_path.resolve()),
        'mechanism_out': str(red_path.resolve()),
        'budget': bud.__dict__,
        'reaction_reduction_pct': 0.0,  # stub until masking applied
        'metrics': metrics
    }
    with open(outdir / 'scorecard.json', 'w') as f:
        json.dump(scorecard, f, indent=2)

    print(f"Reduced mechanism written to: {red_path.resolve().as_posix()}")
    print(f"Scorecard: {str((outdir / 'scorecard.json').resolve().as_posix())}")
    return 0


