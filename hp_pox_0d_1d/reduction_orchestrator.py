from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Tuple

from reduction.ga_bridge import run_ga_gnn
from reduction.evaluator import evaluate, calculate_metrics_penalty


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Hybrid GA–GNN mechanism reduction orchestrator")
    parser.add_argument("--base_mech", required=True, help="Path to base mechanism YAML")
    parser.add_argument("--cases_root", required=True, help="Directory containing validation cases")
    parser.add_argument(
        "--case_list",
        nargs="+",
        required=True,
        help="List of case identifiers to evaluate (PSR/PFR)",
    )
    parser.add_argument("--output_dir", default="runs", help="Directory to store GA outputs")
    parser.add_argument("--seed_mech", default=None, help="Optional seed mechanism for initialization")
    parser.add_argument("--target_min", type=int, default=28, help="Minimum species target")
    parser.add_argument("--target_max", type=int, default=32, help="Maximum species target")
    parser.add_argument("--pop", type=int, default=40, help="GA population size")
    parser.add_argument("--gens", type=int, default=120, help="Number of GA generations")
    parser.add_argument("--elitism", type=int, default=2, help="Number of elites retained each generation")
    parser.add_argument("--mutation", type=float, default=0.15, help="Species mutation rate")
    parser.add_argument(
        "--soft_species",
        nargs="*",
        default=None,
        help="Optional soft-protected species list (space separated)",
    )
    return parser.parse_args()


def summarize_results(summary: dict, case_list: List[str], cases_root: Path) -> None:
    reduced_path = Path(summary["reduced_yaml_path"])
    metrics_path = Path(summary["metrics_json_path"])
    print("=" * 72)
    print("GA–GNN REDUCTION COMPLETE")
    print("=" * 72)
    print(f"Reduced mechanism : {reduced_path}")
    print(f"Metrics JSON      : {metrics_path}")
    print(f"GA trace CSV      : {summary['ga_trace_csv_path']}")
    print(f"Final species / reactions : {summary['closure']['n_species']} / {summary['closure']['n_reactions']}")
    if metrics_path.exists():
        metrics = json.loads(metrics_path.read_text())
        penalty = calculate_metrics_penalty(metrics, n_species=summary['closure']['n_species'])
        print(f"Final KPI penalty : {penalty:.3f}")
        per_case = metrics.get("per_case", {})
        for case in case_list:
            data = per_case.get(case, {})
            rmse = data.get("rmse", {})
            print(
                f"  {case:30s} | RMSE_T={rmse.get('T', float('nan')):.2f} "
                f"RMSE_X_H2={rmse.get('X_H2', float('nan')):.3f} "
                f"RMSE_X_CO={rmse.get('X_CO', float('nan')):.3f}"
            )
    else:
        print("[warn] metrics.json not found; skipping detailed summary")


def main() -> None:
    args = parse_args()
    target_band: Tuple[int, int] = (args.target_min, args.target_max)

    summary = run_ga_gnn(
        base_mech_path=args.base_mech,
        output_dir=args.output_dir,
        cases_root=args.cases_root,
        case_list=args.case_list,
        target_band=target_band,
        pop=args.pop,
        gens=args.gens,
        elitism=args.elitism,
        mutation_rate=args.mutation,
        seed_mech_path=args.seed_mech,
        soft_species=args.soft_species,
    )

    summarize_results(summary, args.case_list, Path(args.cases_root))


if __name__ == "__main__":
    main()
