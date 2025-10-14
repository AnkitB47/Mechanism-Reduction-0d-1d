#!/usr/bin/env python3
"""
HP-POX end-to-end orchestrator.

Examples
--------
Run with the default (discovered) model entry using the bundled GRI30 mechanism::

    python hp_pox_0d_1d/tools/run_full_vs_reduced.py \\
        --mech gri30 \\
        --cases A_case1_richter A_case4_richter \\
        --out comparisons/$(date +%Y%m%d_%H%M) \\
        --resume --skip-full

Run against a custom mechanism path::

    python hp_pox_0d_1d/tools/run_full_vs_reduced.py \\
        --mech mechanisms/aramco/aramco2.yaml \\
        --cases A_case1_richter A_case4_richter \\
        --out comparisons/$(date +%Y%m%d_%H%M) \\
        --force --resume
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import cantera as ct

if __package__ is None or __package__ == "":
    sys.path.append(str(Path(__file__).resolve().parents[1]))
    from tools.mech_resolver import list_available_mechanisms, resolve_mechanism  # type: ignore
    from tools.config_templates import write_full_config  # type: ignore
else:
    from .mech_resolver import list_available_mechanisms, resolve_mechanism
    from .config_templates import write_full_config

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIGS_DIR = REPO_ROOT / "configs"
DEFAULT_FULL_RESULTS_DIR = REPO_ROOT / "hp_pox_results_full"
DEFAULT_REDUCED_RESULTS_DIR = REPO_ROOT / "hp_pox_results_reduced"
DEFAULT_TEMPLATE = REPO_ROOT / "hp_pox_0d_1d" / "config.yaml"

DEFAULT_POP = 40
DEFAULT_GENS = 120
DEFAULT_ELITISM = 2
DEFAULT_MUTATION = 0.15


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run baseline, reduction, and comparison for HP-POX mechanisms.")
    parser.add_argument("--mech", required=True, help="Mechanism alias or YAML path (e.g. 'gri30', 'aramco', 'path/to/mech.yaml').")
    parser.add_argument("--cases", nargs="+", required=True, help="Case identifiers to evaluate (e.g. A_case1_richter).")
    parser.add_argument("--out", default=None, help="Comparison output directory (default: comparisons/<timestamp>).")
    parser.add_argument("--model-entry", default=None, help="Override the model runner script (must accept --config).")
    parser.add_argument("--full-config", default=None, help="Optional baseline config YAML path (auto-generated if missing).")
    parser.add_argument("--full-results-dir", default=str(DEFAULT_FULL_RESULTS_DIR), help="Directory for baseline results.")
    parser.add_argument("--reduced-results-dir", default=str(DEFAULT_REDUCED_RESULTS_DIR), help="Directory for reduced results.")
    parser.add_argument("--force", action="store_true", help="Force rerun of the baseline even if results exist.")
    parser.add_argument("--resume", action="store_true", help="Reuse the latest GA run if available.")
    parser.add_argument("--skip-full", action="store_true", help="Skip the baseline step (assumes results already exist).")
    parser.add_argument("--pop", type=int, default=DEFAULT_POP, help="GA population size.")
    parser.add_argument("--gens", type=int, default=DEFAULT_GENS, help="GA generations.")
    parser.add_argument("--elitism", type=int, default=DEFAULT_ELITISM, help="GA elitism count.")
    parser.add_argument("--mutation", type=float, default=DEFAULT_MUTATION, help="GA mutation rate.")
    parser.add_argument("--dry-run", action="store_true", help="Plan the stages and create stub outputs without running heavy simulations.")
    parser.add_argument("--list-mechanisms", action="store_true", help="List discovered mechanism aliases and exit.")
    return parser.parse_args()


def resolve_path(path_str: str) -> Path:
    path = Path(path_str)
    return path if path.is_absolute() else (REPO_ROOT / path)


def run_command(cmd: List[str], hint: Optional[str] = None, cwd: Optional[Path] = None) -> None:
    workdir = cwd or REPO_ROOT
    cmd_display = " ".join(cmd)
    print(f"[cmd] ({workdir}) {cmd_display}")
    result = subprocess.run(cmd, cwd=workdir, check=False)
    if result.returncode != 0:
        message = f"Command failed (exit {result.returncode}): {cmd_display}"
        if hint:
            message += f"\nHint: {hint}"
        raise RuntimeError(message)


def slugify(value: str) -> str:
    stem = Path(value).stem if value.endswith(".yaml") else value
    cleaned = "".join(ch.lower() if ch.isalnum() else "_" for ch in stem)
    cleaned = cleaned.strip("_")
    return cleaned or "mechanism"


def _candidate_model_entries() -> List[Path]:
    return [
        REPO_ROOT / "hp_pox_0d_1d" / "run_hp_pox.py",
        REPO_ROOT / "hp_pox_model.py",
        REPO_ROOT / "hp_pox_0d_1d" / "hp_pox_model.py",
        REPO_ROOT / "hp_pox_0d_1d" / "hp_pox_physics" / "hp_pox_model.py",
        REPO_ROOT / "hp_pox_0d_1d" / "tools" / "simulate_baseline.py",
    ]


def validate_model_entry(entry: Path) -> bool:
    if not entry.exists():
        return False
    try:
        result = subprocess.run(
            [sys.executable, str(entry), "--help"],
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except Exception:
        return False
    return result.returncode == 0


def _count_mechanism(mech_path: Path) -> Tuple[int, int]:
    try:
        sol = ct.Solution(str(mech_path))
    except Exception as exc:
        print(f"[warn] unable to load mechanism {mech_path}: {exc}")
        return 0, 0
    return int(sol.n_species), int(sol.n_reactions)


def _augment_case_diagnostics(results_dir: Path, cases: List[str], mech_path: Path, n_species: int, n_reactions: int) -> None:
    for case in cases:
        diag_path = results_dir / case / "diagnostics.json"
        if not diag_path.exists():
            continue
        try:
            with diag_path.open("r") as f:
                data = json.load(f) or {}
        except Exception:
            continue
        data["mechanism_path"] = str(mech_path)
        data["n_species_mech"] = int(n_species)
        data["n_reactions_mech"] = int(n_reactions)
        with diag_path.open("w") as f:
            json.dump(data, f, indent=2)


def discover_model_entry(user_entry: Optional[str]) -> Path:
    if user_entry:
        entry = resolve_path(user_entry)
        if not validate_model_entry(entry):
            raise FileNotFoundError(
                f"Model entry script '{entry}' either does not exist or does not accept '--help'. "
                "Pass a valid script with --model-entry."
            )
        return entry

    for candidate in _candidate_model_entries():
        if validate_model_entry(candidate):
            print(f"[model] using entry script {candidate}")
            return candidate

    raise FileNotFoundError(
        "Unable to locate a valid model runner script. "
        "Pass one explicitly via --model-entry."
    )


def ensure_full_run(
    cases: List[str],
    force: bool,
    model_entry: Path,
    config_path: Path,
    results_dir: Path,
    dry_run: bool,
) -> None:
    needs_run = force or baseline_results_missing(cases, results_dir)
    if dry_run:
        print("[full] dry-run: baseline execution skipped.")
        return
    if not needs_run:
        print("[full] results found, reusing existing data")
        return

    cmd = [
        sys.executable,
        str(model_entry),
        "--config",
        str(config_path.resolve()),
    ]
    run_command(cmd, hint="Use --skip-full to reuse existing baseline results.", cwd=REPO_ROOT)

    for case in cases:
        case_dir = results_dir / case
        if not case_dir.is_dir():
            raise FileNotFoundError(f"Expected full results at {case_dir}")
        print(f"[full] done -> {case_dir}")


def baseline_results_missing(cases: List[str], root: Path) -> bool:
    return any(not (root / case).is_dir() for case in cases)


def run_ga_reduction(
    cases: List[str],
    resume: bool,
    base_mech_path: Path,
    full_results_dir: Path,
    pop: int,
    gens: int,
    elitism: int,
    mutation: float,
    dry_run: bool,
) -> Dict[str, object]:
    if dry_run:
        return create_stub_ga_summary(base_mech_path)

    if resume:
        last_dir = latest_run_dir()
        if last_dir:
            summary = load_summary(last_dir)
            reduced_path = Path(summary.get("reduced_yaml_path", ""))
            if reduced_path.exists():
                print(f"[GA] resume -> {last_dir}")
                return summary
            print("[GA] resume requested but reduced YAML missing in last run; re-running GA")

    orchestrator = REPO_ROOT / "hp_pox_0d_1d" / "reduction_orchestrator.py"
    cmd = [
        sys.executable,
        str(orchestrator),
        "--base_mech",
        str(base_mech_path.resolve()),
        "--cases_root",
        str(full_results_dir.resolve()),
        "--case_list",
    ] + cases + [
        "--target_min",
        "28",
        "--target_max",
        "32",
        "--pop",
        str(pop),
        "--gens",
        str(gens),
        "--elitism",
        str(elitism),
        "--mutation",
        str(mutation),
    ]
    run_command(cmd, hint="Use --resume to reuse the latest GA run.", cwd=REPO_ROOT)
    run_dir = latest_run_dir()
    if run_dir is None:
        raise RuntimeError("GA run completed but no runs/<timestamp> directory was created.")
    summary = load_summary(run_dir)
    print(f"[GA] completed -> {run_dir}")
    return summary


def latest_run_dir() -> Optional[Path]:
    runs_root = REPO_ROOT / "runs"
    if not runs_root.is_dir():
        return None
    candidates = [p for p in runs_root.iterdir() if p.is_dir()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def load_summary(run_dir: Path) -> Dict[str, object]:
    summary_path = run_dir / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary not found in {run_dir}")
    with summary_path.open("r") as f:
        return json.load(f)


def ensure_reduced_run(
    cases: List[str],
    model_entry: Path,
    config_path: Path,
    results_dir: Path,
    dry_run: bool,
) -> None:
    if dry_run:
        print("[reduced] dry-run: reduced execution skipped.")
        return
    needs_run = baseline_results_missing(cases, results_dir)
    if not needs_run:
        print("[reduced] results found, reusing existing data")
        return

    cmd = [sys.executable, str(model_entry), "--config", str(config_path.resolve())]
    run_command(cmd, cwd=REPO_ROOT)

    for case in cases:
        case_dir = results_dir / case
        if not case_dir.is_dir():
            raise FileNotFoundError(f"Expected reduced results at {case_dir}")
        print(f"[reduced] done -> {case_dir}")


def run_comparison(
    cases: List[str],
    out_dir: Path,
    reduced_yaml: Path,
    base_mech_path: Path,
    full_results_dir: Path,
    reduced_results_dir: Path,
    dry_run: bool,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return create_stub_comparison(out_dir, cases)

    cmd = [
        sys.executable,
        str(REPO_ROOT / "tools/compare_runs.py"),
        "--full",
        str(full_results_dir.resolve()),
        "--reduced",
        str(reduced_results_dir.resolve()),
        "--cases",
    ] + cases + [
        "--outdir",
        str(out_dir.resolve()),
        "--full-mech",
        str(base_mech_path.resolve()),
        "--reduced-mech",
        str(reduced_yaml.resolve()),
    ]
    run_command(cmd, cwd=REPO_ROOT)
    return out_dir / "summary_metrics.csv"


def load_summary_table(path: Path) -> List[Dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(f"Summary metrics CSV not found: {path}")
    rows: List[Dict[str, object]] = []
    with path.open("r") as f:
        headers = f.readline().strip().split(",")
        for line in f:
            values = line.strip().split(",")
            row = {headers[i]: values[i] if i < len(values) else "" for i in range(len(headers))}
            rows.append(row)
    return rows


def print_digest(rows: List[Dict[str, object]]) -> None:
    if not rows:
        print("[digest] summary table is empty")
        return

    header = [
        "CASE",
        "N_full",
        "N_red",
        "dT_max(K)",
        "RMSE(H2)",
        "RMSE(CO)",
        "CH4 slip(%)",
        "PASS",
    ]
    print()
    print("{:<20s} {:>6s} {:>6s} {:>10s} {:>9s} {:>9s} {:>12s} {:>6s}".format(*header))
    for row in rows:
        pass_flag = row.get("PASS", "").strip().lower() == "true"
        badge = "PASS" if pass_flag else "FAIL"
        print(
            "{:<20s} {:>6s} {:>6s} {:>10s} {:>9s} {:>9s} {:>12s} {:>6s}".format(
                row.get("case", ""),
                row.get("N_species_full", ""),
                row.get("N_species_reduced", ""),
                row.get("max_abs_dT_K", ""),
                row.get("RMSE_X_H2", ""),
                row.get("RMSE_X_CO", ""),
                row.get("CH4_slip_percent", ""),
                badge,
            )
        )


def create_stub_ga_summary(base_mech_path: Path) -> Dict[str, object]:
    runs_root = REPO_ROOT / "runs"
    tag = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = runs_root / f"dry_run_{tag}"
    final_eval_dir = run_dir / "final_eval"
    run_dir.mkdir(parents=True, exist_ok=True)
    final_eval_dir.mkdir(parents=True, exist_ok=True)

    reduced_yaml = run_dir / "reduced.yaml"
    reduced_yaml.write_text(f"# dry-run placeholder for {base_mech_path.name}\n")

    metrics_path = final_eval_dir / "metrics.json"
    metrics_path.write_text(json.dumps({"note": "dry-run metrics placeholder"}, indent=2))

    trace_path = run_dir / "ga_trace.csv"
    trace_path.write_text("generation,penalty,size,psr_pass,pfr_pass\n0,0.0,0,True,True\n")

    summary = {
        "runtime_sec": 0.0,
        "best_penalty": 0.0,
        "final_penalty": 0.0,
        "closure": {"n_species": 0, "n_reactions": 0},
        "reduced_yaml_path": str(reduced_yaml.resolve()),
        "final_metrics_path": str(metrics_path.resolve()),
        "ga_trace_csv_path": str(trace_path.resolve()),
    }
    summary_path = run_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    return summary


def create_stub_comparison(out_dir: Path, cases: List[str]) -> Path:
    rows = []
    for case in cases:
        case_dir = out_dir / case
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "figs").mkdir(parents=True, exist_ok=True)
        metrics = {
            "case": case,
            "max_abs_dT_K": 0.0,
            "RMSE_X_H2": 0.0,
            "RMSE_X_CO": 0.0,
            "PASS": True,
        }
        (case_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
        rows.append(
            {
                "case": case,
                "N_species_full": "",
                "N_reactions_full": "",
                "N_species_reduced": "",
                "N_reactions_reduced": "",
                "max_abs_dT_K": "0.0",
                "RMSE_X_H2": "0.0",
                "RMSE_X_CO": "0.0",
                "RMSE_X_CO2": "0.0",
                "RMSE_X_CH4": "0.0",
                "RMSE_X_H2O": "0.0",
                "RMSE_X_O2": "0.0",
                "RMSE_X_H2CO": "0.0",
                "CH4_slip_percent": "0.0",
                "H2_CO_rel_err": "0.0",
                "PASS": "True",
            }
        )
    summary_path = out_dir / "summary_metrics.csv"
    headers = [
        "case",
        "N_species_full",
        "N_reactions_full",
        "N_species_reduced",
        "N_reactions_reduced",
        "max_abs_dT_K",
        "RMSE_X_H2",
        "RMSE_X_CO",
        "RMSE_X_CO2",
        "RMSE_X_CH4",
        "RMSE_X_H2O",
        "RMSE_X_O2",
        "RMSE_X_H2CO",
        "CH4_slip_percent",
        "H2_CO_rel_err",
        "PASS",
    ]
    with summary_path.open("w") as f:
        f.write(",".join(headers) + "\n")
        for row in rows:
            f.write(",".join(str(row.get(h, "")) for h in headers) + "\n")
    return summary_path


def main() -> None:
    args = parse_args()

    if args.list_mechanisms:
        available = list_available_mechanisms()
        if not available:
            print("No mechanism aliases discovered. Drop YAML files under 'mechanisms/' or pass an explicit path.")
        else:
            print("Discovered mechanism aliases:")
            for family, paths in available.items():
                print(f"  {family}:")
                for path in paths:
                    print(f"    - {path}")
        return

    model_entry = discover_model_entry(args.model_entry)

    mech_path = resolve_mechanism(args.mech).resolve()
    slug = slugify(mech_path.stem)

    full_results_dir = resolve_path(args.full_results_dir)
    reduced_results_dir = resolve_path(args.reduced_results_dir)
    full_results_dir.mkdir(parents=True, exist_ok=True)
    reduced_results_dir.mkdir(parents=True, exist_ok=True)

    if args.full_config:
        full_config_path = resolve_path(args.full_config)
    else:
        full_config_path = CONFIGS_DIR / f"{slug}_full.yaml"
    full_config_path.parent.mkdir(parents=True, exist_ok=True)

    template_src = full_config_path if full_config_path.exists() else DEFAULT_TEMPLATE
    write_full_config(
        template_src if template_src.exists() else None,
        full_config_path,
        mech_path,
        cases_root=full_results_dir,
        output_root=full_results_dir,
    )
    print(f"[config] baseline config ready at {full_config_path}")

    if args.skip_full and args.dry_run:
        print("[orchestrator] --skip-full set -> baseline will be considered complete (dry-run).")
    elif args.skip_full:
        print("[orchestrator] --skip-full set -> baseline step skipped.")
    else:
        ensure_full_run(args.cases, args.force, model_entry, full_config_path, full_results_dir, args.dry_run)

    if not args.dry_run:
        n_species_full, n_reactions_full = _count_mechanism(mech_path)
        _augment_case_diagnostics(full_results_dir, args.cases, mech_path, n_species_full, n_reactions_full)

    ga_summary = run_ga_reduction(
        args.cases,
        args.resume,
        mech_path,
        full_results_dir,
        args.pop,
        args.gens,
        args.elitism,
        args.mutation,
        args.dry_run,
    )
    reduced_yaml = Path(ga_summary["reduced_yaml_path"]).resolve()

    if args.full_config:
        reduced_config_path = full_config_path.with_name(f"{slug}_reduced.yaml")
    else:
        reduced_config_path = CONFIGS_DIR / f"{slug}_reduced.yaml"
    write_full_config(
        full_config_path,
        reduced_config_path,
        reduced_yaml,
        cases_root=reduced_results_dir,
        output_root=reduced_results_dir,
    )
    print(f"[config] generated reduced config at {reduced_config_path}")

    ensure_reduced_run(args.cases, model_entry, reduced_config_path, reduced_results_dir, args.dry_run)

    if not args.dry_run:
        n_species_reduced, n_reactions_reduced = _count_mechanism(reduced_yaml)
        _augment_case_diagnostics(reduced_results_dir, args.cases, reduced_yaml, n_species_reduced, n_reactions_reduced)

    if args.out:
        comparison_dir = resolve_path(args.out)
    else:
        comparison_dir = REPO_ROOT / "comparisons" / dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    summary_csv = run_comparison(
        args.cases,
        comparison_dir,
        reduced_yaml,
        mech_path,
        full_results_dir,
        reduced_results_dir,
        args.dry_run,
    )

    rows = load_summary_table(summary_csv)
    print_digest(rows)
    print(f"[done] comparison artefacts -> {comparison_dir}")


if __name__ == "__main__":
    main()
