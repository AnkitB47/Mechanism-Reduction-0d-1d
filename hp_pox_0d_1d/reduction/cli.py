#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path


def add_parent_to_syspath():
    # Allow importing legacy GA+GNN modules located at the parent repo root
    repo_root = Path(__file__).resolve().parents[2].parent  # D:/Thesis/Master_Thesis
    if repo_root.exists():
        sys.path.insert(0, str(repo_root))
    else:
        print(f"Warning: legacy modules root not found at {repo_root.as_posix()} — some commands may degrade to no-op.")


def cmd_prepare_data(args: argparse.Namespace) -> int:
    from .dataset import prepare_dataset
    return prepare_dataset(args.cases_root, args.cases, Path(args.out))


def cmd_train_gnn(args: argparse.Namespace) -> int:
    add_parent_to_syspath()
    from .gnn_bridge import train_or_load_gnn
    return train_or_load_gnn(Path(args.data), Path(args.config), Path(args.workdir))


def cmd_score(args: argparse.Namespace) -> int:
    add_parent_to_syspath()
    from .gnn_bridge import score_reactions
    return score_reactions(Path(args.data), Path(args.workdir), mech_path=Path(args.mech))


def cmd_reduce(args: argparse.Namespace) -> int:
    add_parent_to_syspath()
    from .ga_bridge import run as ga_run
    from .ga_bridge import parse_budget
    # Cap BLAS threads for process-based parallelism
    import os
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    os.environ.setdefault('MKL_NUM_THREADS', '1')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    # Default cases if not provided
    if not args.cases:
        args.cases = ['A_case1_richter','A_case4_richter','B_case1_134L','B_case4_134L']
    budget = parse_budget(args.budget)
    out_mech = Path(args.out_mech)
    out_mech.parent.mkdir(parents=True, exist_ok=True)
    return ga_run(
        base_mech=Path(args.base_mech),
        out_mech=out_mech,
        cases_root=Path(args.cases_root),
        cases=args.cases,
        scores_path=Path(args.scores) if args.scores else None,
        budget=budget,
        workers=int(args.workers),
        fast_cells=int(args.fast_cells) if args.fast_cells else None,
        pop=int(args.pop) if args.pop else None,
        gens=int(args.gens) if args.gens else None,
        elite=int(args.elite) if args.elite else None,
        mut=float(args.mut) if args.mut else None,
        cx=float(args.cx) if args.cx else None,
        resume=Path(args.resume) if args.resume else None,
        target_species_band=tuple(args.target_species_band) if args.target_species_band else (28, 32),
        size_weight=float(args.alpha) if args.alpha else 30.0
    )


def cmd_validate(args: argparse.Namespace) -> int:
    from .evaluator import validate_mechanism
    return validate_mechanism(Path(args.mech), Path(args.cases_root), args.cases, Path(args.out))


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog='hp_pox_0d_1d.reduction', description='GA+GNN mechanism reduction CLI')
    sub = p.add_subparsers(dest='cmd', required=True)

    # prepare-data
    sp = sub.add_parser('prepare-data', help='Build dataset tensors from axial CSVs and KPIs')
    sp.add_argument('--cases', nargs='+', required=True, help='List of case IDs')
    sp.add_argument('--out', required=True, help='Output dataset directory')
    sp.add_argument('--cases-root', default='hp_pox_results_sandiego', help='Root folder containing <case>/axial_profiles.csv')
    sp.set_defaults(func=cmd_prepare_data)

    # train-gnn
    sp = sub.add_parser('train-gnn', help='Train or load the GNN surrogate')
    sp.add_argument('--data', required=True, help='Dataset directory produced by prepare-data')
    sp.add_argument('--config', required=True, help='GNN config YAML')
    sp.add_argument('--workdir', required=True, help='Workspace directory for checkpoints')
    sp.set_defaults(func=cmd_train_gnn)

    # score
    sp = sub.add_parser('score', help='Predict reaction scores aligned to a mechanism')
    sp.add_argument('--data', required=True, help='Dataset directory produced by prepare-data')
    sp.add_argument('--workdir', required=True, help='Workspace directory with trained GNN')
    sp.add_argument('--mech', required=True, help='Mechanism path to align reaction scores')
    sp.set_defaults(func=cmd_score)

    # reduce
    sp = sub.add_parser('reduce', help='Run GA+GNN reduction loop')
    sp.add_argument('--base-mech', required=True, help='Path to full mechanism YAML')
    sp.add_argument('--scores', required=False, help='Path to rxn_scores.npy (optional)')
    sp.add_argument('--constraints', required=True, help='Constraints YAML')
    sp.add_argument('--out-mech', required=True, help='Output reduced mechanism YAML path')
    sp.add_argument('--cases', nargs='+', required=False, help='Optional list of case directories')
    sp.add_argument('--cases-root', default='hp_pox_results_sandiego', help='Root folder containing <case>/axial_profiles.csv')
    sp.add_argument('--seed', type=int, default=42)
    sp.add_argument('--budget', required=False, default='generations=40,pop=60,topK=3,eval_every=5', help='Budget string')
    sp.add_argument('--workers', required=False, default='1', help='Number of worker processes for evaluation')
    sp.add_argument('--fast-cells', required=False, help='Axial cells for fast GA evaluation (physics unchanged)')
    sp.add_argument('--pop', required=False, help='GA population size override')
    sp.add_argument('--gens', required=False, help='GA generations override')
    sp.add_argument('--elite', required=False, help='GA elitism count')
    sp.add_argument('--mut', required=False, help='GA mutation rate')
    sp.add_argument('--cx', required=False, help='GA crossover rate')
    sp.add_argument('--resume', required=False, help='Directory to resume GA run from (checkpoints)')
    sp.set_defaults(func=cmd_reduce)

    # validate
    sp = sub.add_parser('validate', help='Validate a reduced mechanism via PSR→PFR')
    sp.add_argument('--mech', required=True, help='Path to reduced mechanism YAML')
    sp.add_argument('--cases', nargs='+', required=True, help='List of case IDs')
    sp.add_argument('--cases-root', default='hp_pox_results_sandiego', help='Root folder containing <case>/axial_profiles.csv')
    sp.add_argument('--out', required=True, help='Output directory')
    sp.set_defaults(func=cmd_validate)

    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rc = args.func(args)
    return rc if isinstance(rc, int) else 0


if __name__ == '__main__':
    raise SystemExit(main())


