"""Utilities for running the full reduction pipeline and plotting results."""

import csv
import json
import os
from typing import Callable, List, Tuple, Sequence, Iterable
import numpy as np
import cantera as ct
from mechanism.loader import Mechanism
from reactor.batch import BatchResult, run_constant_pressure
from progress_variable import progress_variable
from metaheuristics.ga import run_ga, GAOptions
from metrics import pv_error, ignition_delay
from graph.construction import build_species_graph
from gnn.models import train_dummy_gnn, predict_scores
from visualizations import (
    plot_mole_fraction,
    plot_ignition_delays,
    plot_convergence,
    plot_pv_errors,
)


def evaluate_selection(
    selection: np.ndarray,
    base_mech: Mechanism,
    Y0: dict,
    tf: float,
    full_res,
    scores: np.ndarray | None = None,
    weights: Sequence[float] | None = None,
) -> float:
    """Return fitness for a given species selection.

    Parameters
    ----------
    selection:
        Binary vector indicating which species to keep.
    base_mech:
        Reference mechanism used as a template.
    Y0, tf:
        Initial composition and final time for the reactor simulation.
    full_res:
        Result from the full mechanism simulation used for comparison.
    scores:
        Optional importance scores for penalising removal of species.
    weights:
        Progress-variable weights.
    """

    idx_keep = [i for i, bit in enumerate(selection) if bit]
    keep = [base_mech.species_names[i] for i in idx_keep]
    mech = Mechanism(base_mech.file_path)
    remove = [s for s in mech.species_names if s not in keep]
    try:
        if remove:
            mech.remove_species(remove)
        res = run_constant_pressure(mech.solution, full_res.temperature[0], ct.one_atm, Y0, tf, nsteps=len(full_res.time))
    except Exception:
        return -1e6
    if weights is None:
        weights_sel = np.ones(len(idx_keep))
    else:
        weights_sel = np.array([weights[i] for i in idx_keep])
    pv_full = progress_variable(full_res.mass_fractions[:, idx_keep], weights_sel)
    pv_red = progress_variable(res.mass_fractions, weights_sel)
    err = pv_error(pv_full, pv_red)
    if scores is None:
        size_penalty = selection.sum() / len(selection)
    else:
        size_penalty = (1 - selection) * scores
        size_penalty = size_penalty.sum() / len(selection)
    return -(err + 0.01 * float(size_penalty))


def run_ga_reduction(
    mech_path: str,
    Y0: dict,
    tf: float,
    steps: int = 200,
) -> Tuple[List[str], List[np.ndarray], List[List[float]], BatchResult, np.ndarray]:
    """Run the GA-based reduction on the given mechanism."""

    mech = Mechanism(mech_path)
    full = run_constant_pressure(mech.solution, 1000.0, ct.one_atm, Y0, tf, nsteps=steps)
    genome_len = len(mech.species_names)

    weights_path = os.path.join("data", "species_weights.json")
    if os.path.exists(weights_path):
        with open(weights_path) as f:
            weight_map = json.load(f)
        weights = np.array([weight_map.get(s, 0.0) for s in mech.species_names])
    else:
        weights = np.ones(genome_len)

    G = build_species_graph(mech.solution)
    gnn_model = train_dummy_gnn(G, epochs=5)
    scores_dict = predict_scores(gnn_model, G)
    scores = np.array([scores_dict[s] for s in mech.species_names])

    eval_fn = lambda sel: evaluate_selection(sel, mech, Y0, tf, full, scores, weights)

    ga_sel, ga_hist = run_ga(
        genome_len,
        eval_fn,
        GAOptions(population_size=8, generations=4),
        return_history=True,
    )

    names = ["GA"]
    solutions = [ga_sel]
    history = [ga_hist]
    return names, solutions, history, full, weights


def save_metrics(
    names: List[str],
    sols: List[np.ndarray],
    full_res: BatchResult,
    mech: Mechanism,
    Y0: dict,
    tf: float,
    out_dir: str,
    weights: Sequence[float] | None = None,
) -> Tuple[List[BatchResult], List[float], List[float], List[int]]:
    """Evaluate reduced mechanisms and store metric CSV files."""

    os.makedirs(out_dir, exist_ok=True)
    results: List[BatchResult] = []
    pv_errors: List[float] = []
    delays: List[float] = []
    sizes: List[int] = []
    if weights is None:
        weights = np.ones(full_res.mass_fractions.shape[1])
    full_delay = ignition_delay(full_res.time, full_res.temperature)

    metrics_rows = []
    for name, sel in zip(names, sols):
        idx_keep = [i for i, b in enumerate(sel) if b]
        keep = [mech.species_names[i] for i in idx_keep]
        red_mech = Mechanism(mech.file_path)
        remove = [s for s in red_mech.species_names if s not in keep]
        try:
            if remove:
                red_mech.remove_species(remove)
            res = run_constant_pressure(
                red_mech.solution,
                1000.0,
                ct.one_atm,
                Y0,
                tf,
                nsteps=len(full_res.time),
            )
            weights_sel = np.array([weights[i] for i in idx_keep])
            pv_full = progress_variable(full_res.mass_fractions[:, idx_keep], weights_sel)
            pv_red = progress_variable(res.mass_fractions, weights_sel)
            err = pv_error(pv_full, pv_red)
            delay = ignition_delay(res.time, res.temperature)
            size_red = int(sel.sum())
            metrics_rows.append([
                name,
                err,
                abs(delay - full_delay) / full_delay,
                size_red,
            ])
            results.append(res)
            pv_errors.append(err)
            delays.append(delay)
            sizes.append(size_red)
        except Exception:
            metrics_rows.append([name, 1.0, 1.0, 0])
            results.append(full_res)
            pv_errors.append(1.0)
            delays.append(full_delay)
            sizes.append(0)

    with open(os.path.join(out_dir, "metrics.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "algorithm",
            "pv_error",
            "ignition_delay_error",
            "species_retained",
        ])
        writer.writerows(metrics_rows)

    with open(os.path.join(out_dir, "ignition_delay.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["algorithm", "ignition_delay_s"])
        writer.writerow(["Full", full_delay])
        writer.writerows(zip(names, delays))

    with open(os.path.join(out_dir, "pv_error.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["algorithm", "pv_error"])
        writer.writerows(zip(names, pv_errors))

    with open(os.path.join(out_dir, "species_retained.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["algorithm", "num_species"])
        writer.writerows(zip(names, sizes))

    return results, pv_errors, delays, sizes


import matplotlib.pyplot as plt

def plot_profiles(
    full_res: BatchResult,
    red_res: BatchResult,
    species_idx: Sequence[int],
    species_names: Sequence[str],
    out_base: str,
) -> None:
    """Plot mole-fraction profiles of selected species."""

    fig, ax = plt.subplots()
    for idx, name in zip(species_idx, species_names):
        ax.plot(full_res.time, full_res.mass_fractions[:, idx], label=f"{name} full")
        ax.plot(red_res.time, red_res.mass_fractions[:, idx], "--", label=f"{name} red")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mass Fraction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    plt.close(fig)

def full_pipeline(mech_path: str, out_dir: str, steps: int = 200):
    """Run the complete testing pipeline and generate plots/CSVs."""

    Y0 = {"CH4": 0.5, "O2": 1.0, "N2": 3.76}
    tf = 0.02

    names, solutions, histories, full, weights = run_ga_reduction(mech_path, Y0, tf, steps)
    mech = Mechanism(mech_path)

    results, pv_err, delays, sizes = save_metrics(names, solutions, full, mech, Y0, tf, out_dir, weights)

    # write fitness history CSVs
    for name, hist in zip(names, histories):
        if not hist:
            continue
        path = os.path.join(out_dir, f"{name.lower()}_fitness.csv")
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["generation", "fitness"])
            for i, val in enumerate(hist):
                writer.writerow([i, val])

    # choose the algorithm with minimum PV error for profile plot
    best_idx = int(np.argmin(pv_err)) if pv_err else 0
    key_species = [s for s in ["CH4", "O2", "CO2"] if s in mech.species_names]
    idxs = [mech.species_names.index(s) for s in key_species]
    if results:
        plot_profiles(full, results[best_idx], idxs, key_species, os.path.join(out_dir, "profiles"))

    plot_ignition_delays(delays, names, os.path.join(out_dir, "ignition_delay"))
    plot_convergence(histories, names, os.path.join(out_dir, "convergence"))
    plot_pv_errors(pv_err, names, os.path.join(out_dir, "pv_error"))
