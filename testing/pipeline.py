import os
import json
import csv
import logging
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

from typing import List, Sequence

from mechanism.loader import Mechanism
from reactor.batch import BatchResult, run_constant_pressure
from metaheuristics.ga import run_ga, GAOptions
from progress_variable import pv_error_aligned, progress_variable
from metrics import ignition_delay
from graph.construction import build_species_graph
from gnn.models import train_gnn, predict_scores
from visualizations import (
    plot_mole_fraction,
    plot_ignition_delays,
    plot_convergence,
    plot_pv_errors,
)

logger = logging.getLogger(__name__)

def evaluate_selection(selection, base_mech, Y0, tf, full_res, weights, critical_idxs):
    size_frac = float(selection.sum()) / len(selection)
    reason = ""

    if critical_idxs and (selection.sum() < max(4, len(critical_idxs)) or np.any(selection[critical_idxs] == 0)):
        return -1e6, 1.0, 1.0, 1.0, "missing_critical"

    keep = [base_mech.species_names[i] for i, bit in enumerate(selection) if bit]
    mech = Mechanism(base_mech.file_path)
    remove = [s for s in mech.species_names if s not in keep]

    try:
        if remove:
            mech.remove_species(remove)
        res = run_constant_pressure(
            mech.solution,
            full_res.temperature[0],
            ct.one_atm,
            Y0,
            tf,
            nsteps=len(full_res.time),
        )
    except Exception as e:
        return -1e6, 1.0, 1.0, 1.0, f"sim_failed:{type(e).__name__}"

    temp_rise = res.temperature.max() - res.temperature[0]
    species_change = np.max(np.abs(res.mass_fractions - res.mass_fractions[0]))
    if temp_rise < 10 or species_change < 1e-6:
        return -1e6, 1.0, 1.0, 1.0, "no_reaction"

    err = pv_error_aligned(
        full_res.mass_fractions,
        res.mass_fractions,
        base_mech.species_names,
        mech.species_names,
        np.asarray(weights)
    )
    delay_full = full_res.ignition_delay or ignition_delay(full_res.time, full_res.temperature)
    delay_red = ignition_delay(res.time, res.temperature)
    delay_diff = abs(delay_red - delay_full) / max(delay_full, 1e-12)
    fitness = -err - 10.0 * delay_diff - 0.1 * size_frac
    return fitness, err, delay_diff, 0.1 * size_frac, reason


def run_ga_reduction(mech_path, Y0, tf, steps):
    mech = Mechanism(mech_path)
    full = run_constant_pressure(mech.solution, 1500.0, ct.one_atm, Y0, tf, nsteps=steps)
    genome_len = len(mech.species_names)

    weights = np.ones(genome_len)
    if os.path.exists("data/species_weights.json"):
        with open("data/species_weights.json") as f:
            wmap = json.load(f)
            weights = np.array([wmap.get(s, 1e-3) for s in mech.species_names])

    # === Construct graph and predict scores ===
    G = build_species_graph(mech.solution)
    
    # Normalize weights before passing to GNN
    weights_norm = weights / (np.max(weights) + 1e-8)

    gnn_model = train_gnn(G, mech.solution, dict(zip(mech.species_names, weights_norm)), epochs=50)
    
    os.makedirs("results", exist_ok=True)

    scores = predict_scores(
        model=gnn_model,
        G=G,
        solution=mech.solution,
        save_path=os.path.join("results", "gnn_scores.csv")
    )

    scores_arr = np.array([scores[s] for s in mech.species_names])
    order = np.argsort(scores_arr)

    critical = [s for s in ["CH4", "O2", "N2"] if s in mech.species_names]
    critical_idxs = [mech.species_names.index(s) for s in critical]

    pop_size = 30
    k = max(1, int(0.5 * genome_len))
    seed = np.zeros(genome_len, dtype=int)
    seed[order[-k:]] = 1
    seed[critical_idxs] = 1

    init_pop = np.zeros((pop_size, genome_len), dtype=int)
    for i in range(pop_size):
        individual = seed.copy()
        flips = np.random.choice(genome_len, int(0.3 * genome_len), replace=False)
        individual[flips] ^= 1
        init_pop[i] = individual


    eval_fn = lambda sel: evaluate_selection(sel, mech, Y0, tf, full, weights, critical_idxs)

    sel, hist, debug = run_ga(
        genome_len,
        eval_fn,
        GAOptions(population_size=pop_size, generations=25, min_species=max(4, len(critical_idxs))),
        return_history=True,
        initial_population=init_pop,
        return_debug=True,
        fixed_indices=critical_idxs,
    )

    with open("best_selection.txt", "w") as f:
        f.write(",".join(map(str, sel.tolist())))

    return ["GA"], [sel], [hist], debug, full, weights


def full_pipeline(mech_path: str, out_dir: str, steps: int = 200, tf: float = 1.0):
    Y0 = {"CH4": 0.5, "O2": 1.0, "N2": 3.76}
    names, sols, hists, debug, full, weights = run_ga_reduction(mech_path, Y0, tf, steps)
    mech = Mechanism(mech_path)
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "ga_fitness.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["generation", "fitness"])
        for i, val in enumerate(hists[0]):
            writer.writerow([i, val])

    with open(os.path.join(out_dir, "debug_fitness.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["generation", "individual", "pv_error", "delay_diff", "size_penalty", "fitness", "reason"])
        for g, gen in enumerate(debug):
            for i, d in enumerate(gen):
                writer.writerow([g, i, *d[:-1], d[-1]])

    best_idx = int(np.argmax(hists[0]))
    best_sel = sols[0]
    keep = [mech.species_names[i] for i, bit in enumerate(best_sel) if bit]

    red_mech = Mechanism(mech_path)
    red_mech.remove_species([s for s in red_mech.species_names if s not in keep])

    red = run_constant_pressure(red_mech.solution, 1500.0, ct.one_atm, Y0, tf, nsteps=len(full.time))

    key_species = [s for s in ["CH4", "O2", "CO2"] if s in mech.species_names]
    idxs = [mech.species_names.index(s) for s in key_species]
    plot_profiles(full, red, idxs, key_species, os.path.join(out_dir, "profiles"))
    plot_ignition_delays([ignition_delay(red.time, red.temperature)], names, os.path.join(out_dir, "ignition_delay"))
    plot_convergence(hists, names, os.path.join(out_dir, "convergence"))
    plot_pv_errors(
        [pv_error_aligned(full.mass_fractions, red.mass_fractions, mech.species_names, red_mech.species_names, weights)],
        names,
        os.path.join(out_dir, "pv_error")
    )


def plot_profiles(full_res: BatchResult, red_res: BatchResult, species_idx: Sequence[int], species_names: Sequence[str], out_base: str) -> None:
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
