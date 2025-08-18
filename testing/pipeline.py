import os
import json
import csv
import logging
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

from typing import List, Sequence

from mechanism.loader import Mechanism
from mechanism.mix import methane_air_mole_fractions, mole_to_mass_fractions
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


def evaluate_selection(
    selection,
    base_mech,
    Y0,
    tf,
    full_res,
    weights,
    critical_idxs,
    T0,
    p0,
    steps,
    log_times,
):
    size_frac = float(selection.sum()) / len(selection)
    reason = ""

    if critical_idxs and (selection.sum() < max(4, len(critical_idxs)) or np.any(selection[critical_idxs] == 0)):
        return -1e6, 1.0, 1.0, 1.0, "missing_critical", selection.sum(), 0.0, 0.0, 0.0

    keep = [base_mech.species_names[i] for i, bit in enumerate(selection) if bit]
    logger.info("keep_count=%d first10=%s", len(keep), keep[:10])
    mech = Mechanism(base_mech.file_path)
    remove = [s for s in mech.species_names if s not in keep]

    try:
        if remove:
            mech.remove_species(remove)
        logger.info(
            "reduced_mech: %d species, %d reactions",
            len(mech.species_names),
            len(mech.reactions()),
        )
        assert set(mech.species_names) == set(keep)
        res = run_constant_pressure(
            mech.solution,
            T0,
            p0,
            Y0,
            tf,
            nsteps=steps,
            use_mole=False,
            log_times=log_times,
        )
    except RuntimeError as e:
        if "no_reaction" in str(e):
            return -1e6, 1.0, 1.0, 1.0, str(e), selection.sum(), 0.0, 0.0, 0.0
        return -1e6, 1.0, 1.0, 1.0, f"sim_failed:{type(e).__name__}", selection.sum(), 0.0, 0.0, 0.0
    except Exception as e:
        return -1e6, 1.0, 1.0, 1.0, f"sim_failed:{type(e).__name__}", selection.sum(), 0.0, 0.0, 0.0

    err = pv_error_aligned(
        full_res.mass_fractions,
        res.mass_fractions,
        base_mech.species_names,
        mech.species_names,
        np.asarray(weights),
    )
    delay_full = full_res.ignition_delay or ignition_delay(full_res.time, full_res.temperature)[0]
    delay_red, _ = ignition_delay(res.time, res.temperature)
    delay_diff = abs(delay_red - delay_full) / max(delay_full, 1e-12)
    delta_T = float(res.temperature[-1] - res.temperature[0])
    delta_Y = float(np.max(np.abs(res.mass_fractions[-1] - res.mass_fractions[0])))
    fitness = -err - 10.0 * delay_diff - 2.0 * max(0.0, size_frac - 0.4)
    logger.info(
        "ΔT=%.3e ΔY=%.3e delay=%.3e pv_err=%.3e fitness=%.3e",
        delta_T,
        delta_Y,
        delay_red,
        err,
        fitness,
    )
    return (
        fitness,
        err,
        delay_diff,
        2.0 * max(0.0, size_frac - 0.4),
        reason,
        int(selection.sum()),
        delta_T,
        delta_Y,
        delay_red,
    )


def run_ga_reduction(mech_path, Y0, tf, steps, T0, p0, log_times):
    mech = Mechanism(mech_path)
    full = run_constant_pressure(
        mech.solution,
        T0,
        p0,
        Y0,
        tf,
        nsteps=steps,
        use_mole=False,
        log_times=log_times,
    )

    delay_full, slope_full = ignition_delay(full.time, full.temperature)
    logger.info(
        "Reference ignition delay %.3e s, max dT/dt %.3e K/s",
        delay_full,
        slope_full,
    )

    if (full.temperature[-1] - full.temperature[0]) < 5.0 or (
        np.max(np.abs(full.mass_fractions[-1] - full.mass_fractions[0])) < 1e-6
    ):
        raise RuntimeError("Reference case looks inert (ΔT or ΔY too small). Check T0/φ/tf.")

    genome_len = len(mech.species_names)

    weights = np.ones(genome_len)
    if os.path.exists("data/species_weights.json"):
        with open("data/species_weights.json") as f:
            wmap = json.load(f)
            weights = np.array([wmap.get(s, 1e-3) for s in mech.species_names])
    logger.info(
        "Species weights (first 5): %s", list(zip(mech.species_names, weights))[:5]
    )

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
    k = max(len(critical_idxs), int(0.2 * genome_len))
    seed = np.zeros(genome_len, dtype=int)
    seed[order[-k:]] = 1
    seed[critical_idxs] = 1

    init_pop = np.zeros((pop_size, genome_len), dtype=int)
    for i in range(pop_size):
        individual = seed.copy()
        flips = np.random.choice(genome_len, int(0.2 * genome_len), replace=False)
        individual[flips] ^= 1
        init_pop[i] = individual


    eval_fn = lambda sel: evaluate_selection(
        sel, mech, Y0, tf, full, weights, critical_idxs, T0, p0, steps, log_times
    )

    min_species = len(critical_idxs) + 5
    max_species = max(min_species + 5, int(0.6 * genome_len))

    sel, hist, debug = run_ga(
        genome_len,
        eval_fn,
        GAOptions(
            population_size=pop_size,
            generations=25,
            min_species=min_species,
            max_species=max_species,
            mutation_rate=0.2,
        ),
        return_history=True,
        initial_population=init_pop,
        return_debug=True,
        fixed_indices=critical_idxs,
    )

    os.makedirs("results", exist_ok=True)
    for g, gen in enumerate(debug):
        best_g = max(gen, key=lambda x: x[-2])  # fitness before genome
        pv_err, delay_diff, size_pen, reason, keep_cnt, dT, dY, delay_red, fit, genome = best_g
        logger.info(
            "gen %02d keep=%d ΔT=%.3e delay=%.3e pv_err=%.3e fitness=%.3e",
            g,
            keep_cnt,
            dT,
            delay_red,
            pv_err,
            fit,
        )
        with open(os.path.join("results", f"best_selection_gen{g:02d}.txt"), "w") as f:
            f.write(",".join(map(str, genome.tolist())))

    with open("best_selection.txt", "w") as f:
        f.write(",".join(map(str, sel.tolist())))

    return ["GA"], [sel], [hist], debug, full, weights


def full_pipeline(
    mech_path: str,
    out_dir: str,
    steps: int = 1000,
    tf: float = 0.5,
    phi: float | None = None,
    preset: str = "methane_air",
    T0: float | None = None,
    p0: float | None = None,
    log_times: bool = False,
):
    mech = Mechanism(mech_path)

    if preset == "methane_air":
        phi = phi or 1.0
        x0 = methane_air_mole_fractions(phi)
        Y0 = mole_to_mass_fractions(mech.solution, x0)
    else:
        raise NotImplementedError(preset)

    T0 = T0 or 1500.0
    p0 = p0 or ct.one_atm
    tf = min(tf, 1.0)

    names, sols, hists, debug, full, weights = run_ga_reduction(
        mech_path, Y0, tf, steps, T0, p0, log_times
    )

    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "ga_fitness.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["generation", "fitness"])
        for i, val in enumerate(hists[0]):
            writer.writerow([i, val])

    with open(os.path.join(out_dir, "debug_fitness.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "generation",
            "individual",
            "pv_error",
            "delay_diff",
            "size_penalty",
            "reason",
            "keep_count",
            "delta_T",
            "delta_Y",
            "delay_red",
            "fitness",
        ])
        for g, gen in enumerate(debug):
            for i, d in enumerate(gen):
                pv_err, delay_diff, size_pen, reason, keep_cnt, dT, dY, delay_red, fit, genome = d
                writer.writerow([g, i, pv_err, delay_diff, size_pen, reason, keep_cnt, dT, dY, delay_red, fit])

    best_idx = int(np.argmax(hists[0]))
    best_sel = sols[0]
    keep = [mech.species_names[i] for i, bit in enumerate(best_sel) if bit]

    red_mech = Mechanism(mech_path)
    red_mech.remove_species([s for s in red_mech.species_names if s not in keep])

    with open(os.path.join(out_dir, "reduced_species.txt"), "w") as f:
        for s in keep:
            f.write(s + "\n")

    red = run_constant_pressure(
        red_mech.solution,
        T0,
        p0,
        Y0,
        tf,
        nsteps=len(full.time) - 1,
        use_mole=False,
        log_times=log_times,
    )

    # choose informative species (like the paper)
    key_species = [s for s in ["CH4", "O2", "CO2", "CO", "H2O", "OH"] if s in mech.species_names]

    # 1) Species profiles (full=solid line, reduced=dots)
    plot_profiles(
        full_res=full,
        red_res=red,
        full_names=mech.species_names,
        red_names=red_mech.species_names,
        species=key_species,
        out_base=os.path.join(out_dir, "profiles"),
    )

    # 2) Ignition delay (two bars: full vs reduced)
    delay_full, _ = ignition_delay(full.time, full.temperature)
    delay_red,  _ = ignition_delay(red.time,  red.temperature)
    plot_ignition_delays(
        delays=[delay_full, delay_red],
        labels=["Full", "Reduced"],
        out_base=os.path.join(out_dir, "ignition_delay"),
    )

    # 3) GA convergence
    plot_convergence(hists, names, os.path.join(out_dir, "convergence"))

    # 4) PV error (aligned)
    pv_err = pv_error_aligned(
        full.mass_fractions,
        red.mass_fractions,
        mech.species_names,
        red_mech.species_names,
        weights,
    )
    plot_pv_errors([pv_err], ["GA"], os.path.join(out_dir, "pv_error"))

    print(
        f"\nSummary:\n{'Mechanism':>10} {'Species':>8} {'Reactions':>10} {'Delay[s]':>12} {'PV err':>8}"
    )
    print(
        f"{'Full':>10} {len(mech.species_names):8d} {len(mech.reactions()):10d} {delay_full:12.3e} {0.0:8.3f}"
    )
    print(
        f"{'Reduced':>10} {len(red_mech.species_names):8d} {len(red_mech.reactions()):10d} {delay_red:12.3e} {pv_err:8.3f}"
    )


def plot_profiles(
    full_res: BatchResult,
    red_res: BatchResult,
    full_names: Sequence[str],
    red_names: Sequence[str],
    species: Sequence[str],
    out_base: str,
) -> None:
    # intersection & indices (align by name)
    common = [s for s in species if s in full_names and s in red_names]
    if not common:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No common species to plot", ha="center")
        fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
        fig.savefig(out_base + ".pdf", bbox_inches="tight")
        plt.close(fig)
        return

    idxF = [full_names.index(s) for s in common]
    idxR = [red_names.index(s)  for s in common]

    fig, ax = plt.subplots()
    for i, s in enumerate(common):
        # full: solid line
        line, = ax.semilogx(
            full_res.time,
            full_res.mass_fractions[:, idxF[i]],
            label=f"{s} full",
            linewidth=2,
        )
        # reduced: markers only, same color
        ax.semilogx(
            red_res.time,
            red_res.mass_fractions[:, idxR[i]],
            linestyle="none",
            marker="o",
            markersize=3.0,
            markerfacecolor="none",
            markeredgewidth=1.2,
            color=line.get_color(),
            label=f"{s} reduced",
        )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mass fraction")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)
