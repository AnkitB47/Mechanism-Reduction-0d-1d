import os
import json
import csv
import logging
import numpy as np
import cantera as ct

from typing import Sequence, Dict, Tuple, List

from mechanism.loader import Mechanism
from mechanism.mix import methane_air_mole_fractions, mole_to_mass_fractions
from reactor.batch import BatchResult, run_constant_pressure
from metaheuristics.ga import run_ga, GAOptions
from progress_variable import PV_SPECIES_DEFAULT, pv_error_aligned
from metrics import ignition_delay
from graph.construction import build_species_graph
from gnn.models import train_gnn, predict_scores
from visualizations import (
    plot_ignition_delays,
    plot_convergence,
    plot_pv_errors,
    plot_species_profiles,
    plot_species_residuals,
    plot_progress_variable,
)

logger = logging.getLogger(__name__)


# -------------------------------
# Helpers: label building for GNN
# -------------------------------

def _pv_weights_for(mech: Mechanism) -> np.ndarray:
    """Build PV weights aligned to mech.species_names."""
    W = np.zeros(len(mech.species_names), dtype=float)
    path = "data/species_weights.json"
    if os.path.exists(path):
        with open(path) as f:
            wmap = json.load(f)
        for i, s in enumerate(mech.species_names):
            W[i] = float(wmap.get(s, 0.0))
    else:
        # If file missing, put unit weights for a sensible default PV set
        for i, s in enumerate(mech.species_names):
            if s in {"CO2", "H2O", "CO", "H2", "O", "H", "OH"}:
                W[i] = 1.0
    # Avoid all-zeros
    if not np.any(W):
        W[:] = 1.0
    return W


def _baseline_short_run(
    mech: Mechanism, T0: float, p0: float, Y0: Dict[str, float], tf_short: float, steps_short: int, log_times: bool
) -> BatchResult:
    return run_constant_pressure(
        mech.solution, T0, p0, Y0, tf_short, nsteps=steps_short, use_mole=False, log_times=log_times
    )


def _loo_scores(
    base_mech: Mechanism,
    baseline: BatchResult,
    T0: float,
    p0: float,
    Y0: Dict[str, float],
    tf_short: float,
    steps_short: int,
    log_times: bool,
    critical: Sequence[str],
) -> Dict[str, float]:
    """
    Leave-one-out scores: remove each species (except critical), rerun short simulation,
    compute PV error vs baseline.
    """
    logger.info("Building LOO scores on short window: tf=%.3e, steps=%d", tf_short, steps_short)
    W = _pv_weights_for(base_mech)

    # Convenience lambdas to compute PV error aligned
    def pv_err_against_baseline(red_res: BatchResult, red_mech: Mechanism) -> float:
        return pv_error_aligned(
            baseline.mass_fractions,
            red_res.mass_fractions,
            base_mech.species_names,
            red_mech.species_names,
            W,
        )

    scores: Dict[str, float] = {}
    for s in base_mech.species_names:
        if s in critical:
            scores[s] = 1.0  # force keep
            continue
        try:
            m = Mechanism(base_mech.file_path)
            m.remove_species([s])
            red = _baseline_short_run(m, T0, p0, Y0, tf_short, steps_short, log_times)
            score = pv_err_against_baseline(red, m)
            if not np.isfinite(score):
                score = 0.0
        except Exception as e:
            logger.debug("LOO failed for %s: %s", s, e)
            score = 0.0
        scores[s] = float(score)

    # Normalize to [0,1]
    vals = np.array(list(scores.values()), dtype=float)
    vmax = float(vals.max()) if vals.size else 0.0
    if vmax > 0:
        for k in scores:
            scores[k] = scores[k] / vmax
    return scores


def _centrality_scores(G, mech: Mechanism) -> Dict[str, float]:
    # Simple degree centrality (fast, stable)
    deg = {n: float(G.degree(n)) for n in G.nodes}
    vmax = max(deg.values()) if deg else 1.0
    return {n: (deg[n] / vmax if vmax > 0 else 0.0) for n in mech.species_names}


def _build_species_labels(
    mech: Mechanism,
    G,
    full_short: BatchResult,
    T0: float,
    p0: float,
    Y0: Dict[str, float],
    tf_short: float,
    steps_short: int,
    log_times: bool,
    alpha: float,
    critical: Sequence[str],
    cache_labels: bool,
    cache_dir: str = "results",
) -> Dict[str, float]:
    """
    Blend LOO and centrality to produce labels in [0,1]. Cache to JSON.
    """
    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, "species_labels.json")
    if cache_labels and os.path.exists(cache_path):
        with open(cache_path) as f:
            labels = json.load(f)
        logger.info("Loaded cached labels from %s", cache_path)
        return {k: float(v) for k, v in labels.items()}

    # Compute fresh
    loo = _loo_scores(mech, full_short, T0, p0, Y0, tf_short, steps_short, log_times, critical)
    cen = _centrality_scores(G, mech)

    labels: Dict[str, float] = {}
    for s in mech.species_names:
        labels[s] = float(alpha * loo.get(s, 0.0) + (1.0 - alpha) * cen.get(s, 0.0))

    # Ensure critical are high
    for s in critical:
        labels[s] = 1.0

    with open(cache_path, "w") as f:
        json.dump(labels, f, indent=2)
    logger.info("Wrote labels to %s", cache_path)
    return labels


# -------------------------------
# GA evaluation
# -------------------------------

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

    # hard guard
    if critical_idxs and (selection.sum() < max(4, len(critical_idxs)) or np.any(selection[critical_idxs] == 0)):
        return 1.0, 1.0, 0.0, "missing_critical", int(selection.sum()), 0.0, 0.0, 0.0, -1e6

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
            return 1.0, 1.0, 0.0, str(e), int(selection.sum()), 0.0, 0.0, 0.0, -1e6
        return 1.0, 1.0, 0.0, f"sim_failed:{type(e).__name__}", int(selection.sum()), 0.0, 0.0, 0.0, -1e6
    except Exception as e:
        return 1.0, 1.0, 0.0, f"sim_failed:{type(e).__name__}", int(selection.sum()), 0.0, 0.0, 0.0, -1e6

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

    # penalty if > 40% species removed too aggressively
    size_pen = 2.0 * max(0.0, size_frac - 0.4)
    fitness = -err - 10.0 * delay_diff - size_pen

    logger.info(
        "ΔT=%.3e ΔY=%.3e delay=%.3e pv_err=%.3e fitness=%.3e",
        delta_T,
        delta_Y,
        delay_red,
        err,
        fitness,
    )

    # RETURN ORDER expected by GA debug writer:
    # (pv_err, delay_diff, size_penalty, reason, keep_count, delta_T, delta_Y, delay_red, fitness)
    return (
        float(err),
        float(delay_diff),
        float(size_pen),
        reason,
        int(selection.sum()),
        float(delta_T),
        float(delta_Y),
        float(delay_red),
        float(fitness),
    )


# -------------------------------
# GA driver
# -------------------------------

def run_ga_reduction(
    mech_path: str,
    Y0: Dict[str, float],
    tf: float,
    steps: int,
    T0: float,
    p0: float,
    log_times: bool,
    alpha: float,
    tf_short: float,
    steps_short: int,
    cache_labels: bool,
):
    mech = Mechanism(mech_path)

    # Reference run (full window)
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

    # PV weights for alignment
    weights = _pv_weights_for(mech)
    logger.info("Species weights (first 5): %s", list(zip(mech.species_names, weights))[:5])

    # --- Short-window baseline for LOO
    full_short = _baseline_short_run(mech, T0, p0, Y0, tf_short, steps_short, log_times)

    # Graph + labels
    G = build_species_graph(mech.solution)
    critical = [s for s in ["CH4", "O2", "N2"] if s in mech.species_names]

    labels = _build_species_labels(
        mech, G, full_short, T0, p0, Y0, tf_short, steps_short, log_times, alpha, critical, cache_labels
    )

    # GNN training (labels are already normalized 0..1)
    gnn_model = train_gnn(
        G,
        mech.solution,
        labels,
        epochs=200,
    )

    os.makedirs("results", exist_ok=True)

    scores = predict_scores(
        model=gnn_model,
        G=G,
        solution=mech.solution,
        save_path=os.path.join("results", "gnn_scores.csv"),
    )

    scores_arr = np.array([scores[s] for s in mech.species_names])
    std = float(scores_arr.std())
    logger.info(
        "GNN score stats min=%.3f mean=%.3f max=%.3f std=%.3e",
        float(scores_arr.min()),
        float(scores_arr.mean()),
        float(scores_arr.max()),
        std,
    )

    # Fallback if the GNN is degenerate
    if np.allclose(scores_arr, 0.0) or std < 1e-8:
        logger.warning("Using degree centrality as fallback for seeding.")
        seed_scores = np.array([G.degree(n) for n in mech.species_names], dtype=float)
    else:
        seed_scores = scores_arr

    order = np.argsort(seed_scores)

    # GA seeding & constraints
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

    # Debug drops per-gen
    os.makedirs("results", exist_ok=True)
    for g, gen in enumerate(debug):
        # fitness is at -2 index (last before genome)
        best_g = max(gen, key=lambda x: x[-2])
        pv_err, delay_diff, size_pen, reason, keep_cnt, dT, dY, delay_red, fit, genome = best_g
        logger.info(
            "gen %02d keep=%d ΔT=%.3e delay=%.3e pv_err=%.3e fitness=%.3e",
            g, keep_cnt, dT, delay_red, pv_err, fit,
        )
        with open(os.path.join("results", f"best_selection_gen{g:02d}.txt"), "w") as f:
            f.write(",".join(map(str, genome.tolist())))

    with open("best_selection.txt", "w") as f:
        f.write(",".join(map(str, sel.tolist())))

    return ["GA"], [sel], [hist], debug, full, weights


# -------------------------------
# Full pipeline entry point
# -------------------------------

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
    alpha: float = 0.8,
    tf_short: float | None = None,
    steps_short: int | None = None,
    cache_labels: bool = True,
):
    mech = Mechanism(mech_path)

    # Mixture preset
    if preset == "methane_air":
        phi = phi or 1.0
        x0 = methane_air_mole_fractions(phi)
        Y0 = mole_to_mass_fractions(mech.solution, x0)
    else:
        raise NotImplementedError(preset)

    T0 = T0 or 1500.0
    p0 = p0 or ct.one_atm
    tf = min(tf, 1.0)

    # Defaults for short LOO runs
    if tf_short is None:
        tf_short = 0.25 * tf
    if steps_short is None:
        steps_short = max(50, int(0.25 * steps))

    names, sols, hists, debug, full, weights = run_ga_reduction(
        mech_path, Y0, tf, steps, T0, p0, log_times, alpha, tf_short, steps_short, cache_labels
    )

    os.makedirs(out_dir, exist_ok=True)

    # GA fitness history
    with open(os.path.join(out_dir, "ga_fitness.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["generation", "fitness"])
        for i, val in enumerate(hists[0]):
            writer.writerow([i, val])

    # Per-individual debug
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

    # Build reduced mechanism from best selection
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
    key_species = [s for s in ["O2", "CO2", "H2O", "CO", "CH4", "OH"] if s in mech.species_names]

    # 1) Species profiles (full=solid line, reduced=markers)
    plot_species_profiles(
        full.time,
        full.mass_fractions,
        mech.species_names,
        red.time,
        red.mass_fractions,
        red_mech.species_names,
        key_species,
        os.path.join(out_dir, "profiles"),
    )

    # Residuals (full - reduced)
    plot_species_residuals(
        full.time,
        full.mass_fractions,
        red.mass_fractions,
        mech.species_names,
        red_mech.species_names,
        key_species,
        os.path.join(out_dir, "profiles_residual"),
    )

    # 2) Ignition delay (two bars: full vs reduced)
    delay_full, _ = ignition_delay(full.time, full.temperature)
    delay_red, _ = ignition_delay(red.time, red.temperature)
    plot_ignition_delays(
        delays=[delay_full, delay_red],
        labels=["Full", "Reduced"],
        out_base=os.path.join(out_dir, "ignition_delay"),
    )

    # 3) GA convergence
    plot_convergence(hists, names, os.path.join(out_dir, "convergence"))

    # 4) PV error (aligned) + overlay PV
    pv_err = pv_error_aligned(
        full.mass_fractions,
        red.mass_fractions,
        mech.species_names,
        red_mech.species_names,
        weights,
    )
    plot_pv_errors([pv_err], ["GA"], os.path.join(out_dir, "pv_error"))
    with open(os.path.join(out_dir, "pv_error.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pv_error"])
        writer.writerow([pv_err])

    # PV overlay
    map_full = {s: i for i, s in enumerate(mech.species_names)}
    map_red = {s: i for i, s in enumerate(red_mech.species_names)}
    use = [s for s in PV_SPECIES_DEFAULT if s in map_full and s in map_red]
    w = np.array([weights[map_full[s]] for s in use])
    pv_full = (full.mass_fractions[:, [map_full[s] for s in use]] * w).sum(axis=1)
    pv_red = (red.mass_fractions[:, [map_red[s] for s in use]] * w).sum(axis=1)
    plot_progress_variable(
        full.time,
        pv_full,
        red.time,
        pv_red,
        os.path.join(out_dir, "pv_overlay"),
    )

    print(
        f"\nSummary:\n{'Mechanism':>10} {'Species':>8} {'Reactions':>10} {'Delay[s]':>12} {'PV err':>8}"
    )
    print(
        f"{'Full':>10} {len(mech.species_names):8d} {len(mech.reactions()):10d} {delay_full:12.3e} {0.0:8.3f}"
    )
    print(
        f"{'Reduced':>10} {len(red_mech.species_names):8d} {len(red_mech.reactions()):10d} {delay_red:12.3e} {pv_err:8.3f}"
    )
    if pv_err < 0.05:
        print(
            "Summary:\n"
            f"  Species: {len(mech.species_names)} -> {len(red_mech.species_names)}\n"
            f"  Reactions: {len(mech.reactions())} -> {len(red_mech.reactions())}\n"
            f"  Delay_full/red: {delay_full:.3e} / {delay_red:.3e} s\n"
            f"  PV_error: {pv_err*100:.1f}%"
        )
