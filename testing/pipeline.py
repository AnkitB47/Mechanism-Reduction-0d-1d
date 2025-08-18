import os
import json
import csv
import logging
import numpy as np
import cantera as ct

from typing import Sequence, Dict, Tuple, List

from mechanism.loader import Mechanism
from mechanism.mix import methane_air_mole_fractions, mole_to_mass_fractions, HR_PRESETS
from reactor.batch import BatchResult, run_constant_pressure, run_isothermal_const_p
from metaheuristics.ga import run_ga, GAOptions
from progress_variable import PV_SPECIES_DEFAULT, pv_error_aligned
from metrics import ignition_delay
from graph.construction import build_species_graph
from gnn.models import train_gnn, predict_scores
from visualizations import (
    plot_ignition_delays,
    plot_convergence,
    plot_species_profiles,
    plot_species_residuals,
    plot_progress_variable,
    plot_timescales,
)
from timescales import pv_timescale, spts

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
    mech: Mechanism,
    T0: float,
    p0: float,
    Y0: Dict[str, float],
    tf_short: float,
    steps_short: int,
    log_times: bool,
    runner,
) -> BatchResult:
    return runner(mech.solution, T0, p0, Y0, tf_short, nsteps=steps_short, use_mole=False, log_times=log_times)


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
    runner,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    Leave-one-out scores and post-ignition residual integrals for each species.
    Returns two dictionaries keyed by species name.
    """
    logger.info("Building LOO scores on short window: tf=%.3e, steps=%d", tf_short, steps_short)
    W = _pv_weights_for(base_mech)

    tau_short = baseline.ignition_delay or ignition_delay(baseline.time, baseline.temperature)[0]
    mask_post = baseline.time > tau_short

    scores: Dict[str, float] = {}
    resid: Dict[str, float] = {}
    for s in base_mech.species_names:
        if s in critical:
            scores[s] = 1.0
            resid[s] = 0.0
            continue
        try:
            m = Mechanism(base_mech.file_path)
            m.remove_species([s])
            red = _baseline_short_run(m, T0, p0, Y0, tf_short, steps_short, log_times, runner)

            score = pv_error_aligned(
                baseline.mass_fractions,
                red.mass_fractions,
                base_mech.species_names,
                m.species_names,
                W,
            )

            Y_red_interp = np.zeros_like(baseline.mass_fractions)
            for j, name in enumerate(base_mech.species_names):
                if name in m.species_names:
                    jr = m.species_names.index(name)
                    Y_red_interp[:, j] = np.interp(baseline.time, red.time, red.mass_fractions[:, jr])

            diff = np.abs(baseline.mass_fractions[mask_post] - Y_red_interp[mask_post])
            resid_val = float(
                np.trapz((diff * W[None, :]).sum(axis=1), baseline.time[mask_post])
            )

            if not np.isfinite(score):
                score = 0.0
            if not np.isfinite(resid_val):
                resid_val = 0.0
        except Exception as e:
            logger.debug("LOO failed for %s: %s", s, e)
            score = 0.0
            resid_val = 0.0
        scores[s] = float(score)
        resid[s] = resid_val

    # Normalize to [0,1]
    vals = np.array(list(scores.values()), dtype=float)
    vmax = float(vals.max()) if vals.size else 0.0
    if vmax > 0:
        for k in scores:
            scores[k] = scores[k] / vmax

    rv = np.array(list(resid.values()), dtype=float)
    rmax = float(rv.max()) if rv.size else 0.0
    if rmax > 0:
        for k in resid:
            resid[k] = resid[k] / rmax

    return scores, resid


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
    runner,
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
    loo, resid = _loo_scores(
        mech, full_short, T0, p0, Y0, tf_short, steps_short, log_times, critical, runner
    )
    cen = _centrality_scores(G, mech)

    labels: Dict[str, float] = {}
    for s in mech.species_names:
        base = alpha * loo.get(s, 0.0) + (1.0 - alpha) * cen.get(s, 0.0)
        labels[s] = float(base * (1.0 + 0.5 * resid.get(s, 0.0)))

    # Renormalize to [0,1]
    vals = np.array(list(labels.values()), dtype=float)
    vmax = float(vals.max()) if vals.size else 0.0
    if vmax > 0:
        for k in labels:
            labels[k] = labels[k] / vmax

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
    runner,
    tau_full,
    tau_pv_full,
    tau_spts_full,
    target_species: int | None,
    zeta: float = 20.0,
):
    reason = ""

    if critical_idxs and (selection.sum() < max(4, len(critical_idxs)) or np.any(selection[critical_idxs] == 0)):
        info = (1.0, 1.0, 0.0, 0.0, 0.0, "missing_critical", int(selection.sum()), 0.0, 0.0, 0.0)
        return -1e6, *info

    sel = selection.copy()
    attempt = 0
    while True:
        keep = [base_mech.species_names[i] for i, bit in enumerate(sel) if bit]
        mech = Mechanism(base_mech.file_path)
        remove = [s for s in mech.species_names if s not in keep]

        try:
            if remove:
                mech.remove_species(remove)
            res = runner(
                mech.solution,
                T0,
                p0,
                Y0,
                tf,
                nsteps=steps,
                use_mole=False,
                log_times=log_times,
                time_grid=full_res.time,
            )
        except RuntimeError as e:
            if "no_reaction" in str(e) and attempt == 0:
                sel = sel.copy()
                sel[critical_idxs] = 1
                attempt += 1
                continue
            info = (1.0, 1.0, 0.0, 0.0, 0.0, str(e), int(sel.sum()), 0.0, 0.0, 0.0)
            return -1e6, *info
        except Exception as e:
            info = (1.0, 1.0, 0.0, 0.0, 0.0, f"sim_failed:{type(e).__name__}", int(sel.sum()), 0.0, 0.0, 0.0)
            return -1e6, *info
        break

    # interpolate reduced onto full time grid for penalties
    Y_red_interp = np.zeros((len(full_res.time), len(mech.species_names)))
    for j, s in enumerate(mech.species_names):
        Y_red_interp[:, j] = np.interp(full_res.time, res.time, res.mass_fractions[:, j])

    err = pv_error_aligned(
        full_res.mass_fractions,
        Y_red_interp,
        base_mech.species_names,
        mech.species_names,
        np.asarray(weights),
    )
    delay_red, _ = ignition_delay(res.time, res.temperature)
    delay_diff = abs(delay_red - tau_full) / max(tau_full, 1e-12)
    delta_T = float(res.temperature[-1] - res.temperature[0])
    delta_Y = float(np.max(np.abs(res.mass_fractions[-1] - res.mass_fractions[0])))

    # timescale mismatch
    _, tau_pv_red = pv_timescale(res.time, res.mass_fractions, mech.species_names)
    tau_spts_red = spts(res.time, res.mass_fractions)
    log_full_pv = np.log10(tau_pv_full + 1e-30)
    log_red_pv = np.log10(tau_pv_red + 1e-30)
    log_full_spts = np.log10(tau_spts_full + 1e-30)
    log_red_spts = np.log10(tau_spts_red + 1e-30)
    tau_mis = np.linalg.norm(log_full_pv - log_red_pv) + np.linalg.norm(log_full_spts - log_red_spts)

    # post-ignition species penalty
    map_full = {s: i for i, s in enumerate(base_mech.species_names)}
    map_red = {s: i for i, s in enumerate(mech.species_names)}
    Y_red_full = np.zeros_like(full_res.mass_fractions)
    for s, j in map_full.items():
        if s in map_red:
            Y_red_full[:, j] = Y_red_interp[:, map_red[s]]
    mask_post = full_res.time > tau_full
    diff = np.abs(full_res.mass_fractions[mask_post] - Y_red_full[mask_post])
    pen_species = float(np.sum(np.asarray(weights) * diff.mean(axis=0)))

    keep_cnt = int(sel.sum())
    size_frac = keep_cnt / len(selection)
    size_pen = 3.0 * max(0.0, size_frac - 0.4)
    if target_species is not None:
        qpen = ((keep_cnt - target_species) / len(selection)) ** 2
        size_pen += 3.0 * qpen

    fitness = -(1.0 * err + 12.0 * delay_diff + 1.0 * tau_mis + zeta * pen_species + size_pen)

    info = (
        float(err),
        float(delay_diff),
        float(size_pen),
        float(tau_mis),
        float(pen_species),
        reason,
        keep_cnt,
        float(delta_T),
        float(delta_Y),
        float(delay_red),
    )
    logger.info(
        "ΔT=%.3e ΔY=%.3e delay=%.3e pv_err=%.3e τ_mis=%.3e pen=%.3e fitness=%.3e",
        delta_T,
        delta_Y,
        delay_red,
        err,
        tau_mis,
        pen_species,
        fitness,
    )

    return float(fitness), *info


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
    isothermal: bool,
    min_species: int | None = None,
    max_species: int | None = None,
    target_species: int | None = None,
    generations: int = 60,
    population_size: int = 40,
    mutation_rate: float = 0.25,
):
    mech = Mechanism(mech_path)

    runner = run_isothermal_const_p if isothermal else run_constant_pressure

    # Reference run (full window)
    full = runner(
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
    if log_times:
        t0 = max(1e-12, 0.8 * delay_full)
        t1 = min(tf, 1.6 * delay_full)
        dense = np.linspace(t0, t1, 200)
        coarse = np.geomspace(1e-12, tf, steps)
        time_grid = np.unique(np.concatenate((coarse, dense)))
        time_grid = np.insert(time_grid, 0, 0.0)
        full = runner(
            mech.solution,
            T0,
            p0,
            Y0,
            tf,
            nsteps=len(time_grid) - 1,
            use_mole=False,
            log_times=False,
            time_grid=time_grid,
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

    pv_full, tau_pv_full = pv_timescale(full.time, full.mass_fractions, mech.species_names)
    tau_spts_full = spts(full.time, full.mass_fractions)

    # PV weights for alignment
    weights = _pv_weights_for(mech)
    logger.info("Species weights (first 5): %s", list(zip(mech.species_names, weights))[:5])

    # --- Short-window baseline for LOO
    full_short = _baseline_short_run(mech, T0, p0, Y0, tf_short, steps_short, log_times, runner)

    # Graph + labels
    G = build_species_graph(mech.solution)
    critical = [
        s
        for s in [
            "CH4",
            "O2",
            "N2",
            "H2O",
            "CO",
            "CO2",
            "OH",
            "H",
            "O",
            "HO2",
            "H2O2",
            "H2",
        ]
        if s in mech.species_names
    ]

    labels = _build_species_labels(
        mech,
        G,
        full_short,
        T0,
        p0,
        Y0,
        tf_short,
        steps_short,
        log_times,
        alpha,
        critical,
        cache_labels,
        runner,
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
    if std < 1e-6:
        logger.warning("Using degree centrality as fallback for seeding.")
        seed_scores = np.array([G.degree(n) for n in mech.species_names], dtype=float)
    else:
        seed_scores = scores_arr

    order = np.argsort(seed_scores)

    # GA seeding & constraints
    critical_idxs = [mech.species_names.index(s) for s in critical]
    pop_size = population_size
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
        sel,
        mech,
        Y0,
        tf,
        full,
        weights,
        critical_idxs,
        T0,
        p0,
        len(full.time) - 1,
        log_times,
        runner,
        delay_full,
        tau_pv_full,
        tau_spts_full,
        target_species,
    )

    ms = min_species or (len(critical_idxs) + 5)
    M = max_species or max(ms + 5, int(0.6 * genome_len))

    sel, hist, debug = run_ga(
        genome_len,
        eval_fn,
        GAOptions(
            population_size=pop_size,
            generations=generations,
            min_species=ms,
            max_species=M,
            mutation_rate=mutation_rate,
        ),
        return_history=True,
        initial_population=init_pop,
        return_debug=True,
        fixed_indices=critical_idxs,
    )

    # Debug drops per-gen
    os.makedirs("results", exist_ok=True)
    for g, gen in enumerate(debug):
        best_g = max(gen, key=lambda x: x[-2])
        pv_err, delay_diff, size_pen, tau_mis, pen_species, reason, keep_cnt, dT, dY, delay_red, fit, genome = best_g
        logger.info(
            "gen %02d keep=%d ΔT=%.3e delay=%.3e pv_err=%.3e τ_mis=%.3e pen=%.3e fitness=%.3e",
            g,
            keep_cnt,
            dT,
            delay_red,
            pv_err,
            tau_mis,
            pen_species,
            fit,
        )
        with open(os.path.join("results", f"best_selection_gen{g:02d}.txt"), "w") as f:
            f.write(",".join(map(str, genome.tolist())))

    # Selection report
    kept = []
    dropped = []
    for i, s in enumerate(mech.species_names):
        if sel[i]:
            kept.append((s, labels.get(s, 0.0), scores.get(s, 0.0)))
        else:
            dropped.append((s, labels.get(s, 0.0), scores.get(s, 0.0)))
    kept.sort(key=lambda x: x[1], reverse=True)
    dropped.sort(key=lambda x: x[1], reverse=True)
    with open(os.path.join("results", "selection_report.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["status", "species", "label", "score"])
        for status, data in [("kept", kept[:20]), ("dropped", dropped[:20])]:
            for s, lab, sc in data:
                writer.writerow([status, s, lab, sc])

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
    isothermal: bool = False,
    min_species: int | None = None,
    max_species: int | None = None,
    target_species: int | None = None,
    generations: int = 60,
    population: int = 40,
    mutation: float = 0.25,
    focus: str = "auto",
    focus_window: Tuple[float, float] | None = None,
):
    mech = Mechanism(mech_path)

    # Mixture preset
    if preset == "methane_air":
        phi = phi or 1.0
        x0 = methane_air_mole_fractions(phi)
        Y0 = mole_to_mass_fractions(mech.solution, x0)
        T0 = T0 or 1500.0
        p0 = p0 or ct.one_atm
    elif preset in HR_PRESETS:
        T0, p0, Y0 = HR_PRESETS[preset](mech.solution)
    else:
        raise NotImplementedError(preset)

    tf = min(tf, 1.0)

    # Defaults for short LOO runs
    if tf_short is None:
        tf_short = 0.25 * tf
    if steps_short is None:
        steps_short = max(50, int(0.25 * steps))

    use_iso = isothermal and preset in HR_PRESETS
    runner = run_isothermal_const_p if use_iso else run_constant_pressure

    names, sols, hists, debug, full, weights = run_ga_reduction(
        mech_path,
        Y0,
        tf,
        steps,
        T0,
        p0,
        log_times,
        alpha,
        tf_short,
        steps_short,
        cache_labels,
        use_iso,
        min_species=min_species,
        max_species=max_species,
        target_species=target_species,
        generations=generations,
        population_size=population,
        mutation_rate=mutation,
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
            "tau_mismatch",
            "pen_species",
            "reason",
            "keep_count",
            "delta_T",
            "delta_Y",
            "delay_red",
            "fitness",
        ])
        for g, gen in enumerate(debug):
            for i, d in enumerate(gen):
                pv_err, delay_diff, size_pen, tau_mis, pen_species, reason, keep_cnt, dT, dY, delay_red, fit, genome = d
                writer.writerow([g, i, pv_err, delay_diff, size_pen, tau_mis, pen_species, reason, keep_cnt, dT, dY, delay_red, fit])

        # Build reduced mechanism from best selection
    best_sel = sols[0]
    keep = [mech.species_names[i] for i, bit in enumerate(best_sel) if bit]

    red_mech = Mechanism(mech_path)
    red_mech.remove_species([s for s in red_mech.species_names if s not in keep])

    with open(os.path.join(out_dir, "reduced_species.txt"), "w") as f:
        for s in keep:
            f.write(s + "\n")

    # Re-run reduced on same window/grid choice
    red = runner(
        red_mech.solution,
        T0,
        p0,
        Y0,
        tf,
        nsteps=len(full.time) - 1,
        use_mole=False,
        log_times=log_times,
        time_grid=full.time,
    )

    # Informative species (like paper)
    key_species = [s for s in ["O2", "CO2", "H2O", "CO", "CH4", "OH"] if s in mech.species_names]

    # --- ignition delays (both)
    delay_full, _ = ignition_delay(full.time, full.temperature)
    delay_red,  _ = ignition_delay(red.time,  red.temperature)

    # Residuals per species after ignition
    mask = full.time > delay_full
    map_full = {s: i for i, s in enumerate(mech.species_names)}
    map_red = {s: i for i, s in enumerate(red_mech.species_names)}
    with open(os.path.join(out_dir, "residual_species.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "mean_res", "max_res"])
        for s in mech.species_names:
            if s not in map_red:
                continue
            yf = full.mass_fractions[mask, map_full[s]]
            yr = np.interp(full.time[mask], red.time, red.mass_fractions[:, map_red[s]])
            diff = np.abs(yf - yr)
            writer.writerow([s, float(diff.mean()), float(diff.max())])

    if focus == "auto":
        fw = (0.8 * delay_full, 1.6 * delay_full)
    elif focus == "window" and focus_window is not None:
        fw = focus_window
    else:
        fw = None

    plot_species_profiles(
        full.time,
        full.mass_fractions,
        mech.species_names,
        red.time,
        red.mass_fractions,
        red_mech.species_names,
        key_species,
        os.path.join(out_dir, "profiles"),
        tau_full=delay_full,
        tau_red=delay_red,
        focus=focus,
        focus_window=fw,
    )

    plot_species_residuals(
        full.time,
        full.mass_fractions,
        red.time,
        red.mass_fractions,
        mech.species_names,
        red_mech.species_names,
        key_species,
        os.path.join(out_dir, "profiles_residual"),
        tau_full=delay_full,
        focus=focus,
        focus_window=fw,
    )

    # 2) Ignition delay bars
    plot_ignition_delays(
        delays=[delay_full, delay_red],
        labels=["Full", "Reduced"],
        out_base=os.path.join(out_dir, "ignition_delay"),
    )

    # 3) GA convergence
    plot_convergence(hists, names, os.path.join(out_dir, "convergence"))

    # 4) PV error (aligned) + PV overlay
    pv_err = pv_error_aligned(
        full.mass_fractions,
        red.mass_fractions,
        mech.species_names,
        red_mech.species_names,
        weights,
    )
    with open(os.path.join(out_dir, "pv_error.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pv_error"])
        writer.writerow([pv_err])

    # PV overlay (compute BEFORE plotting)
    map_full = {s: i for i, s in enumerate(mech.species_names)}
    map_red  = {s: i for i, s in enumerate(red_mech.species_names)}
    use = [s for s in PV_SPECIES_DEFAULT if s in map_full and s in map_red]
    w = np.array([weights[map_full[s]] for s in use])
    pv_full = (full.mass_fractions[:, [map_full[s] for s in use]] * w).sum(axis=1)
    pv_red  = (red.mass_fractions[:,  [map_red[s]  for s in use]] * w).sum(axis=1)

    plot_progress_variable(
        full.time,
        pv_full,
        red.time,
        pv_red,
        os.path.join(out_dir, "pv_overlay"),
        tau_full=delay_full,
        focus=focus,
        focus_window=fw,
    )

    # 5) Time-scales overlay (PVTS / SPTS)
    _, tau_pv_full = pv_timescale(full.time, full.mass_fractions, mech.species_names)
    _, tau_pv_red  = pv_timescale(red.time,  red.mass_fractions,  red_mech.species_names)
    tau_spts_full  = spts(full.time, full.mass_fractions)
    tau_spts_red   = spts(red.time,  red.mass_fractions)

    plot_timescales(
        full.time,
        tau_pv_full,
        tau_spts_full,
        red.time,
        tau_pv_red,
        tau_spts_red,
        os.path.join(out_dir, "timescales"),
    )

    with open(os.path.join(out_dir, "timescales.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time","tau_pv_full","tau_pv_red","tau_spts_full","tau_spts_red"])
        for t, tpvf, tpvr, tsf, tsr in zip(full.time, tau_pv_full, tau_pv_red, tau_spts_full, tau_spts_red):
            writer.writerow([t, tpvf, tpvr, tsf, tsr])

    # --- Console summary
    print(f"\nSummary:\n{'Mechanism':>10} {'Species':>8} {'Reactions':>10} {'Delay[s]':>12} {'PV err':>8}")
    print(f"{'Full':>10} {len(mech.species_names):8d} {len(mech.reactions()):10d} {delay_full:12.3e} {0.0:8.3f}")
    print(f"{'Reduced':>10} {len(red_mech.species_names):8d} {len(red_mech.reactions()):10d} {delay_red:12.3e} {pv_err:8.3f}")
    if pv_err < 0.05:
        print("Summary:\n"
              f"  Species: {len(mech.species_names)} -> {len(red_mech.species_names)}\n"
              f"  Reactions: {len(mech.reactions())} -> {len(red_mech.reactions())}\n"
              f"  Delay_full/red: {delay_full:.3e} / {delay_red:.3e} s\n"
              f"  PV_error: {pv_err*100:.1f}%")
