"""Utility functions for progress variable calculations."""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np
from metrics import pv_error

PV_SPECIES_DEFAULT = ["CO2", "H2O", "CO", "H2", "O", "H", "OH"]


def progress_variable(Y: np.ndarray, weights: Sequence[float]) -> np.ndarray:
    """Return :math:`PV(t)=\sum_i w_i Y_i(t)`.

    Parameters
    ----------
    Y:
        Matrix of mass fractions with shape ``(n_steps, n_species)``.
    weights:
        Sequence of species weights ``w_i``.

    Returns
    -------
    np.ndarray
        Progress variable as a one-dimensional array of length ``n_steps``.
    """

    W = np.asarray(weights, dtype=float)
    if Y.shape[1] != W.size:
        raise ValueError("Weight length must match number of species")
    if not np.any(W):
        logging.getLogger(__name__).warning("All PV species weights are zero; returning zeros")
        return np.zeros(Y.shape[0])
    pv = (Y * W[None, :]).sum(axis=1)
    if pv[-1] <= pv[0]:
        pv[-1] = pv[0] + 1e-6
    return pv


def pv_error_aligned(
    full_Y: np.ndarray,
    red_Y: np.ndarray,
    full_names: Sequence[str],
    red_names: Sequence[str],
    weights: Sequence[float] | None = None,
) -> float:
    """Compute PV error aligning species by name.

    Parameters
    ----------
    full_Y, red_Y:
        Mass-fraction matrices for the full and reduced simulations.
    full_names, red_names:
        Species names corresponding to the columns of ``full_Y`` and
        ``red_Y`` respectively.
    weights:
        Optional weights for the PV species.
    """

    logger = logging.getLogger(__name__)
    map_full = {s: i for i, s in enumerate(full_names)}
    map_red = {s: i for i, s in enumerate(red_names)}
    use = [s for s in PV_SPECIES_DEFAULT if s in map_full and s in map_red]
    if not use:
        return 1e9

    F = full_Y[:, [map_full[s] for s in use]]
    R = red_Y[:, [map_red[s] for s in use]]
    if weights is not None:
        w = np.asarray([weights[map_full[s]] for s in use])
    else:
        w = np.ones(len(use))

    logger.info("PV species/weights (first 5): %s", list(zip(use, w))[:5])
    logger.info(
        "PV alignment (first 5): %s",
        [f"{s}->{s}" for s in use[:5]],
    )

    pvF = (F * w).sum(axis=1)
    pvR = (R * w).sum(axis=1)
    err = np.linalg.norm(pvF - pvR) / (np.linalg.norm(pvF) + 1e-16)
    return float(err)


def optimise_weights(
    Y: np.ndarray,
    pv_ref: np.ndarray,
    eval_options: "GAOptions | None" = None,
) -> np.ndarray:
    """Optimise progress-variable weights via a genetic algorithm.

    The optimisation minimises the L2 error between ``progress_variable(Y, w)``
    and a reference progress variable ``pv_ref``.

    Parameters
    ----------
    Y:
        Mass-fraction matrix ``(n_steps, n_species)``.
    pv_ref:
        Reference progress variable of length ``n_steps``.
    eval_options:
        Optional :class:`~metaheuristics.ga.GAOptions` to control the GA.

    Returns
    -------
    np.ndarray
        Optimised weight vector ``w``.
    """

    from metaheuristics.ga import GAOptions, run_ga

    options = eval_options or GAOptions()

    genome_length = Y.shape[1]

    def eval_fn(genome: np.ndarray) -> float:
        # map binary genome to weights in [0, 1]
        w = genome.astype(float)
        pv = progress_variable(Y, w)
        err = np.linalg.norm(pv - pv_ref) / (np.linalg.norm(pv_ref) + 1e-12)
        # GA maximises fitness, so return negative error
        return -err

    best = run_ga(genome_length, eval_fn, options)
    return best.astype(float)

