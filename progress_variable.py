"""Utility functions for progress variable calculations."""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np
from metrics import pv_error


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
    Y_full: np.ndarray,
    Y_red: np.ndarray,
    names_full: Sequence[str],
    names_red: Sequence[str],
    weights_full: Sequence[float],
) -> float:
    """Compute PV error aligning species by name.

    This utility builds a weight vector for the reduced mechanism that matches
    the ordering in ``Y_red`` while extracting the corresponding mass-fraction
    columns from ``Y_full``.

    Parameters
    ----------
    Y_full, Y_red:
        Mass-fraction matrices for the full and reduced simulations.
    names_full, names_red:
        Species names corresponding to the columns of ``Y_full`` and
        ``Y_red`` respectively.
    weights_full:
        Progress-variable weights defined for the full mechanism ordering.
    """

    name_to_idx = {n: i for i, n in enumerate(names_full)}
    idx_full = [name_to_idx[n] for n in names_red]
    weights = np.array([weights_full[name_to_idx[n]] for n in names_red])
    pv_f = progress_variable(Y_full[:, idx_full], weights)
    pv_r = progress_variable(Y_red, weights)
    return pv_error(pv_f, pv_r)


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

