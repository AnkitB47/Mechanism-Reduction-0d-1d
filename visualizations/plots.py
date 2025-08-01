from __future__ import annotations

"""Plotting utilities for mechanism reduction results."""

import os
from typing import Sequence, Iterable

import matplotlib.pyplot as plt
import numpy as np


def _save(fig: plt.Figure, out_base: str) -> None:
    """Save figure as PNG and PDF."""
    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    plt.close(fig)


def plot_mole_fraction(
    time_full: np.ndarray,
    y_full: np.ndarray,
    time_red: np.ndarray,
    y_red: np.ndarray,
    species: Sequence[str],
    out_base: str,
) -> None:
    """Plot mole-fraction profiles for selected species."""
    fig, ax = plt.subplots()
    for i, sp in enumerate(species):
        ax.plot(time_full, y_full[:, i], label=f"{sp} full")
        ax.plot(time_red, y_red[:, i], "--", label=f"{sp} red")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mole fraction")
    ax.legend()
    _save(fig, out_base)


def plot_ignition_delays(delays: Sequence[float], labels: Sequence[str], out_base: str) -> None:
    """Bar chart of ignition delays."""
    fig, ax = plt.subplots()
    ax.bar(labels, delays)
    ax.set_ylabel("Ignition delay [s]")
    _save(fig, out_base)


def plot_convergence(histories: Sequence[Iterable[float]], labels: Sequence[str], out_base: str) -> None:
    """Plot convergence histories for metaheuristics."""
    fig, ax = plt.subplots()
    for hist, label in zip(histories, labels):
        if hist:
            ax.plot(list(range(len(hist))), list(hist), label=label)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best fitness")
    ax.legend()
    _save(fig, out_base)


def plot_pv_errors(errors: Sequence[float], labels: Sequence[str], out_base: str) -> None:
    """Bar chart of progress variable errors."""
    fig, ax = plt.subplots()
    ax.bar(labels, errors)
    ax.set_ylabel("PV error")
    _save(fig, out_base)
