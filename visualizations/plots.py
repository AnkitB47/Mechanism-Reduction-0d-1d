from __future__ import annotations
"""Plotting utilities for mechanism reduction results."""

import os
from typing import Sequence, Iterable

import matplotlib.pyplot as plt
import numpy as np


def _save(fig: plt.Figure, out_base: str) -> None:
    """Save figure as PNG and PDF (pub-ready)."""
    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)


def plot_mole_fraction(
    time_full: np.ndarray,
    y_full: np.ndarray,
    time_red: np.ndarray,
    y_red: np.ndarray,
    species: Sequence[str],
    out_base: str,
) -> None:
    """(Kept for backward-compat; simple overlay)."""
    fig, ax = plt.subplots()
    for i, sp in enumerate(species):
        ax.plot(time_full, y_full[:, i], label=f"{sp} full")
        ax.plot(time_red,  y_red[:,  i], "--", label=f"{sp} red")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mole fraction")
    ax.legend(frameon=False, ncol=2)
    _save(fig, out_base)


def plot_ignition_delays(delays: Sequence[float], labels: Sequence[str], out_base: str) -> None:
    """Bar chart of ignition delays (e.g., ['Full','Reduced'])."""
    fig, ax = plt.subplots()
    bars = ax.bar(labels, delays)
    ax.set_ylabel("Ignition delay [s]")
    ax.grid(axis="y", alpha=0.3)
    # annotate values
    for b, v in zip(bars, delays):
        ax.annotate(f"{v:.2e}", (b.get_x() + b.get_width()/2, b.get_height()),
                    ha="center", va="bottom", fontsize=9, xytext=(0, 3), textcoords="offset points")
    _save(fig, out_base)


def plot_convergence(histories: Sequence[Iterable[float]], labels: Sequence[str], out_base: str) -> None:
    """Plot convergence histories for metaheuristics."""
    fig, ax = plt.subplots()
    for hist, label in zip(histories, labels):
        if hist:
            ax.plot(list(range(len(hist))), list(hist), label=label, linewidth=2)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Best fitness")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    _save(fig, out_base)


def plot_pv_errors(errors: Sequence[float], labels: Sequence[str], out_base: str) -> None:
    """Bar chart of progress variable errors."""
    fig, ax = plt.subplots()
    ax.bar(labels, errors)
    ax.set_ylabel("PV error")
    ax.grid(axis="y", alpha=0.3)
    _save(fig, out_base)
