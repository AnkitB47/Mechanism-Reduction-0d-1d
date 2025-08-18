from __future__ import annotations
"""Plotting utilities for mechanism reduction results."""

import os
from typing import Iterable, Sequence

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


def plot_species_profiles(
    time_full: np.ndarray,
    Y_full: np.ndarray,
    names_full: Sequence[str],
    time_red: np.ndarray,
    Y_red: np.ndarray,
    names_red: Sequence[str],
    species: Sequence[str],
    out_base: str,
) -> None:
    """Overlay species profiles for full vs. reduced mechanisms.

    Full profiles are drawn as solid lines, while reduced-problem profiles are
    plotted using unfilled circle markers. The time axis uses a logarithmic
    scale for clarity.
    """

    common = [s for s in species if s in names_full and s in names_red]
    if not common:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No common species to plot", ha="center")
        _save(fig, out_base)
        return

    idxF = [names_full.index(s) for s in common]
    idxR = [names_red.index(s) for s in common]

    fig, ax = plt.subplots()
    for i, s in enumerate(common):
        line, = ax.semilogx(
            time_full,
            Y_full[:, idxF[i]],
            linewidth=2,
            label=f"{s} full",
        )
        ax.semilogx(
            time_red,
            Y_red[:, idxR[i]],
            linestyle="none",
            marker="o",
            markersize=4,
            fillstyle="none",
            color=line.get_color(),
            markevery=20,
            label=f"{s} reduced",
        )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mass fraction")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, frameon=False)
    _save(fig, out_base)


def plot_species_residuals(
    time: np.ndarray,
    Y_full: np.ndarray,
    Y_red: np.ndarray,
    names_full: Sequence[str],
    names_red: Sequence[str],
    species: Sequence[str],
    out_base: str,
) -> None:
    """Plot residuals ``Y_red - Y_full`` for selected species."""

    common = [s for s in species if s in names_full and s in names_red]
    if not common:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No common species", ha="center")
        _save(fig, out_base)
        return

    idxF = [names_full.index(s) for s in common]
    idxR = [names_red.index(s) for s in common]

    fig, ax = plt.subplots()
    for i, s in enumerate(common):
        ax.semilogx(
            time,
            Y_red[:, idxR[i]] - Y_full[:, idxF[i]],
            label=s,
        )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("ΔY (red - full)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False)
    _save(fig, out_base)


def plot_progress_variable(
    time_full: np.ndarray,
    pv_full: np.ndarray,
    time_red: np.ndarray,
    pv_red: np.ndarray,
    out_base: str,
) -> None:
    """Overlay progress variables for full vs. reduced mechanisms."""

    fig, ax = plt.subplots()
    ax.semilogx(time_full, pv_full, label="PV full", linewidth=2)
    ax.semilogx(
        time_red,
        pv_red,
        linestyle="none",
        marker="o",
        fillstyle="none",
        markevery=20,
        label="PV reduced",
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Progress variable")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False)
    _save(fig, out_base)
