from __future__ import annotations
"""Plotting utilities for mechanism reduction results (paper-ready)."""

import os
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np


# -------------------------
# small helpers
# -------------------------
def _save(fig: plt.Figure, out_base: str) -> None:
    """Save figure as PNG and PDF (pub-ready)."""
    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)


def _downsample_markevery(n: int, target: int = 30) -> int:
    """Pick a sensible 'markevery' so ~target markers appear."""
    if n <= target:
        return 1
    step = int(np.ceil(n / target))
    return max(step, 2)


# -------------------------
# legacy/simple plots
# -------------------------
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


# -------------------------
# paper-style profiles
# -------------------------
# (fixed colors per species so legend पढ़ते ही पहचान बने)
_COLOR = {
    "O2":  "#1f77b4",  # blue
    "CO2": "#2ca02c",  # green
    "H2O": "#d62728",  # red
    "CO":  "#17becf",  # teal
    "CH4": "#9467bd",  # purple
    "OH":  "#8c564b",  # brown
}

def _style_color(name: str, fallback: str) -> str:
    return _COLOR.get(name, fallback)


def plot_species_profiles(
    time_full: np.ndarray,
    Y_full: np.ndarray,
    names_full: Sequence[str],
    time_red: np.ndarray,
    Y_red: np.ndarray,
    names_red: Sequence[str],
    species: Sequence[str],
    out_base: str,
    *,
    tau_full: float | None = None,
    tau_red: float | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    markevery: int = 20,
) -> None:
    """
    Paper-style overlay of species profiles.
    - full = solid line
    - reduced = hollow circle markers (same color)
    - log x-axis; optional zoom around ignition delay `tau_full`
    - optional y-limits for a fixed vertical range (e.g., (-0.02, 0.45))
    """
    # common species in the requested order
    common = [s for s in species if s in names_full and s in names_red]
    if not common:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No common species to plot", ha="center")
        _save(fig, out_base)
        return

    idxF = [names_full.index(s) for s in common]
    idxR = [names_red.index(s) for s in common]

    fig, ax = plt.subplots()

    if tau_full is not None and np.isfinite(tau_full) and tau_full > 0:
        ax.axvline(tau_full, color="0.2", ls="--", lw=1, zorder=0)
    if tau_red is not None and np.isfinite(tau_red) and tau_red > 0:
        ax.axvline(tau_red, color="0.6", ls="--", lw=1, zorder=0)

    mk_every = markevery
    for i, s in enumerate(common):
        color = _style_color(s, fallback=f"C{i}")
        ax.semilogx(time_full, Y_full[:, idxF[i]], linewidth=2.2, color=color, label=f"{s} full")
        ax.semilogx(
            time_red,
            Y_red[:, idxR[i]],
            linestyle="none",
            marker="o",
            markersize=4.0,
            fillstyle="none",
            color=color,
            markevery=mk_every,
            label=f"{s} reduced",
        )

    if xlim is not None:
        ax.set_xlim(*xlim)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mass fraction")
    ax.grid(True, which="both", alpha=0.3)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.legend(ncol=2, frameon=False)
    _save(fig, out_base)


def plot_species_residuals(
    time_full: np.ndarray,
    Y_full: np.ndarray,
    time_red: np.ndarray,
    Y_red: np.ndarray,
    names_full: Sequence[str],
    names_red: Sequence[str],
    species: Sequence[str],
    out_base: str,
    xlim: tuple[float, float] | None = None,
) -> None:
    """Plot residuals ``Y_full(t) - Y_red(t)`` on full time grid (reduced is interpolated)."""
    common = [s for s in species if s in names_full and s in names_red]
    if not common:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No common species", ha="center")
        _save(fig, out_base)
        return

    # interpolate reduced to full grid so ΔY aligns
    Y_red_interp = np.empty((len(time_full), len(common)))
    for j, s in enumerate(common):
        jF = names_full.index(s)
        jR = names_red.index(s)
        Y_red_interp[:, j] = np.interp(time_full, time_red, Y_red[:, jR])

    fig, ax = plt.subplots()
    for j, s in enumerate(common):
        color = _style_color(s, fallback=f"C{j}")
        ax.semilogx(time_full, Y_full[:, names_full.index(s)] - Y_red_interp[:, j], label=s, color=color)

    ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("ΔY (full - red)")
    ax.grid(True, which="both", alpha=0.3)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.legend(frameon=False, ncol=2)
    _save(fig, out_base)


def plot_progress_variable(
    time_full: np.ndarray,
    pv_full: np.ndarray,
    time_red: np.ndarray,
    pv_red: np.ndarray,
    out_base: str,
    *,
    tau: float | None = None,
    xlim: tuple[float, float] | None = None,
) -> None:
    """Overlay progress variables for full vs. reduced mechanisms."""
    fig, ax = plt.subplots()
    ax.semilogx(time_full, pv_full, label="PV full", linewidth=2)
    mk_every = _downsample_markevery(len(time_red))
    ax.semilogx(time_red, pv_red, linestyle="none", marker="o",
                fillstyle="none", markevery=mk_every, label="PV reduced")
    if tau is not None and np.isfinite(tau) and tau > 0:
        ax.axvline(tau, color="0.4", ls="--", lw=1)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Progress variable")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False)
    _save(fig, out_base)


def plot_timescales(
    time_full: np.ndarray,
    tau_pv_full: np.ndarray,
    tau_spts_full: np.ndarray,
    time_red: np.ndarray,
    tau_pv_red: np.ndarray,
    tau_spts_red: np.ndarray,
    out_base: str,
) -> None:
    """Overlay PVTS and SPTS time scales for full and reduced cases."""
    fig, ax = plt.subplots()
    mk_every = _downsample_markevery(len(time_red))
    ax.semilogx(time_full, tau_pv_full, label="PVTS full", linewidth=2)
    ax.semilogx(time_red, tau_pv_red, linestyle="--", marker="o",
                markevery=mk_every, fillstyle="none", label="PVTS reduced")
    ax.semilogx(time_full, tau_spts_full, label="SPTS full", linewidth=2)
    ax.semilogx(time_red, tau_spts_red, linestyle="--", marker="s",
                markevery=mk_every, fillstyle="none", label="SPTS reduced")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Time scale [s]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, ncol=2)
    _save(fig, out_base)
