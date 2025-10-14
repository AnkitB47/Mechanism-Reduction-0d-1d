#!/usr/bin/env python3
"""
Quick baseline vs reduced comparison helper.

Given two result roots (each containing per-case folders with
`axial_profiles.csv`, `unified_profile.csv`, `diagnostics.json`, `kpis.json`,
etc.), this script builds a KPI summary table and generates a couple of
overlay plots per case (temperature and dry-basis species).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Any, Iterable, List

import matplotlib.pyplot as plt
import pandas as pd


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing JSON: {path}")
    with path.open("r") as fh:
        return json.load(fh)


def _load_axial(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing axial CSV: {path}")
    df = pd.read_csv(path)
    # Normalise to expected columns for plotting
    if "z_m" not in df.columns:
        raise ValueError(f"'z_m' column not found in {path}")
    if "T_C" not in df.columns and "T_K" in df.columns:
        df["T_C"] = df["T_K"] - 273.15
    return df


def _dry_species_from_kpi(data: Dict[str, Any]) -> Dict[str, float]:
    sim_dry = data.get("sim_dry") or {}
    out = {k: float(sim_dry.get(k, 0.0)) for k in ("H2", "CO", "CO2", "CH4")}
    out["H2_CO"] = (
        out["H2"] / max(out["CO"], 1e-12) if out["CO"] > 0 else 0.0
    )
    out["T_out_K"] = float(data.get("pfr_outlet_temperature_C", 0.0)) + 273.15
    return out


def _percent_diff(ref: float, test: float) -> float:
    if abs(ref) < 1e-12:
        return float("nan")
    return 100.0 * (test - ref) / ref


def _make_case_plot(
    case: str,
    base_df: pd.DataFrame,
    red_df: pd.DataFrame,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    def plot_temperature():
        if "T_C" not in base_df.columns or "T_C" not in red_df.columns:
            return
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(base_df["z_m"], base_df["T_C"], label="Baseline", lw=2)
        ax.plot(red_df["z_m"], red_df["T_C"], label="Reduced", lw=2, ls="--")
        ax.set_xlabel("Axial position z (m)")
        ax.set_ylabel("Temperature (°C)")
        ax.set_title(f"{case}: Temperature profile")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "T_vs_z.png", dpi=180)
        plt.close(fig)

    def plot_species_vs_z():
        species = ["X_H2", "X_CO", "X_CO2", "X_CH4"]
        available = [s for s in species if s in base_df.columns and s in red_df.columns]
        if not available:
            return
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for col in available:
            ax.plot(base_df["z_m"], base_df[col], label=f"{col} (Baseline)", lw=2)
            ax.plot(
                red_df["z_m"],
                red_df[col],
                label=f"{col} (Reduced)",
                lw=2,
                ls="--",
            )
        ax.set_xlabel("Axial position z (m)")
        ax.set_ylabel("Dry mole fraction")
        ax.set_title(f"{case}: Dry species profile")
        ax.grid(True, alpha=0.3)
        ax.legend(ncol=2, fontsize="small")
        fig.tight_layout()
        fig.savefig(out_dir / "species_vs_z.png", dpi=180)
        plt.close(fig)

    def plot_species_vs_time():
        time_col = "residence_time_s"
        if time_col not in base_df.columns or time_col not in red_df.columns:
            return
        species = ["X_H2", "X_CO", "X_CO2", "X_CH4"]
        available = [s for s in species if s in base_df.columns and s in red_df.columns]
        if not available:
            return
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for col in available:
            ax.plot(base_df[time_col], base_df[col], label=f"{col} (Baseline)", lw=2)
            ax.plot(
                red_df[time_col],
                red_df[col],
                label=f"{col} (Reduced)",
                lw=2,
                ls="--",
            )
        ax.set_xlabel("Residence time (s)")
        ax.set_ylabel("Dry mole fraction")
        ax.set_title(f"{case}: Dry species vs time")
        ax.grid(True, alpha=0.3)
        ax.legend(ncol=2, fontsize="small")
        fig.tight_layout()
        fig.savefig(out_dir / "species_vs_time.png", dpi=180)
        plt.close(fig)

    def plot_velocity():
        col = "u_m_s"
        if col not in base_df.columns or col not in red_df.columns:
            return
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(base_df["z_m"], base_df[col], label="Baseline", lw=2)
        ax.plot(red_df["z_m"], red_df[col], label="Reduced", lw=2, ls="--")
        ax.set_xlabel("Axial position z (m)")
        ax.set_ylabel("Velocity (m/s)")
        ax.set_title(f"{case}: Axial velocity")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "u_vs_z.png", dpi=180)
        plt.close(fig)

    def plot_h2co_ratio():
        if "X_H2" not in base_df.columns or "X_CO" not in base_df.columns:
            return
        if "X_H2" not in red_df.columns or "X_CO" not in red_df.columns:
            return
        base_ratio = base_df["X_H2"] / base_df["X_CO"].clip(lower=1e-12)
        red_ratio = red_df["X_H2"] / red_df["X_CO"].clip(lower=1e-12)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(base_df["z_m"], base_ratio, label="Baseline", lw=2)
        ax.plot(red_df["z_m"], red_ratio, label="Reduced", lw=2, ls="--")
        ax.set_xlabel("Axial position z (m)")
        ax.set_ylabel("H₂/CO (mole ratio)")
        ax.set_title(f"{case}: H₂/CO ratio profile")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "H2CO_vs_z.png", dpi=180)
        plt.close(fig)

    plot_temperature()
    plot_species_vs_z()
    plot_species_vs_time()
    plot_velocity()
    plot_h2co_ratio()


def compare_cases(
    base_root: Path,
    reduced_root: Path,
    cases: Iterable[str],
    out_dir: Path,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    summaries: List[Dict[str, Any]] = []

    for case in cases:
        base_case = base_root / case
        red_case = reduced_root / case
        if not base_case.exists():
            raise FileNotFoundError(f"Baseline case directory missing: {base_case}")
        if not red_case.exists():
            raise FileNotFoundError(f"Reduced case directory missing: {red_case}")

        base_kpi = _load_json(base_case / "kpis.json")
        red_kpi = _load_json(red_case / "kpis.json")
        base_stats = _dry_species_from_kpi(base_kpi)
        red_stats = _dry_species_from_kpi(red_kpi)

        summary = {"case": case}
        for key in ("H2", "CO", "CO2", "CH4", "H2_CO", "T_out_K"):
            base_val = base_stats.get(key, float("nan"))
            red_val = red_stats.get(key, float("nan"))
            summary[f"{key}_base"] = base_val
            summary[f"{key}_red"] = red_val
            summary[f"{key}_err_%"] = _percent_diff(base_val, red_val)
        summaries.append(summary)

        base_axial = _load_axial(base_case / "axial_profiles.csv")
        red_axial = _load_axial(red_case / "axial_profiles.csv")
        _make_case_plot(case, base_axial, red_axial, out_dir / case)

    df = pd.DataFrame(summaries)
    csv_path = out_dir / "kpi_comparison.csv"
    df.to_csv(csv_path, index=False)
    print(f"[compare] KPI summary written to {csv_path.as_posix()}")
    print(f"[compare] Plots stored under {out_dir.as_posix()}/<case>/")
    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare baseline vs reduced HP-POX results (KPI + overlay plots)."
    )
    parser.add_argument(
        "--baseline",
        required=True,
        help="Baseline results root (folder containing per-case outputs).",
    )
    parser.add_argument(
        "--reduced",
        required=True,
        help="Reduced results root (folder containing per-case outputs).",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        required=True,
        help="Case names to compare (e.g. A_case1_richter A_case4_richter).",
    )
    parser.add_argument(
        "--out",
        default="comparisons/simple_compare",
        help="Output directory for summary and plots.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compare_cases(
        base_root=Path(args.baseline),
        reduced_root=Path(args.reduced),
        cases=args.cases,
        out_dir=Path(args.out),
    )


if __name__ == "__main__":
    main()
