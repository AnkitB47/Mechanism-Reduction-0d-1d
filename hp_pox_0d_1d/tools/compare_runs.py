#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import cantera as ct  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[1]
SPECIES_PLOT_ORDER = ["H2", "CO", "CO2", "H2O", "CH4", "O2", "H2CO"]
MAJOR_SPECIES_FOR_RMSE = ["H2", "CO", "CO2", "CH4", "H2O", "O2", "H2CO"]


@dataclass
class RunData:
    case: str
    full_dir: Path
    reduced_dir: Path
    out_dir: Path
    figures_dir: Path
    full_df: pd.DataFrame
    reduced_df: pd.DataFrame
    full_kpis: Dict[str, float]
    reduced_kpis: Dict[str, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare full vs reduced HP-POX runs")
    parser.add_argument("--full", required=True, help="Directory with full-run results")
    parser.add_argument("--reduced", required=True, help="Directory with reduced-run results")
    parser.add_argument("--cases", nargs="+", required=True, help="Case identifiers to compare")
    parser.add_argument("--outdir", required=True, help="Comparison output directory")
    parser.add_argument("--full-mech", default="mechanisms/full_hp_pox.yaml", help="Path to full mechanism YAML")
    parser.add_argument("--reduced-mech", default=None, help="Path to reduced mechanism YAML")
    return parser.parse_args()


def resolve_path(path_str: str) -> Path:
    path = Path(path_str)
    return path if path.is_absolute() else (REPO_ROOT / path)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_json(path: Path) -> Dict[str, float]:
    if not path.exists():
        return {}
    with path.open("r") as f:
        return json.load(f)


def detect_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    for name in candidates:
        if name in df.columns:
            return name
    return None


def align_series(
    coord_full: np.ndarray,
    values_full: np.ndarray,
    coord_red: np.ndarray,
    values_red: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    if len(coord_full) == 0 or len(coord_red) == 0:
        return np.array([]), np.array([])
    if np.allclose(coord_full[: len(coord_red)], coord_red[: len(coord_full)], rtol=1e-6, atol=1e-9):
        n = min(len(values_full), len(values_red))
        return values_full[:n], values_red[:n]
    if np.all(np.diff(coord_red) >= 0):
        interp_red = np.interp(coord_full, coord_red, values_red, left=np.nan, right=np.nan)
        mask = ~np.isnan(interp_red)
        return values_full[mask], interp_red[mask]
    n = min(len(values_full), len(values_red))
    return values_full[:n], values_red[:n]


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    diff = a - b
    return float(np.sqrt(np.mean(diff * diff)))


def load_run_data(case: str, args: argparse.Namespace, out_root: Path) -> RunData:
    full_case_dir = resolve_path(args.full) / case
    reduced_case_dir = resolve_path(args.reduced) / case
    if not full_case_dir.is_dir():
        raise FileNotFoundError(f"Full run results missing for case {case}: {full_case_dir}")
    if not reduced_case_dir.is_dir():
        raise FileNotFoundError(f"Reduced run results missing for case {case}: {reduced_case_dir}")

    axial_name = "axial_profiles.csv"
    full_csv = full_case_dir / axial_name
    reduced_csv = reduced_case_dir / axial_name
    if not full_csv.exists() or not reduced_csv.exists():
        raise FileNotFoundError(f"Expected axial_profiles.csv in both runs for case {case}")

    full_df = pd.read_csv(full_csv)
    reduced_df = pd.read_csv(reduced_csv)
    full_kpis = load_json(full_case_dir / "kpis.json")
    reduced_kpis = load_json(reduced_case_dir / "kpis.json")

    case_out = out_root / case
    figs_dir = case_out / "figs"
    ensure_dir(figs_dir)
    ensure_dir(case_out / "full_src")
    ensure_dir(case_out / "reduced_src")

    return RunData(
        case=case,
        full_dir=full_case_dir,
        reduced_dir=reduced_case_dir,
        out_dir=case_out,
        figures_dir=figs_dir,
        full_df=full_df,
        reduced_df=reduced_df,
        full_kpis=full_kpis,
        reduced_kpis=reduced_kpis,
    )


def copy_original_figs(run_data: RunData) -> None:
    for label, src_dir in (("full_src", run_data.full_dir / "figs"), ("reduced_src", run_data.reduced_dir / "figs")):
        if not src_dir.is_dir():
            print(f"[compare] {run_data.case}: skipping copy of {'full' if label=='full_src' else 'reduced'} figs (missing {src_dir})")
            continue
        dst_dir = run_data.out_dir / label
        for src_file in src_dir.glob("*.png"):
            shutil.copy2(src_file, dst_dir / src_file.name)


def plot_overlay(
    run_data: RunData,
    x_full: np.ndarray,
    y_full: np.ndarray,
    x_red: np.ndarray,
    y_red: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    filename: str,
) -> None:
    if len(x_full) == 0 or len(x_red) == 0:
        print(f"[compare] {run_data.case}: skipping {filename} (insufficient data)")
        return

    plt.figure(figsize=(7, 4))
    plt.plot(x_full, y_full, label="Full", linewidth=2)
    plt.plot(x_red, y_red, label="Reduced", linestyle="--", linewidth=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    diff = np.abs(np.interp(x_full, x_red, y_red, left=np.nan, right=np.nan) - y_full)
    max_diff = np.nanmax(diff) if diff.size else float("nan")
    plt.text(0.02, 0.95, f"max |Δ| = {max_diff:.3g}", transform=plt.gca().transAxes, fontsize=9, va="top")
    plt.legend()
    plt.tight_layout()
    plt.savefig(run_data.figures_dir / filename, dpi=200)
    plt.close()


def collect_species_column(df: pd.DataFrame, species: str) -> Optional[str]:
    return detect_column(df, [f"X_{species}", species])


def build_plots(run_data: RunData) -> None:
    df_full = run_data.full_df
    df_red = run_data.reduced_df

    z_col_full = detect_column(df_full, ["z_m", "z", "distance_m"])
    z_col_red = detect_column(df_red, ["z_m", "z", "distance_m"])
    t_col_full = detect_column(df_full, ["t_s", "time_s", "time"])
    t_col_red = detect_column(df_red, ["t_s", "time_s", "time"])

    z_full = df_full[z_col_full].to_numpy() if z_col_full else np.array([])
    z_red = df_red[z_col_red].to_numpy() if z_col_red else np.array([])
    t_full = df_full[t_col_full].to_numpy() if t_col_full else np.array([])
    t_red = df_red[t_col_red].to_numpy() if t_col_red else np.array([])

    # Species vs z
    available_species = []
    for sp in SPECIES_PLOT_ORDER:
        col_full = collect_species_column(df_full, sp)
        col_red = collect_species_column(df_red, sp)
        if col_full and col_red:
            available_species.append((sp, col_full, col_red))
    if available_species and z_full.size and z_red.size:
        plt.figure(figsize=(8, 5))
        for sp, col_f, col_r in available_species:
            aligned_full, aligned_red = align_series(z_full, df_full[col_f].to_numpy(), z_red, df_red[col_r].to_numpy())
            if aligned_full.size == 0:
                continue
            plt.plot(z_full[: len(aligned_full)], aligned_full, label=f"{sp} (Full)")
            plt.plot(z_full[: len(aligned_red)], aligned_red, linestyle="--", label=f"{sp} (Reduced)")
        plt.xlabel("z [m]")
        plt.ylabel("Mole fraction")
        plt.title(f"{run_data.case} — Species profiles")
        plt.legend(ncol=2, fontsize=9)
        plt.tight_layout()
        plt.savefig(run_data.figures_dir / "species_vs_z.png", dpi=200)
        plt.close()

    # Species vs time
    if available_species and t_full.size and t_red.size:
        plt.figure(figsize=(8, 5))
        for sp, col_f, col_r in available_species:
            aligned_full, aligned_red = align_series(t_full, df_full[col_f].to_numpy(), t_red, df_red[col_r].to_numpy())
            if aligned_full.size == 0:
                continue
            plt.plot(t_full[: len(aligned_full)], aligned_full, label=f"{sp} (Full)")
            plt.plot(t_full[: len(aligned_red)], aligned_red, linestyle="--", label=f"{sp} (Reduced)")
        plt.xlabel("time [s]")
        plt.ylabel("Mole fraction")
        plt.title(f"{run_data.case} — Species vs time")
        plt.legend(ncol=2, fontsize=9)
        plt.tight_layout()
        plt.savefig(run_data.figures_dir / "species_vs_time.png", dpi=200)
        plt.close()

    # Temperature vs z
    temp_col_full = detect_column(df_full, ["T_K", "T_C"])
    temp_col_red = detect_column(df_red, ["T_K", "T_C"])
    if z_full.size and z_red.size and temp_col_full and temp_col_red:
        temp_full = df_full[temp_col_full].to_numpy()
        temp_red = df_red[temp_col_red].to_numpy()
        if temp_col_full.endswith("_C"):
            temp_full = temp_full + 273.15
        if temp_col_red.endswith("_C"):
            temp_red = temp_red + 273.15
        plot_overlay(
            run_data,
            z_full,
            temp_full,
            z_red,
            temp_red,
            xlabel="z [m]",
            ylabel="Temperature [K]",
            title=f"{run_data.case} — Temperature vs z",
            filename="T_vs_z.png",
        )

    # Velocity vs z
    vel_full_col = detect_column(df_full, ["u_m_per_s", "velocity_m_per_s", "u"])
    vel_red_col = detect_column(df_red, ["u_m_per_s", "velocity_m_per_s", "u"])
    if z_full.size and z_red.size and vel_full_col and vel_red_col:
        plot_overlay(
            run_data,
            z_full,
            df_full[vel_full_col].to_numpy(),
            z_red,
            df_red[vel_red_col].to_numpy(),
            xlabel="z [m]",
            ylabel="Velocity [m/s]",
            title=f"{run_data.case} — Velocity vs z",
            filename="u_vs_z.png",
        )

    # Velocity vs time
    if t_full.size and t_red.size and vel_full_col and vel_red_col:
        plot_overlay(
            run_data,
            t_full,
            df_full[vel_full_col].to_numpy(),
            t_red,
            df_red[vel_red_col].to_numpy(),
            xlabel="time [s]",
            ylabel="Velocity [m/s]",
            title=f"{run_data.case} — Velocity vs time",
            filename="u_vs_time.png",
        )

    # Velocity & pressure vs z
    pressure_full_col = detect_column(df_full, ["p_Pa", "pressure_Pa", "p"])
    pressure_red_col = detect_column(df_red, ["p_Pa", "pressure_Pa", "p"])
    if z_full.size and z_red.size and vel_full_col and vel_red_col and pressure_full_col and pressure_red_col:
        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax2 = ax1.twinx()
        aligned_vel_full, aligned_vel_red = align_series(
            z_full, df_full[vel_full_col].to_numpy(), z_red, df_red[vel_red_col].to_numpy()
        )
        aligned_p_full, aligned_p_red = align_series(
            z_full, df_full[pressure_full_col].to_numpy(), z_red, df_red[pressure_red_col].to_numpy()
        )
        if aligned_vel_full.size:
            ax1.plot(z_full[: len(aligned_vel_full)], aligned_vel_full, label="Velocity (Full)", color="tab:blue")
            ax1.plot(
                z_full[: len(aligned_vel_red)],
                aligned_vel_red,
                linestyle="--",
                label="Velocity (Reduced)",
                color="tab:blue",
            )
        if aligned_p_full.size:
            ax2.plot(
                z_full[: len(aligned_p_full)], aligned_p_full / 1e5, label="Pressure (Full)", color="tab:red"
            )
            ax2.plot(
                z_full[: len(aligned_p_red)],
                aligned_p_red / 1e5,
                linestyle="--",
                label="Pressure (Reduced)",
                color="tab:red",
            )
        ax1.set_xlabel("z [m]")
        ax1.set_ylabel("Velocity [m/s]", color="tab:blue")
        ax2.set_ylabel("Pressure [bar]", color="tab:red")
        ax1.set_title(f"{run_data.case} — Velocity & Pressure vs z")
        handles, labels = [], []
        for ax in (ax1, ax2):
            h, l = ax.get_legend_handles_labels()
            handles.extend(h)
            labels.extend(l)
        if handles:
            fig.legend(handles, labels, loc="upper right")
        fig.tight_layout()
        fig.savefig(run_data.figures_dir / "u_p_vs_z.png", dpi=200)
        plt.close(fig)

    # H2CO vs z
    h2co_full_col = collect_species_column(df_full, "H2CO")
    h2co_red_col = collect_species_column(df_red, "H2CO")
    if z_full.size and z_red.size and h2co_full_col and h2co_red_col:
        plot_overlay(
            run_data,
            z_full,
            df_full[h2co_full_col].to_numpy(),
            z_red,
            df_red[h2co_red_col].to_numpy(),
            xlabel="z [m]",
            ylabel="Mole fraction",
            title=f"{run_data.case} — H2CO vs z",
            filename="H2CO_vs_z.png",
        )


def compute_species_counts(full_mech: Optional[str], reduced_mech: Optional[str]) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
    def _load_counts(path: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
        if path is None:
            return None, None
        mech_path = resolve_path(path)
        if not mech_path.exists():
            return None, None
        try:
            sol = ct.Solution(str(mech_path))
            return sol.n_species, sol.n_reactions
        except Exception:
            return None, None

    return (*_load_counts(full_mech), *_load_counts(reduced_mech))


def compute_case_metrics(run_data: RunData) -> Dict[str, float]:
    df_full = run_data.full_df
    df_red = run_data.reduced_df
    z_full_col = detect_column(df_full, ["z_m", "z", "distance_m"])
    z_red_col = detect_column(df_red, ["z_m", "z", "distance_m"])
    temp_col_full = detect_column(df_full, ["T_K", "T_C"])
    temp_col_red = detect_column(df_red, ["T_K", "T_C"])
    metrics: Dict[str, float] = {}

    if z_full_col and z_red_col and temp_col_full and temp_col_red:
        z_full = df_full[z_full_col].to_numpy()
        z_red = df_red[z_red_col].to_numpy()
        temp_full = df_full[temp_col_full].to_numpy()
        temp_red = df_red[temp_col_red].to_numpy()
        if temp_col_full.endswith("_C"):
            temp_full += 273.15
        if temp_col_red.endswith("_C"):
            temp_red += 273.15
        aligned_full, aligned_red = align_series(z_full, temp_full, z_red, temp_red)
        metrics["max_abs_dT_K"] = float(np.max(np.abs(aligned_full - aligned_red))) if aligned_full.size else float("nan")
    else:
        metrics["max_abs_dT_K"] = float("nan")

    coord_full = df_full[z_full_col].to_numpy() if z_full_col else np.arange(len(df_full))
    coord_red = df_red[z_red_col].to_numpy() if z_red_col else np.arange(len(df_red))

    for sp in MAJOR_SPECIES_FOR_RMSE:
        col_full = collect_species_column(df_full, sp)
        col_red = collect_species_column(df_red, sp)
        key = f"RMSE_X_{sp}"
        if not col_full or not col_red:
            metrics[key] = float("nan")
            continue
        aligned_full, aligned_red = align_series(coord_full, df_full[col_full].to_numpy(), coord_red, df_red[col_red].to_numpy())
        metrics[key] = rmse(aligned_full, aligned_red)

    full_slip = run_data.full_kpis.get("CH4_slip_pct") or run_data.full_kpis.get("CH4_slip_percent")
    reduced_slip = run_data.reduced_kpis.get("CH4_slip_pct") or run_data.reduced_kpis.get("CH4_slip_percent")
    if full_slip is not None and reduced_slip is not None:
        metrics["CH4_slip_percent"] = abs(float(reduced_slip) - float(full_slip))
    elif reduced_slip is not None:
        metrics["CH4_slip_percent"] = float(reduced_slip)
    else:
        metrics["CH4_slip_percent"] = float("nan")

    h2co_full = run_data.full_kpis.get("H2_CO_ratio") or run_data.full_kpis.get("H2_CO_out") or run_data.full_kpis.get("H2_CO_out_avg")
    h2co_red = run_data.reduced_kpis.get("H2_CO_ratio") or run_data.reduced_kpis.get("H2_CO_out") or run_data.reduced_kpis.get("H2_CO_out_avg")
    if h2co_full and h2co_red and float(h2co_full) != 0.0:
        metrics["H2_CO_rel_err"] = abs(float(h2co_red) - float(h2co_full)) / abs(float(h2co_full))
    else:
        metrics["H2_CO_rel_err"] = float("nan")

    return metrics


def determine_pass(metrics: Dict[str, float]) -> bool:
    temp_ok = metrics.get("max_abs_dT_K", float("inf")) <= 50.0
    rmse_ok = all(metrics.get(f"RMSE_X_{sp}", float("inf")) <= 0.10 for sp in ["H2", "CO", "CO2", "CH4", "H2O", "O2"])
    h2co_ok = metrics.get("H2_CO_rel_err", float("inf")) <= 0.15
    return temp_ok and rmse_ok and h2co_ok


def main() -> None:
    args = parse_args()
    out_root = resolve_path(args.outdir)
    ensure_dir(out_root)

    full_mech_path = args.full_mech
    reduced_mech_path = args.reduced_mech

    summary_rows: List[Dict[str, object]] = []

    for case in args.cases:
        run_data = load_run_data(case, args, out_root)
        copy_original_figs(run_data)
        build_plots(run_data)
        metrics = compute_case_metrics(run_data)
        metrics["PASS"] = determine_pass(metrics)

        metrics_path = run_data.out_dir / "metrics.json"
        with metrics_path.open("w") as f:
            json.dump(metrics, f, indent=2)

        n_species_full, n_reactions_full, n_species_red, n_reactions_red = compute_species_counts(
            full_mech_path, reduced_mech_path
        )

        summary_rows.append(
            {
                "case": case,
                "N_species_full": n_species_full,
                "N_reactions_full": n_reactions_full,
                "N_species_reduced": n_species_red,
                "N_reactions_reduced": n_reactions_red,
                "max_abs_dT_K": metrics.get("max_abs_dT_K"),
                "RMSE_X_H2": metrics.get("RMSE_X_H2"),
                "RMSE_X_CO": metrics.get("RMSE_X_CO"),
                "RMSE_X_CO2": metrics.get("RMSE_X_CO2"),
                "RMSE_X_CH4": metrics.get("RMSE_X_CH4"),
                "RMSE_X_H2O": metrics.get("RMSE_X_H2O"),
                "RMSE_X_O2": metrics.get("RMSE_X_O2"),
                "RMSE_X_H2CO": metrics.get("RMSE_X_H2CO"),
                "CH4_slip_percent": metrics.get("CH4_slip_percent"),
                "H2_CO_rel_err": metrics.get("H2_CO_rel_err"),
                "PASS": metrics.get("PASS"),
            }
        )

        status = "PASS" if metrics["PASS"] else "FAIL"
        print(f"[compare] {case} → {status}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = out_root / "summary_metrics.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[compare] Summary written to {summary_path}")


if __name__ == "__main__":
    main()
