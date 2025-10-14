#!/usr/bin/env python3
"""
Plot unified axial profiles (physics + diagnostics) per case.

Generates 6 figures per CSV:
  1) T_vs_z.png
  2) species_vs_z.png
  3) H2CO_vs_z.png
  4) fallback_map.png
  5) stepsizes.png
  6) convergence.png

Output directory defaults to hp_pox_results/<case>/figs/ for each CSV.
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


def ensure_out_dir(csv_path: Path, base_out: Path | None) -> Path:
    case_dir = csv_path.parent
    out_dir = (base_out / case_dir.name / 'figs') if base_out else (case_dir / 'figs')
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def safe_col(df: pd.DataFrame, col: str) -> bool:
    return col in df.columns and df[col].notna().any()


def plot_temperature(df: pd.DataFrame, out_dir: Path):
    if not (safe_col(df, 'z_m') and (safe_col(df, 'T') or safe_col(df, 'T_C') or safe_col(df, 'temperature_K'))):
        return
    z = df['z_m'] if 'z_m' in df else df['z']
    if 'T' in df:
        T = df['T']
    elif 'temperature_K' in df:
        T = df['temperature_K']
    else:
        T = df['T_C'] + 273.15
    plt.figure(figsize=(10, 6))
    plt.plot(z, T, label='T (K)', color='tab:red')

    # Add target temperature markers
    if 'case1' in str(out_dir).lower():
        target_T = 1201 + 273.15  # Case 1 target
        plt.axhline(y=target_T, color='red', linestyle='--', alpha=0.7, label=f'Target: {target_T:.1f}K')
    elif 'case4' in str(out_dir).lower():
        target_T = 1401 + 273.15  # Case 4 target
        plt.axhline(y=target_T, color='red', linestyle='--', alpha=0.7, label=f'Target: {target_T:.1f}K')

    # Add kernel temperature marker at z=0
    kernel_T = 3000  # Default - should be read from diagnostics
    try:
        # Try to read kernel temperature from diagnostics
        import json
        diag_file = out_dir.parent / 'diagnostics.json'
        if diag_file.exists():
            with open(diag_file, 'r') as f:
                diag_data = json.load(f)
                kernel_T = diag_data.get('kernel_T', 3000)
    except Exception:
        pass

    plt.scatter([0], [kernel_T], marker='v', color='orange', s=100, label=f'Kernel T: {kernel_T:.0f}K', zorder=10)
    plt.annotate(f'Kernel: {kernel_T:.0f} K',
                xy=(0, kernel_T),
                xytext=(0.02, kernel_T),
                textcoords='data',
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="none", alpha=0.7),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"))

    plt.xlabel('z (m)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs z')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / 'T_vs_z.png', dpi=300)
    plt.close()


def plot_species(df: pd.DataFrame, out_dir: Path):
    z = df['z_m'] if 'z_m' in df else df.get('z')
    if z is None:
        return
    species_cols = ['X_H2', 'X_CO', 'X_CO2', 'X_CH4', 'X_H2O', 'X_O2', 'X_N2']
    present = [c for c in species_cols if safe_col(df, c)]
    if not present:
        return
    plt.figure(figsize=(10, 6))
    for c in present:
        plt.plot(z, 100.0 * df[c], label=c)

    # Add vertical line at z=0 to show the "jump"
    plt.axvline(z.iloc[0], color='k', linestyle='-', alpha=0.3, linewidth=1)

    # Add target composition markers
    if 'case1' in str(out_dir).lower():
        # Case 1 targets (dry basis)
        plt.axhline(y=48.27, color='red', linestyle='--', alpha=0.7, label='H2 Target: 48.27%')
        plt.axhline(y=23.79, color='blue', linestyle='--', alpha=0.7, label='CO Target: 23.79%')
        plt.axhline(y=4.19, color='green', linestyle='--', alpha=0.7, label='CO2 Target: 4.19%')
        plt.axhline(y=3.76, color='orange', linestyle='--', alpha=0.7, label='CH4 Target: 3.76%')
    elif 'case4' in str(out_dir).lower():
        # Case 4 targets (dry basis)
        plt.axhline(y=48.06, color='red', linestyle='--', alpha=0.7, label='H2 Target: 48.06%')
        plt.axhline(y=25.61, color='blue', linestyle='--', alpha=0.7, label='CO Target: 25.61%')
        plt.axhline(y=3.89, color='green', linestyle='--', alpha=0.7, label='CO2 Target: 3.89%')
        plt.axhline(y=0.06, color='orange', linestyle='--', alpha=0.7, label='CH4 Target: 0.06%')

    plt.xlabel('z (m)')
    plt.ylabel('Mole fraction (%)')
    plt.title('Species vs z')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'species_vs_z.png', dpi=300)
    plt.close()


def plot_h2co(df: pd.DataFrame, out_dir: Path):
    z = df['z_m'] if 'z_m' in df else df.get('z')
    if z is None:
        return
    # Prefer precomputed H2_CO if present
    if safe_col(df, 'H2_CO'):
        vals = df['H2_CO'].replace(0, np.nan)
    elif safe_col(df, 'X_H2') and safe_col(df, 'X_CO'):
        vals = df['X_H2'] / df['X_CO'].replace(0, np.nan)
    else:
        return
    if not np.isfinite(vals).any():
        return
    plt.figure(figsize=(10, 6))
    plt.plot(z, vals, color='tab:purple', label='H2/CO')
    plt.axhline(2.0, color='k', linestyle='--', label='Target 2.0')
    plt.yscale('log')
    plt.xlabel('z (m)')
    plt.ylabel('H2/CO (log)')
    plt.title('H2/CO vs z')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'H2CO_vs_z.png', dpi=300)
    plt.close()


def plot_species_vs_time(df: pd.DataFrame, out_dir: Path):
    if not safe_col(df, 'residence_time_s'):
        return
    t = df['residence_time_s']
    species_cols = ['X_H2', 'X_CO', 'X_CO2', 'X_CH4', 'X_H2O', 'X_O2', 'X_N2']
    present = [c for c in species_cols if safe_col(df, c)]
    if not present:
        return
    plt.figure(figsize=(10, 6))
    for c in present:
        plt.plot(t, 100.0 * df[c], label=c)

    # Add vertical line at t=0 to show the "jump"
    plt.axvline(t.iloc[0], color='k', linestyle='-', alpha=0.3, linewidth=1)

    plt.xlabel('Residence time (s)')
    plt.ylabel('Mole fraction (%)')
    plt.title('Species vs residence time')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'species_vs_time.png', dpi=300)
    plt.close()


def plot_velocity_pressure(df: pd.DataFrame, out_dir: Path):
    if not safe_col(df, 'z_m'):
        return
    z = df['z_m']
    plt.figure(figsize=(10, 6))
    # Velocity
    if safe_col(df, 'u_m_s'):
        plt.plot(z, df['u_m_s'], label='u (m/s)', color='tab:blue')

        # Add jet velocity annotation at z=0
        if len(z) > 0:
            jet_vel = df['u_m_s'].iloc[0]  # First value at inlet
            # Determine case-specific jet velocity annotation
            case_text = ""
            if 'case1' in str(out_dir).lower():
                case_text = "Case 1"
            elif 'case4' in str(out_dir).lower():
                case_text = "Case 4"

            if case_text:
                plt.scatter([0], [jet_vel], marker='o', color='red', s=50, label=f'Jet velocity: {jet_vel:.1f} m/s', zorder=10)
                plt.annotate(f'Jet: {jet_vel:.1f} m/s',
                            xy=(0, jet_vel),
                            xytext=(0.02, jet_vel),
                            textcoords='data',
                            fontsize=10,
                            bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="none", alpha=0.7),
                            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"))

    # Pressure (prefer p -> Pa)
    if safe_col(df, 'p'):
        pbar = df['p'] / 1e5
        plt.plot(z, pbar, label='P (bar)', color='tab:green')
    elif safe_col(df, 'pressure_Pa'):
        pbar = df['pressure_Pa'] / 1e5
        plt.plot(z, pbar, label='P (bar)', color='tab:green')
    else:
        # If not logged, assume constant 51 bar for context
        plt.plot(z, np.full_like(z, 51.0), '--', label='P (bar, assumed 51)', color='tab:green')
    plt.xlabel('z (m)')
    plt.ylabel('u (m/s), P (bar)')
    plt.title('Velocity and Pressure vs z')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'u_p_vs_z.png', dpi=300)
    plt.close()

def plot_fallback(df: pd.DataFrame, out_dir: Path):
    if not (safe_col(df, 'z_m') and safe_col(df, 'used_fallback')):
        return
    z = df['z_m']
    fb = df['used_fallback'].astype(float)
    mask = fb == 1.0
    plt.figure(figsize=(10, 3.5))
    if mask.any():
        plt.stem(z[mask], np.ones(mask.sum()), basefmt=' ', linefmt='r-', markerfmt='ro', label='Fallback')
    plt.ylim(-0.1, 1.1)
    plt.xlabel('z (m)')
    plt.ylabel('Fallback')
    plt.title('Fallback map')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'fallback_map.png', dpi=300)
    plt.close()


def plot_stepsizes(df: pd.DataFrame, out_dir: Path):
    if not (safe_col(df, 'z_m') and safe_col(df, 'dt_min') and safe_col(df, 'dt_effective')):
        return
    z = df['z_m']
    plt.figure(figsize=(10, 6))
    plt.fill_between(z, df['dt_min'], df['dt_effective'], alpha=0.25, color='tab:blue', label='dt range')
    plt.plot(z, df['dt_min'], 'b-', linewidth=1.0, label='dt_min')
    plt.plot(z, df['dt_effective'], 'b--', linewidth=1.5, label='dt_effective')
    plt.yscale('log')
    plt.xlabel('z (m)')
    plt.ylabel('Time step (s, log)')
    plt.title('Adaptive chemistry step sizes')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / 'stepsizes.png', dpi=300)
    plt.close()


def plot_convergence(df: pd.DataFrame, out_dir: Path):
    if not (safe_col(df, 'z_m') and (safe_col(df, 'substeps') or safe_col(df, 'halvings'))):
        return
    z = df['z_m']
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    if safe_col(df, 'substeps'):
        ax1.plot(z, df['substeps'], color='tab:orange', label='substeps')
        ax1.set_ylabel('Substeps')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    if safe_col(df, 'halvings'):
        ax2.plot(z, df['halvings'], color='tab:red', label='halvings')
        ax2.set_ylabel('Halvings')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    ax2.set_xlabel('z (m)')
    plt.tight_layout()
    plt.savefig(out_dir / 'convergence.png', dpi=300)
    plt.close()


def process_one(csv_path: Path, base_out: Path | None):
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
        return
    out_dir = ensure_out_dir(csv_path, base_out)
    # Allow alternate column names
    if 'z' in df.columns and 'z_m' not in df.columns:
        df = df.rename(columns={'z': 'z_m'})
    plot_temperature(df, out_dir)
    plot_species(df, out_dir)
    plot_h2co(df, out_dir)
    plot_fallback(df, out_dir)
    plot_stepsizes(df, out_dir)
    plot_convergence(df, out_dir)
    plot_species_vs_time(df, out_dir)
    plot_velocity_pressure(df, out_dir)
    # Extra: velocity-only plots
    try:
        if safe_col(df, 'z_m') and safe_col(df, 'u_m_s'):
            plt.figure(figsize=(10, 6))
            plt.plot(df['z_m'], df['u_m_s'], color='tab:blue', label='u (m/s)')
            plt.xlabel('z (m)')
            plt.ylabel('u (m/s)')
            plt.title('Velocity vs z')
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig(out_dir / 'u_vs_z.png', dpi=300)
            plt.close()
        if safe_col(df, 'residence_time_s') and safe_col(df, 'u_m_s'):
            plt.figure(figsize=(10, 6))
            plt.plot(df['residence_time_s'], df['u_m_s'], color='tab:blue', label='u (m/s)')
            plt.xlabel('Residence time (s)')
            plt.ylabel('u (m/s)')
            plt.title('Velocity vs residence time')
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig(out_dir / 'u_vs_time.png', dpi=300)
            plt.close()
    except Exception:
        pass
    # Print concise summary
    summary = {
        'cells': len(df),
        'fallbacks': int(df['used_fallback'].sum()) if 'used_fallback' in df else 0,
        'mean_substeps': float(df['substeps'].mean()) if 'substeps' in df else float('nan'),
        'mean_halvings': float(df['halvings'].mean()) if 'halvings' in df else float('nan'),
    }
    if 'H2_CO' in df:
        vals = pd.to_numeric(df['H2_CO'], errors='coerce')
        summary['H2_CO_min'] = float(np.nanmin(vals)) if np.isfinite(vals).any() else None
        summary['H2_CO_max'] = float(np.nanmax(vals)) if np.isfinite(vals).any() else None
    species_cols = [c for c in df.columns if c.startswith('X_')]
    print(f"Summary for {csv_path.parent.name}: {summary} | species: {species_cols}")


def main():
    ap = argparse.ArgumentParser(description='Plot unified axial profiles for one or more cases')
    ap.add_argument('paths', nargs='+', help='One or more CSVs (axial_profiles.csv) or case directories')
    ap.add_argument('--out', default=None, help='Base output directory (default: per case figs/ under each case)')
    args = ap.parse_args()
    base_out = Path(args.out) if args.out else None
    csvs: list[Path] = []
    for p in args.paths:
        path = Path(p)
        if path.is_dir():
            cand = path / 'axial_profiles.csv'
            if cand.exists():
                csvs.append(cand)
            else:
                print(f"Warning: {cand} not found; skipping {path}")
        elif path.is_file():
            csvs.append(path)
        else:
            print(f"Warning: {p} not found; skipping")
    if not csvs:
        print('No valid CSVs found.')
        sys.exit(1)
    for csv_path in csvs:
        process_one(csv_path, base_out)


if __name__ == '__main__':
    main()


