#!/usr/bin/env python3
"""
Overlay plots comparing two mechanisms for the same case using axial_profiles.csv.

Usage:
  python tools/plot_overlay_unified.py --case A_case1_richter \
      --gri hp_pox_results_gri/A_case1_richter/axial_profiles.csv \
      --sd  hp_pox_results_sandiego/A_case1_richter/axial_profiles.csv \
      --out report/figs
"""

from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
from .common_io import ensure_noninteractive_backend, read_axial_csv, safe_mkdir


def safe_col(df: pd.DataFrame, col: str) -> bool:
    return col in df.columns and df[col].notna().any()


def plot_overlay(case: str, gri_csv: Path, sd_csv: Path, out_root: Path):
    out_dir = out_root / case
    out_dir.mkdir(parents=True, exist_ok=True)

    df_g = pd.read_csv(gri_csv)
    df_s = pd.read_csv(sd_csv)

    # T vs z
    z_col = 'z_m' if 'z_m' in df_g.columns else ('z' if 'z' in df_g.columns else None)
    if z_col and safe_col(df_g, z_col) and safe_col(df_g, 'T_C'):
        plt.figure(figsize=(10,6))
        plt.plot(df_g[z_col], df_g['T_C'], label='GRI T (°C)', color='tab:blue')
        plt.plot(df_s[z_col], df_s['T_C'], label='SanDiego T (°C)', color='tab:red')
        plt.xlabel('z (m)')
        plt.ylabel('Temperature (°C)')
        plt.title(f'Temperature vs z — {case}')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / 'T_vs_z_overlay.png', dpi=300)
        plt.close()

    # Species vs z (X_)
    species = ['H2','CO','CO2','CH4','H2O','O2','N2']
    plt.figure(figsize=(11,7))
    for sp in species:
        col = f'X_{sp}'
        if safe_col(df_g, col) and safe_col(df_s, col) and z_col:
            plt.plot(df_g[z_col], df_g[col], label=f'GRI {sp}', linestyle='-')
            plt.plot(df_s[z_col], df_s[col], label=f'SD {sp}', linestyle='--')
    if plt.gca().has_data():
        plt.xlabel('z (m)')
        plt.ylabel('Mole fraction (-)')
        plt.title(f'Species vs z — {case}')
        plt.grid(True, alpha=0.3)
        plt.legend(ncol=2, fontsize=9)
        plt.tight_layout()
        plt.savefig(out_dir / 'species_vs_z_overlay.png', dpi=300)
    plt.close()

    # H2/CO vs z
    have_ratio = False
    if 'H2_CO' in df_g.columns and 'H2_CO' in df_s.columns and z_col:
        plt.figure(figsize=(10,6))
        plt.semilogy(df_g[z_col], df_g['H2_CO'], label='GRI H2/CO', color='tab:blue')
        plt.semilogy(df_s[z_col], df_s['H2_CO'], label='SanDiego H2/CO', color='tab:red')
        plt.axhline(2.0, color='k', linestyle='--', linewidth=1.0, label='Target 2.0')
        plt.xlabel('z (m)')
        plt.ylabel('H2/CO (log)')
        plt.title(f'H2/CO vs z — {case}')
        plt.grid(True, which='both', alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / 'H2CO_vs_z_overlay.png', dpi=300)
        plt.close()
        have_ratio = True

    # Summary print
    print(f"Overlay plots saved: {out_dir.as_posix()} | ratio_plotted={have_ratio}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True)
    ap.add_argument('--gri', required=True)
    ap.add_argument('--sd', required=True)
    ap.add_argument('--out', default='report/figs')
    args = ap.parse_args()

    plot_overlay(args.case, Path(args.gri), Path(args.sd), Path(args.out))


def plot_T_overlay(base_csv: Union[str, Path], red_csv: Union[str, Path], 
                   out_png: Union[str, Path], title: str) -> None:
    """
    Baseline solid vs Reduced dashed, T vs z.
    
    Args:
        base_csv: Path to baseline CSV
        red_csv: Path to reduced CSV  
        out_png: Output PNG path
        title: Plot title
    """
    ensure_noninteractive_backend()
    
    try:
        df_base, z_col = read_axial_csv(base_csv)
        df_red, _ = read_axial_csv(red_csv)
        
        plt.figure(figsize=(10, 6))
        plt.plot(df_base[z_col], df_base['T_C'], 'k-', label='baseline', linewidth=2)
        plt.plot(df_red[z_col], df_red['T_C'], 'k--', label='reduced', linewidth=2)
        plt.xlabel('z (m)')
        plt.ylabel('T (°C)')
        plt.title(title)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        # Ensure output directory exists
        safe_mkdir(Path(out_png).parent)
        plt.savefig(out_png, dpi=200)
        plt.close()
        
    except Exception as e:
        print(f"Error plotting T overlay: {e}")
        plt.close()


def plot_species_overlay(base_csv: Union[str, Path], red_csv: Union[str, Path], 
                        out_png: Union[str, Path], species_cols: list, title: str) -> None:
    """
    Baseline solid vs Reduced dashed, log-y, legend.
    
    Args:
        base_csv: Path to baseline CSV
        red_csv: Path to reduced CSV
        out_png: Output PNG path
        species_cols: List of species column names (e.g., ['X_H2', 'X_CO'])
        title: Plot title
    """
    ensure_noninteractive_backend()
    
    try:
        df_base, z_col = read_axial_csv(base_csv)
        df_red, _ = read_axial_csv(red_csv)
        
        plt.figure(figsize=(11, 7))
        
        for col in species_cols:
            if col in df_base.columns and col in df_red.columns:
                species_name = col[2:] if col.startswith('X_') else col  # Remove X_ prefix
                plt.plot(df_base[z_col], df_base[col], label=f"{species_name} base", linewidth=2)
                plt.plot(df_red[z_col], df_red[col], '--', label=f"{species_name} red", linewidth=2)
        
        plt.xlabel('z (m)')
        plt.ylabel('Mole fraction')
        plt.yscale('log')
        plt.title(title)
        plt.grid(True, which='both', alpha=0.3)
        plt.legend(ncol=2, fontsize=9)
        plt.tight_layout()
        
        # Ensure output directory exists
        safe_mkdir(Path(out_png).parent)
        plt.savefig(out_png, dpi=200)
        plt.close()
        
    except Exception as e:
        print(f"Error plotting species overlay: {e}")
        plt.close()


if __name__ == '__main__':
    main()


