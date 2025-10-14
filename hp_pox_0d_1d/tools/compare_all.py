#!/usr/bin/env python3
"""
Compare baseline vs reduced mechanism results across all cases.
Generates overlay plots and KPI comparison CSV.

Usage:
  python tools/compare_all.py --base-root hp_pox_results_gri \
      --red-root reduction/validate/gri --out comparison
"""

import argparse
import shutil
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from tools.common_io import ensure_noninteractive_backend, safe_mkdir
from tools.plot_overlay_unified import plot_T_overlay, plot_species_overlay
from tools.generate_kpi_compare import build_kpi_row


def main():
    parser = argparse.ArgumentParser(description='Compare baseline vs reduced results')
    parser.add_argument('--base-root', default='hp_pox_results_gri',
                        help='Root folder containing baseline results')
    parser.add_argument('--red-root', default='reduction/validate/gri',
                        help='Root folder containing reduced results')
    parser.add_argument('--out', default='comparison',
                        help='Output directory for comparison files')
    parser.add_argument('--cases', nargs='+', 
                        default=['A_case1_richter', 'A_case4_richter', 'B_case1_134L', 'B_case4_134L'],
                        help='List of case names to process')
    
    args = parser.parse_args()
    
    # Set up paths
    base_root = Path(args.base_root)
    red_root = Path(args.red_root)
    out_root = Path(args.out)
    
    # Ensure non-interactive backend
    ensure_noninteractive_backend()
    
    # Create output directory
    safe_mkdir(out_root)
    
    # Species columns to plot
    species_cols = ["X_H2", "X_CO", "X_CO2", "X_CH4", "X_H2O", "X_O2", "X_N2"]
    
    rows = []
    
    print(f"Comparing {len(args.cases)} cases...")
    print(f"Base root: {base_root}")
    print(f"Red root: {red_root}")
    print(f"Output: {out_root}")
    
    for case in args.cases:
        print(f"Processing {case}...")
        
        # File paths
        base_csv = base_root / case / 'axial_profiles.csv'
        red_csv = red_root / case / 'axial_profiles.csv'
        base_kpi = base_root / case / 'kpis.json'
        red_kpi = red_root / case / 'kpis.json'
        
        # Check if files exist
        if not base_csv.exists():
            print(f"  Warning: {base_csv} not found, skipping")
            continue
        if not red_csv.exists():
            print(f"  Warning: {red_csv} not found, skipping")
            continue
            
        # Create case figures directory
        figs_dir = red_root / case / 'figs'
        safe_mkdir(figs_dir)
        
        # Generate overlay plots
        try:
            plot_T_overlay(base_csv, red_csv, figs_dir / "T_vs_z.png", f"T vs z — {case}")
            plot_species_overlay(base_csv, red_csv, figs_dir / "species_vs_z.png", 
                               species_cols, f"Species vs z — {case}")
            
            # Copy plots to comparison directory
            shutil.copy(figs_dir / "T_vs_z.png", out_root / f"{case}_T_vs_z.png")
            shutil.copy(figs_dir / "species_vs_z.png", out_root / f"{case}_species_vs_z.png")
            
            print(f"  ✓ Plots generated and copied")
            
        except Exception as e:
            print(f"  Error generating plots for {case}: {e}")
        
        # Collect KPI data
        if base_kpi.exists() and red_kpi.exists():
            try:
                row = build_kpi_row(case, base_kpi, red_kpi)
                if row:
                    rows.append(row)
                    print(f"  ✓ KPIs collected")
                else:
                    print(f"  Warning: No KPI data for {case}")
            except Exception as e:
                print(f"  Error collecting KPIs for {case}: {e}")
        else:
            print(f"  Warning: KPI files missing for {case}")
    
    # Save KPI comparison CSV
    if rows:
        df_out = pd.DataFrame(rows)
        out_csv = out_root / 'kpi_compare_gri.csv'
        df_out.to_csv(out_csv, index=False)
        print(f"\nKPI comparison saved: {out_csv}")
        print(df_out.to_string(index=False))
    else:
        print("\nWarning: No KPI data collected")
    
    print(f"\nComparison complete!")
    print(f"Per-case figs: {red_root}/<case>/figs/")
    print(f"Comparison figs + CSV: {out_root}/")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
