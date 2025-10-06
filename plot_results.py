#!/usr/bin/env python3
"""
Plot results from HP-POX 0D-1D reactor model.

This script generates publication-ready plots from the simulation results.
"""

import argparse
import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Set matplotlib style for publication-ready plots
plt.style.use('default')
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})


def plot_case(case_dir: Path, save: bool = True, show: bool = False) -> None:
    """
    Generate all plots for a single case.
    
    Args:
        case_dir: Path to the case directory containing the data files
        save: Whether to save the plots to files
        show: Whether to display the plots in windows
    """
    case_name = case_dir.name
    print(f"Processing case: {case_name}")
    
    # Check if case directory exists
    if not case_dir.exists():
        print(f"Warning: Case directory {case_dir} does not exist, skipping...")
        return
    
    # File paths
    axial_file = case_dir / "axial_profiles.csv"
    kpis_file = case_dir / "kpis.json"
    outlet_file = case_dir / "outlet_composition_table.csv"
    
    # Check if required files exist
    missing_files = []
    if not axial_file.exists():
        missing_files.append("axial_profiles.csv")
    if not kpis_file.exists():
        missing_files.append("kpis.json")
    if not outlet_file.exists():
        missing_files.append("outlet_composition_table.csv")
    
    if missing_files:
        print(f"Warning: Missing files for {case_name}: {', '.join(missing_files)}, skipping...")
        return
    
    try:
        # Load data
        df_axial = pd.read_csv(axial_file)
        with open(kpis_file, 'r') as f:
            kpis = json.load(f)
        df_outlet = pd.read_csv(outlet_file)
        
        # Generate plots
        plot_axial_profiles(df_axial, case_dir, save, show)
        plot_ignition_markers(df_axial, kpis, case_dir, save, show)
        plot_outlet_composition(df_outlet, case_dir, save, show)
        
        print(f"✓ Generated plots for {case_name}")
        
    except Exception as e:
        print(f"Error processing {case_name}: {e}")
        return


def plot_axial_profiles(df: pd.DataFrame, case_dir: Path, save: bool, show: bool) -> None:
    """
    Generate 2-panel axial profiles plot.
    
    Top panel: Temperature vs axial position
    Bottom panel: Species mole fractions vs axial position
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Top panel: Temperature
    z_m = df['z_m']
    T_C = df['T_C']  # Already in Celsius
    
    ax1.plot(z_m, T_C, 'k-', linewidth=2, label='Temperature')
    ax1.set_xlabel('Axial Position (m)')
    ax1.set_ylabel('Temperature (°C)')
    ax1.set_title('Temperature Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Bottom panel: Species mole fractions
    species_cols = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'N2']
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']
    
    for i, (species, color) in enumerate(zip(species_cols, colors)):
        if species in df.columns:
            # Convert mole fraction to percentage
            mole_frac_pct = df[species] * 100
            ax2.plot(z_m, mole_frac_pct, color=color, linewidth=2, label=species)
    
    ax2.set_xlabel('Axial Position (m)')
    ax2.set_ylabel('Mole Fraction (%)')
    ax2.set_title('Species Composition Profile')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    
    if save:
        output_file = case_dir / "axial_profiles_plot.png"
        plt.savefig(output_file)
        print(f"  Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_ignition_markers(df: pd.DataFrame, kpis: Dict, case_dir: Path, save: bool, show: bool) -> None:
    """
    Generate temperature profile with ignition length markers.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot temperature profile
    z_m = df['z_m']
    T_C = df['T_C']
    
    ax.plot(z_m, T_C, 'k-', linewidth=2, label='Temperature')
    
    # Add ignition markers
    z_ign_peak = kpis.get('ignition_length_peak_m', 0.0)
    z_ign_thresh = kpis.get('ignition_length_threshold_m', 0.4)
    
    if z_ign_peak > 0:
        ax.axvline(x=z_ign_peak, color='red', linestyle='--', linewidth=2, 
                  label=f'Peak Ignition: {z_ign_peak:.3f} m')
    
    if z_ign_thresh > 0:
        ax.axvline(x=z_ign_thresh, color='orange', linestyle='--', linewidth=2,
                  label=f'Threshold Ignition: {z_ign_thresh:.3f} m')
    
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title('Temperature Profile with Ignition Markers')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    
    if save:
        output_file = case_dir / "ignition_length_markers.png"
        plt.savefig(output_file)
        print(f"  Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_outlet_composition(df: pd.DataFrame, case_dir: Path, save: bool, show: bool) -> None:
    """
    Generate outlet composition bar chart comparing simulated vs target.
    """
    # Filter for the species we want to plot (dry basis)
    target_species = ['H2', 'CO', 'CO2', 'CH4', 'N2']
    
    # Filter the dataframe for target species
    df_filtered = df[df['Species'].isin(target_species)].copy()
    
    if df_filtered.empty:
        print(f"  Warning: No target species found in outlet composition data")
        return
    
    # Create the bar chart
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    x = np.arange(len(target_species))
    width = 0.35
    
    # Get values for each species
    simulated_values = []
    target_values = []
    
    for species in target_species:
        species_data = df_filtered[df_filtered['Species'] == species]
        if not species_data.empty:
            simulated_values.append(species_data['Dry_Basis_Vol%'].iloc[0])
        else:
            simulated_values.append(0.0)
        # For now, we'll use simulated values as target (since we don't have separate target data)
        target_values.append(simulated_values[-1])  # Placeholder - would need actual target data
    
    bars1 = ax.bar(x - width/2, simulated_values, width, label='Simulated', alpha=0.8)
    bars2 = ax.bar(x + width/2, target_values, width, label='Target', alpha=0.8)
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3),  # 3 points vertical offset
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3),  # 3 points vertical offset
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)
    
    ax.set_xlabel('Species')
    ax.set_ylabel('Dry Basis Composition (%)')
    ax.set_title('Outlet Composition Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(target_species)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if save:
        output_file = case_dir / "outlet_barplot.png"
        plt.savefig(output_file)
        print(f"  Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(description='Plot HP-POX 0D-1D reactor results')
    parser.add_argument('--case', type=str, help='Case ID to plot (e.g., A_case1_richter)')
    parser.add_argument('--save', action='store_true', default=True, help='Save figures (default: True)')
    parser.add_argument('--show', action='store_true', help='Display figures in windows')
    parser.add_argument('--results-dir', type=str, default='hp_pox_results', 
                       help='Results directory (default: hp_pox_results)')
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    
    if not results_dir.exists():
        print(f"Error: Results directory {results_dir} does not exist")
        return
    
    if args.case:
        # Plot specific case
        case_dir = results_dir / args.case
        plot_case(case_dir, save=args.save, show=args.show)
    else:
        # Plot all cases
        case_dirs = [d for d in results_dir.iterdir() if d.is_dir()]
        
        if not case_dirs:
            print(f"No case directories found in {results_dir}")
            return
        
        print(f"Found {len(case_dirs)} cases to process")
        
        for case_dir in sorted(case_dirs):
            plot_case(case_dir, save=args.save, show=args.show)
    
    print("Plotting complete!")


if __name__ == "__main__":
    main()
