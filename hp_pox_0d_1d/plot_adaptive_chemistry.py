#!/usr/bin/env python3
"""
Comprehensive plotting script for adaptive chemistry PFR results.
Generates T(x), H₂/CO, fallback map, and adaptive step size visualizations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import json
import argparse
from typing import Dict, Any, Optional

def load_results(case_dir: str) -> Dict[str, Any]:
    """Load PFR results from case directory."""
    case_path = Path(case_dir)
    
    # Load axial profiles
    axial_file = case_path / "axial_profiles.csv"
    if not axial_file.exists():
        raise FileNotFoundError(f"Axial profiles not found: {axial_file}")
    
    axial_data = pd.read_csv(axial_file)
    
    # Load KPIs
    kpis_file = case_path / "kpis.json"
    kpis = {}
    if kpis_file.exists():
        with open(kpis_file, 'r') as f:
            kpis = json.load(f)
    
    # Load chemistry diagnostics
    chem_file = case_path / "chemistry_diagnostics.csv"
    chem_data = None
    if chem_file.exists():
        chem_data = pd.read_csv(chem_file)
    
    return {
        'axial': axial_data,
        'kpis': kpis,
        'chemistry': chem_data,
        'case_path': case_path
    }

def plot_temperature_profile(ax, data: pd.DataFrame, kpis: Dict[str, Any]):
    """Plot temperature profile with ignition markers."""
    z = data['z_m']
    T = data['T_C']
    
    ax.plot(z, T, 'b-', linewidth=2, label='Temperature')
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Temperature (°C)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Temperature Profile')
    
    # Add ignition markers if available
    if 'z_ign_peak_m' in kpis and kpis['z_ign_peak_m'] > 0:
        ax.axvline(kpis['z_ign_peak_m'], color='red', linestyle='--', 
                  label=f'Peak dT/dz: {kpis["z_ign_peak_m"]:.3f} m')
    
    if 'z_ign_thresh_m' in kpis and kpis['z_ign_thresh_m'] > 0:
        ax.axvline(kpis['z_ign_thresh_m'], color='orange', linestyle='--', 
                  label=f'Threshold: {kpis["z_ign_thresh_m"]:.3f} m')
    
    ax.legend()

def plot_species_profiles(ax, data: pd.DataFrame):
    """Plot major species profiles."""
    z = data['z_m']
    
    # Major species to plot
    species_to_plot = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'O2', 'N2']
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'gray']
    
    for species, color in zip(species_to_plot, colors):
        if species in data.columns:
            ax.plot(z, data[species] * 100, color=color, linewidth=2, label=species)
    
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Mole Fraction (%)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Species Profiles')
    ax.legend()

def plot_h2_co_ratio(ax, data: pd.DataFrame, kpis: Dict[str, Any]):
    """Plot H₂/CO ratio profile."""
    z = data['z_m']
    
    # Calculate H₂/CO ratio
    if 'H2' in data.columns and 'CO' in data.columns:
        h2_co = data['H2'] / np.maximum(data['CO'], 1e-10)
        ax.plot(z, h2_co, 'r-', linewidth=2, label='H₂/CO Ratio')
        
        # Add target line
        ax.axhline(2.0, color='black', linestyle='--', alpha=0.7, label='Target: 2.0')
        
        # Add outlet value
        if len(h2_co) > 0:
            outlet_ratio = h2_co.iloc[-1]
            ax.text(0.02, 0.98, f'Outlet H₂/CO: {outlet_ratio:.2f}', 
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('H₂/CO Ratio')
    ax.grid(True, alpha=0.3)
    ax.set_title('H₂/CO Ratio Profile')
    ax.legend()
    ax.set_yscale('log')

def plot_chemistry_fallback_map(ax, chem_data: pd.DataFrame):
    """Plot chemistry fallback usage map."""
    if chem_data is None:
        ax.text(0.5, 0.5, 'No Chemistry Diagnostics Available', 
               transform=ax.transAxes, ha='center', va='center')
        ax.set_title('Chemistry Fallback Map')
        return
    
    z = chem_data['z']
    fallback_used = chem_data['fallback_used']
    
    # Create fallback map
    ax.scatter(z[fallback_used], np.ones(fallback_used.sum()), 
              c='red', s=20, alpha=0.7, label='Fallback Used')
    ax.scatter(z[~fallback_used], np.ones((~fallback_used).sum()), 
              c='green', s=20, alpha=0.7, label='Direct Chemistry')
    
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Chemistry Method')
    ax.set_ylim(0.5, 1.5)
    ax.set_yticks([1])
    ax.set_yticklabels(['Chemistry'])
    ax.grid(True, alpha=0.3)
    ax.set_title('Chemistry Fallback Map')
    ax.legend()

def plot_adaptive_step_sizes(ax, chem_data: pd.DataFrame):
    """Plot adaptive step size evolution."""
    if chem_data is None:
        ax.text(0.5, 0.5, 'No Chemistry Diagnostics Available', 
               transform=ax.transAxes, ha='center', va='center')
        ax.set_title('Adaptive Step Sizes')
        return
    
    z = chem_data['z']
    dt_min = chem_data['dt_min']
    dt_max = chem_data['dt_max']
    substeps = chem_data['substeps']
    
    # Plot step size range
    ax.fill_between(z, dt_min, dt_max, alpha=0.3, color='blue', label='Step Size Range')
    ax.plot(z, dt_min, 'b-', linewidth=1, label='Min Step')
    ax.plot(z, dt_max, 'b--', linewidth=1, label='Max Step')
    
    # Plot substeps (scaled)
    ax2 = ax.twinx()
    ax2.plot(z, substeps, 'r-', linewidth=2, label='Substeps')
    ax2.set_ylabel('Number of Substeps', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Time Step (s)')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Adaptive Step Sizes')
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')

def plot_chemistry_convergence(ax, chem_data: pd.DataFrame):
    """Plot chemistry convergence statistics."""
    if chem_data is None:
        ax.text(0.5, 0.5, 'No Chemistry Diagnostics Available', 
               transform=ax.transAxes, ha='center', va='center')
        ax.set_title('Chemistry Convergence')
        return
    
    z = chem_data['z']
    halvings = chem_data['halvings']
    converged = chem_data['converged']
    residual = chem_data['residual']
    
    # Plot halvings
    ax.plot(z, halvings, 'r-', linewidth=2, label='Halvings')
    
    # Plot convergence status
    ax2 = ax.twinx()
    ax2.scatter(z[converged], np.ones(converged.sum()), 
               c='green', s=20, alpha=0.7, label='Converged')
    ax2.scatter(z[~converged], np.ones((~converged).sum()), 
               c='red', s=20, alpha=0.7, label='Not Converged')
    ax2.set_ylabel('Convergence Status', color='green')
    ax2.set_ylim(0.5, 1.5)
    ax2.set_yticks([1])
    ax2.set_yticklabels(['Convergence'])
    
    ax.set_xlabel('Axial Position (m)')
    ax.set_ylabel('Number of Halvings')
    ax.grid(True, alpha=0.3)
    ax.set_title('Chemistry Convergence')
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')

def create_comprehensive_plot(results: Dict[str, Any], save_path: Optional[str] = None, show: bool = True):
    """Create comprehensive adaptive chemistry visualization."""
    data = results['axial']
    kpis = results['kpis']
    chem_data = results['chemistry']
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    # Temperature profile
    ax1 = fig.add_subplot(gs[0, 0])
    plot_temperature_profile(ax1, data, kpis)
    
    # Species profiles
    ax2 = fig.add_subplot(gs[0, 1])
    plot_species_profiles(ax2, data)
    
    # H₂/CO ratio
    ax3 = fig.add_subplot(gs[1, 0])
    plot_h2_co_ratio(ax3, data, kpis)
    
    # Chemistry fallback map
    ax4 = fig.add_subplot(gs[1, 1])
    plot_chemistry_fallback_map(ax4, chem_data)
    
    # Adaptive step sizes
    ax5 = fig.add_subplot(gs[2, 0])
    plot_adaptive_step_sizes(ax5, chem_data)
    
    # Chemistry convergence
    ax6 = fig.add_subplot(gs[2, 1])
    plot_chemistry_convergence(ax6, chem_data)
    
    # Add title
    case_name = results['case_path'].name
    fig.suptitle(f'Adaptive Chemistry PFR Results: {case_name}', fontsize=16, fontweight='bold')
    
    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Comprehensive plot saved to: {save_path}")
    
    # Show if requested
    if show:
        plt.show()
    else:
        plt.close()

def print_chemistry_summary(results: Dict[str, Any]):
    """Print chemistry diagnostics summary."""
    chem_data = results['chemistry']
    kpis = results['kpis']
    
    print("\n" + "="*60)
    print("ADAPTIVE CHEMISTRY SUMMARY")
    print("="*60)
    
    if chem_data is not None:
        total_cells = len(chem_data)
        fallback_cells = chem_data['fallback_used'].sum()
        converged_cells = chem_data['converged'].sum()
        
        print(f"Total Chemistry Cells: {total_cells}")
        print(f"Cells Using Fallback: {fallback_cells} ({100*fallback_cells/total_cells:.1f}%)")
        print(f"Converged Cells: {converged_cells} ({100*converged_cells/total_cells:.1f}%)")
        print(f"Average Substeps: {chem_data['substeps'].mean():.1f}")
        print(f"Average Halvings: {chem_data['halvings'].mean():.1f}")
        print(f"Min Step Size: {chem_data['dt_min'].min():.2e} s")
        print(f"Max Step Size: {chem_data['dt_max'].max():.2e} s")
    else:
        print("No chemistry diagnostics available")
    
    # KPIs
    if kpis:
        print(f"\nKey Performance Indicators:")
        if 'H2_CO_ratio' in kpis:
            print(f"H₂/CO Ratio: {kpis['H2_CO_ratio']:.2f}")
        if 'z_ign_peak_m' in kpis:
            print(f"Ignition Length: {kpis['z_ign_peak_m']:.3f} m")
        if 'T_out_C' in kpis:
            print(f"Outlet Temperature: {kpis['T_out_C']:.1f} °C")
    
    print("="*60)

def main():
    """Main plotting function."""
    parser = argparse.ArgumentParser(description='Plot adaptive chemistry PFR results')
    parser.add_argument('--case', type=str, required=True, 
                       help='Case directory path (e.g., hp_pox_results/A_case1_richter)')
    parser.add_argument('--save', type=str, default=None,
                       help='Save plot to file (optional)')
    parser.add_argument('--show', action='store_true', default=True,
                       help='Show plot (default: True)')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not show plot')
    
    args = parser.parse_args()
    
    # Load results
    try:
        results = load_results(args.case)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1
    
    # Print summary
    print_chemistry_summary(results)
    
    # Create plot
    save_path = args.save
    if save_path is None:
        save_path = Path(args.case) / "adaptive_chemistry_comprehensive.png"
    
    show_plot = args.show and not args.no_show
    create_comprehensive_plot(results, save_path, show_plot)
    
    return 0

if __name__ == "__main__":
    exit(main())
