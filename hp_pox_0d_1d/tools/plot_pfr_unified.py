#!/usr/bin/env python3
"""
Unified PFR plotting script for physics + diagnostics data.
Generates 6 comprehensive plots from unified CSV files.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import sys

def load_unified_csv(csv_path: str) -> pd.DataFrame:
    """Load unified CSV file with validation."""
    try:
        df = pd.read_csv(csv_path)
        
        # Validate required columns
        required_cols = ['z', 'T', 'used_fallback', 'substeps', 'halvings', 'dt_min', 'dt_effective']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        print(f"Loaded unified CSV: {len(df)} rows, {len(df.columns)} columns")
        return df
        
    except Exception as e:
        print(f"Error loading CSV {csv_path}: {e}")
        sys.exit(1)

def plot_temperature_vs_z(df: pd.DataFrame, save_path: str):
    """Plot temperature vs axial position."""
    plt.figure(figsize=(10, 6))
    plt.plot(df['z'], df['T'], 'b-', linewidth=2, label='Temperature')
    plt.xlabel('Axial Position (m)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature Profile')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Temperature plot saved: {save_path}")

def plot_species_vs_z(df: pd.DataFrame, save_path: str):
    """Plot species mole fractions vs axial position."""
    plt.figure(figsize=(10, 6))
    
    # Available species
    species_names = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'O2', 'N2']
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'gray']
    
    for species, color in zip(species_names, colors):
        if species in df.columns and df[species].max() > 0:
            plt.plot(df['z'], df[species] * 100, color=color, linewidth=2, label=species)
    
    plt.xlabel('Axial Position (m)')
    plt.ylabel('Mole Fraction (%)')
    plt.title('Species Profiles')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Species plot saved: {save_path}")

def plot_h2co_vs_z(df: pd.DataFrame, save_path: str):
    """Plot H2/CO ratio vs axial position."""
    plt.figure(figsize=(10, 6))
    
    if 'H2_CO' in df.columns and df['H2_CO'].max() > 0:
        plt.plot(df['z'], df['H2_CO'], 'r-', linewidth=2, label='H₂/CO Ratio')
        
        # Add target line
        plt.axhline(2.0, color='black', linestyle='--', alpha=0.7, label='Target: 2.0')
        
        # Add outlet value annotation
        if len(df) > 0:
            outlet_ratio = df['H2_CO'].iloc[-1]
            plt.text(0.02, 0.98, f'Outlet H₂/CO: {outlet_ratio:.2f}', 
                    transform=plt.gca().transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('Axial Position (m)')
    plt.ylabel('H₂/CO Ratio')
    plt.title('H₂/CO Ratio Profile')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  H2/CO ratio plot saved: {save_path}")

def plot_fallback_map(df: pd.DataFrame, save_path: str):
    """Plot fallback usage map."""
    plt.figure(figsize=(10, 6))
    
    z = df['z']
    fallback_used = df['used_fallback']
    
    # Create stem plot
    fallback_mask = fallback_used == 1
    plt.stem(z[fallback_mask], np.ones(fallback_mask.sum()), 
             basefmt=' ', linefmt='red-', markerfmt='ro', label='Fallback Used')
    
    plt.xlabel('Axial Position (m)')
    plt.ylabel('Fallback Usage')
    plt.title('Chemistry Fallback Map')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.ylim(-0.1, 1.1)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Fallback map saved: {save_path}")

def plot_stepsizes(df: pd.DataFrame, save_path: str):
    """Plot step sizes vs axial position."""
    plt.figure(figsize=(10, 6))
    
    z = df['z']
    dt_min = df['dt_min']
    dt_effective = df['dt_effective']
    
    # Plot step size range
    plt.fill_between(z, dt_min, dt_effective, alpha=0.3, color='blue', label='Step Size Range')
    plt.plot(z, dt_min, 'b-', linewidth=1, label='dt_min')
    plt.plot(z, dt_effective, 'b--', linewidth=2, label='dt_effective')
    
    plt.xlabel('Axial Position (m)')
    plt.ylabel('Time Step (s)')
    plt.title('Adaptive Step Sizes')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Step sizes plot saved: {save_path}")

def plot_convergence(df: pd.DataFrame, save_path: str):
    """Plot convergence metrics vs axial position."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    z = df['z']
    substeps = df['substeps']
    halvings = df['halvings']
    
    # Plot substeps
    ax1.plot(z, substeps, 'b-', linewidth=2, label='Substeps')
    ax1.set_ylabel('Number of Substeps')
    ax1.set_title('Chemistry Convergence Metrics')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot halvings
    ax2.plot(z, halvings, 'r-', linewidth=2, label='Halvings')
    ax2.set_xlabel('Axial Position (m)')
    ax2.set_ylabel('Number of Halvings')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Convergence plot saved: {save_path}")

def create_all_plots(csv_path: str, output_dir: str):
    """Create all 6 plots for a unified CSV file."""
    # Load data
    df = load_unified_csv(csv_path)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    plot_temperature_vs_z(df, str(output_path / "T_vs_z.png"))
    plot_species_vs_z(df, str(output_path / "species_vs_z.png"))
    plot_h2co_vs_z(df, str(output_path / "H2CO_vs_z.png"))
    plot_fallback_map(df, str(output_path / "fallback_map.png"))
    plot_stepsizes(df, str(output_path / "stepsizes.png"))
    plot_convergence(df, str(output_path / "convergence.png"))
    
    print(f"All plots saved to: {output_path}")

def main():
    """Main plotting function."""
    parser = argparse.ArgumentParser(description='Plot unified PFR results')
    parser.add_argument('csv_file', type=str, help='Unified CSV file path')
    parser.add_argument('--out', type=str, default='figures', help='Output directory (default: figures)')
    
    args = parser.parse_args()
    
    # Extract basename for output directory
    csv_basename = Path(args.csv_file).stem
    output_dir = Path(args.out) / csv_basename
    
    print(f"Creating plots for: {args.csv_file}")
    print(f"Output directory: {output_dir}")
    
    try:
        create_all_plots(args.csv_file, str(output_dir))
        print("✅ All plots generated successfully!")
    except Exception as e:
        print(f"❌ Error generating plots: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
