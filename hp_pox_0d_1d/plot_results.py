#!/usr/bin/env python3
"""
Plot results from HP-POX 0D-1D reactor model.

This script generates publication-ready plots from the simulation results.
"""

import argparse
import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import os
from datetime import datetime

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


def detect_species_columns(df: pd.DataFrame) -> List[str]:
    """
    Detect species columns in the dataframe by excluding known non-species columns.
    
    Args:
        df: Input dataframe
        
    Returns:
        List of species column names
    """
    exclude_cols = {
        'z_m', 'T_C', 'z', 'T', 'index', 'Unnamed: 0',
        'rho_kg_m3', 'u_m_s', 'residence_time_s', 'P_Pa', 'P_bar',
        'cp_J_kg_K', 'mu_Pa_s', 'k_W_m_K', 'Re', 'Pr', 'Nu'
    }
    
    # Look for species columns (typically X_* or Y_* patterns)
    species_cols = []
    for col in df.columns:
        if col not in exclude_cols:
            # Check if it looks like a species column
            if (col.startswith('X_') or col.startswith('Y_') or 
                col in ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'N2', 'O2', 'H', 'O', 'OH', 'HO2', 'H2O2',
                       'C', 'CH', 'CH2', 'CH2(S)', 'CH3', 'HCO', 'CH2O', 'CH2OH', 'CH3O', 'CH3OH',
                       'C2H', 'C2H2', 'C2H3', 'C2H4', 'C2H5', 'C2H6', 'HCCO', 'CH2CO', 'HCCOH',
                       'N', 'NH', 'NH2', 'NH3', 'NNH', 'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN',
                       'H2CN', 'HCNN', 'HCNO', 'HOCN', 'HNCO', 'NCO', 'AR', 'C3H7', 'C3H8',
                       'CH2CHO', 'CH3CHO']):
                species_cols.append(col)
    
    return species_cols


def extract_species_mapping(species_cols: List[str]) -> Dict[str, str]:
    """
    Extract species name mapping from column names.
    
    Args:
        species_cols: List of species column names
        
    Returns:
        Dictionary mapping species names to column names
    """
    species_mapping = {}
    
    for col in species_cols:
        if col.startswith('X_') or col.startswith('Y_'):
            # Extract species name by removing prefix
            species_name = col[2:]  # Remove 'X_' or 'Y_'
            species_mapping[species_name] = col
        elif col in ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'N2', 'O2', 'H', 'O', 'OH', 'HO2', 'H2O2',
                     'C', 'CH', 'CH2', 'CH2(S)', 'CH3', 'HCO', 'CH2O', 'CH2OH', 'CH3O', 'CH3OH',
                     'C2H', 'C2H2', 'C2H3', 'C2H4', 'C2H5', 'C2H6', 'HCCO', 'CH2CO', 'HCCOH',
                     'N', 'NH', 'NH2', 'NH3', 'NNH', 'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN',
                     'H2CN', 'HCNN', 'HCNO', 'HOCN', 'HNCO', 'NCO', 'AR', 'C3H7', 'C3H8',
                     'CH2CHO', 'CH3CHO']:
            species_mapping[col] = col
    
    return species_mapping


def compute_dry_basis_composition(df: pd.DataFrame, species_mapping: Dict[str, str]) -> Dict[str, float]:
    """
    Compute dry basis composition from the last row of axial profiles.
    
    Args:
        df: Axial profiles dataframe
        species_mapping: Dictionary mapping species names to column names
        
    Returns:
        Dictionary of species mole fractions on dry basis
    """
    if df.empty:
        return {}
    
    # Get the last row (outlet)
    outlet_row = df.iloc[-1]
    
    # Build species composition dict using mapping
    composition = {}
    for species_name, col_name in species_mapping.items():
        if col_name in df.columns:
            composition[species_name] = float(outlet_row[col_name])
        else:
            composition[species_name] = 0.0
    
    # Compute dry basis by excluding H2O
    h2o_fraction = 0.0
    if 'H2O' in composition:
        h2o_fraction = composition['H2O']
    
    dry_total = 1.0 - h2o_fraction
    if dry_total > 0:
        # Only normalize non-H2O species
        for species in composition:
            if species != 'H2O':
                composition[species] = composition[species] / dry_total
    # If no H2O or dry_total <= 0, treat as already dry
    
    return composition


def extract_target_values(kpis: Dict) -> Tuple[Dict[str, float], bool]:
    """
    Extract target values from KPIs with support for multiple key formats.
    
    Args:
        kpis: KPIs dictionary
        
    Returns:
        Tuple of (target_values_dict, targets_available)
    """
    target_values = {}
    targets_available = False
    
    # Try different possible key names
    target_keys = ['targets', 'outlet_target', 'target_dry_basis']
    
    for key in target_keys:
        if key in kpis and kpis[key] is not None:
            targets = kpis[key]
            
            # Handle both dict and list formats
            if isinstance(targets, dict):
                target_values = targets.copy()
                targets_available = True
                break
            elif isinstance(targets, list):
                # Assume list format with species names as keys
                # This would need to be adapted based on actual format
                target_values = {f'species_{i}': val for i, val in enumerate(targets)}
                targets_available = True
                break
    
    return target_values, targets_available


def plot_case(case_dir: Path, save: bool = True, show: bool = False, force: bool = False, debug: bool = False) -> None:
    """
    Generate all plots for a single case.
    
    Args:
        case_dir: Path to the case directory containing the data files
        save: Whether to save the plots to files
        show: Whether to display the plots in windows
        force: Whether to force regeneration of existing plots
        debug: Whether to print debug information
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
    if not axial_file.exists():
        print(f"Error: Required file {axial_file} does not exist, skipping...")
        return
    
    if not kpis_file.exists():
        print(f"Warning: KPIs file {kpis_file} does not exist, using defaults...")
        kpis = {}
    else:
        try:
            with open(kpis_file, 'r') as f:
                kpis = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load KPIs file: {e}, using defaults...")
            kpis = {}
    
    # Check if outlet composition file exists
    outlet_exists = outlet_file.exists()
    if not outlet_exists:
        print(f"Info: Outlet composition file {outlet_file} does not exist, will compute from axial profiles...")
    
    try:
        # Load axial profiles data
        df_axial = pd.read_csv(axial_file)
        if df_axial.empty:
            print(f"Error: Axial profiles file {axial_file} is empty, skipping...")
            return
        
        # Detect species columns and create mapping
        species_cols = detect_species_columns(df_axial)
        species_mapping = extract_species_mapping(species_cols)
        
        if debug:
            print(f"  DEBUG: Detected species columns: {species_cols}")
            print(f"  DEBUG: Species mapping: {species_mapping}")
        
        # Load outlet composition if available
        df_outlet = None
        if outlet_exists:
            try:
                df_outlet = pd.read_csv(outlet_file)
            except Exception as e:
                print(f"Warning: Could not load outlet composition file: {e}")
                df_outlet = None
        
        # Compute dry basis composition
        sim_dry = compute_dry_basis_composition(df_axial, species_mapping)
        target_dry, targets_available = extract_target_values(kpis)
        
        if debug:
            print(f"  DEBUG: outlet_composition_table.csv found: {outlet_exists}")
            print(f"  DEBUG: sim_dry: {sim_dry}")
            print(f"  DEBUG: target_dry: {target_dry}")
            print(f"  DEBUG: targets_available: {targets_available}")
        
        # Generate plots
        plot_log = {
            'case_name': case_name,
            'timestamp': datetime.now().isoformat(),
            'species_columns': species_cols,
            'species_mapping': species_mapping,
            'outlet_csv_existed': outlet_exists,
            'dry_basis_computed': not outlet_exists,
            'sim_dry': sim_dry,
            'target_dry': target_dry,
            'targets_available': targets_available,
            'files_used': {
                'axial_profiles': str(axial_file),
                'kpis': str(kpis_file),
                'outlet_composition': str(outlet_file) if outlet_exists else None
            }
        }
        
        plot_axial_profiles(df_axial, species_cols, case_dir, save, show, force)
        plot_ignition_markers(df_axial, kpis, case_dir, save, show, force)
        
        # Generate outlet composition plot
        plot_outlet_composition(df_axial, species_mapping, df_outlet, kpis, case_dir, save, show, force, debug)
        
        # Save plot log
        if save:
            log_file = case_dir / "plot_log.json"
            with open(log_file, 'w') as f:
                json.dump(plot_log, f, indent=2)
            print(f"  Saved: {log_file}")
        
        print(f"✓ Generated plots for {case_name}")
        
    except Exception as e:
        print(f"Error processing {case_name}: {e}")
        import traceback
        traceback.print_exc()
        return


def plot_axial_profiles(df: pd.DataFrame, species_cols: List[str], case_dir: Path, save: bool, show: bool, force: bool) -> None:
    """
    Generate 2-panel axial profiles plot.
    
    Top panel: Temperature vs axial position
    Bottom panel: Species mole fractions vs axial position
    """
    # Check if we should skip due to existing file
    output_file = case_dir / "axial_profiles_plot.png"
    if not force and output_file.exists():
        print(f"  Skipping {output_file} (already exists, use --force to overwrite)")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Top panel: Temperature
    # Handle different temperature column names
    if 'z_m' in df.columns:
        z_m = df['z_m']
    elif 'z' in df.columns:
        z_m = df['z']
    else:
        print("Error: No axial position column found (z_m or z)")
        return
    
    if 'T_C' in df.columns:
        T_C = df['T_C']  # Already in Celsius
    elif 'T' in df.columns:
        T_C = df['T'] - 273.15  # Convert from K to C
    else:
        print("Error: No temperature column found (T_C or T)")
        return
    
    ax1.plot(z_m, T_C, 'k-', linewidth=2, label='Temperature')
    ax1.set_xlabel('Axial Position (m)')
    ax1.set_ylabel('Temperature (°C)')
    ax1.set_title('Temperature Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Bottom panel: Species mole fractions
    # Select species to plot
    target_species = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'N2']
    available_species = [s for s in target_species if s in species_cols]
    
    # If fewer than 3 target species available, use top 6 by max mole fraction
    if len(available_species) < 3:
        species_max_values = {}
        for species in species_cols:
            if species in df.columns:
                species_max_values[species] = df[species].max()
        
        # Sort by max value and take top 6
        sorted_species = sorted(species_max_values.items(), key=lambda x: x[1], reverse=True)
        available_species = [s[0] for s in sorted_species[:6]]
        print(f"  Using top species by max mole fraction: {available_species}")
    
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    
    for i, species in enumerate(available_species):
        if species in df.columns:
            color = colors[i % len(colors)]
            # Convert mole fraction to percentage
            mole_frac_pct = df[species] * 100
            ax2.plot(z_m, mole_frac_pct, color=color, linewidth=2, label=species)
    
    ax2.set_xlabel('Axial Position (m)')
    ax2.set_ylabel('Mole Fraction (%)')
    ax2.set_title('Species Composition Profile')
    ax2.grid(True, alpha=0.3)
    # Place legend outside the plot
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(output_file)
        print(f"  Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_ignition_markers(df: pd.DataFrame, kpis: Dict, case_dir: Path, save: bool, show: bool, force: bool) -> None:
    """
    Generate temperature profile with ignition length markers.
    """
    # Check if we should skip due to existing file
    output_file = case_dir / "ignition_length_markers.png"
    if not force and output_file.exists():
        print(f"  Skipping {output_file} (already exists, use --force to overwrite)")
        return
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot temperature profile
    # Handle different column names
    if 'z_m' in df.columns:
        z_m = df['z_m']
    elif 'z' in df.columns:
        z_m = df['z']
    else:
        print("Error: No axial position column found (z_m or z)")
        return
    
    if 'T_C' in df.columns:
        T_C = df['T_C']  # Already in Celsius
    elif 'T' in df.columns:
        T_C = df['T'] - 273.15  # Convert from K to C
    else:
        print("Error: No temperature column found (T_C or T)")
        return
    
    ax.plot(z_m, T_C, 'k-', linewidth=2, label='Temperature')
    
    # Add ignition markers
    z_ign_peak = kpis.get('z_ign_peak_m', kpis.get('ignition_length_peak_m', 0.0))
    z_ign_thresh = kpis.get('z_ign_thresh_m', kpis.get('ignition_length_threshold_m', 0.0))
    
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
        plt.savefig(output_file)
        print(f"  Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_outlet_composition(df_axial: pd.DataFrame, species_mapping: Dict[str, str], df_outlet: Optional[pd.DataFrame], 
                           kpis: Dict, case_dir: Path, save: bool, show: bool, force: bool, debug: bool = False) -> None:
    """
    Generate outlet composition bar chart comparing simulated vs target.
    """
    # Check if we should skip due to existing file
    output_file = case_dir / "outlet_barplot.png"
    if not force and output_file.exists():
        print(f"  Skipping {output_file} (already exists, use --force to overwrite)")
        return
    
    # Target species for plotting (always show all 5)
    target_species = ['H2', 'CO', 'CO2', 'CH4', 'N2']
    
    # Get simulated values
    simulated_values = {}
    if df_outlet is not None:
        # Use outlet composition table if available
        if 'Species' in df_outlet.columns and 'Dry_Basis_Vol%' in df_outlet.columns:
            for species in target_species:
                species_data = df_outlet[df_outlet['Species'] == species]
                if not species_data.empty:
                    simulated_values[species] = float(species_data['Dry_Basis_Vol%'].iloc[0])
                else:
                    simulated_values[species] = 0.0
        else:
            if debug:
                print("  DEBUG: Outlet composition table has unexpected format, computing from axial profiles...")
            df_outlet = None
    
    if df_outlet is None:
        # Compute from axial profiles using species mapping
        composition = compute_dry_basis_composition(df_axial, species_mapping)
        for species in target_species:
            simulated_values[species] = composition.get(species, 0.0) * 100  # Convert to percentage
    
    # Get target values from KPIs
    target_values, targets_available = extract_target_values(kpis)
    
    # Ensure we have values for all target species (fill with 0.0 if missing)
    for species in target_species:
        if species not in simulated_values:
            simulated_values[species] = 0.0
        if species not in target_values:
            target_values[species] = 0.0
    
    if debug:
        print(f"  DEBUG: Final sim_dry: {simulated_values}")
        print(f"  DEBUG: Final target_dry: {target_values}")
        print(f"  DEBUG: targets_available: {targets_available}")
    
    # Create the bar chart
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    x = np.arange(len(target_species))
    width = 0.35
    
    # Prepare data arrays
    sim_values = [simulated_values.get(species, 0.0) for species in target_species]
    targ_values = [target_values.get(species, 0.0) for species in target_species] if targets_available else []
    
    # Plot bars
    bars1 = ax.bar(x - width/2, sim_values, width, label='Simulated', alpha=0.8, color='blue')
    
    if targets_available:
        bars2 = ax.bar(x + width/2, targ_values, width, label='Target', alpha=0.8, color='orange')
    else:
        bars2 = None
    
    # Add value labels on bars (always show labels, even for 0.0)
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax.annotate(f'{height:.1f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3),  # 3 points vertical offset
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)
    
    if bars2 is not None:
        for i, bar in enumerate(bars2):
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
    
    if targets_available:
        ax.legend()
    else:
        ax.legend([bars1[0]], ['Simulated (Target n/a)'])
    
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(output_file)
        print(f"  Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def run_unit_tests():
    """Run lightweight unit tests for species mapping and dry-basis computation."""
    print("Running unit tests...")
    
    # Test 1: Species mapping with X_* columns
    test_cols = ['z_m', 'T_C', 'X_H2', 'X_CO', 'X_CO2', 'X_CH4', 'X_H2O', 'X_N2', 'rho_kg_m3']
    species_cols = detect_species_columns(pd.DataFrame(columns=test_cols))
    species_mapping = extract_species_mapping(species_cols)
    
    expected_species = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'N2']
    for species in expected_species:
        assert species in species_mapping, f"Species {species} not found in mapping"
        assert species_mapping[species] == f'X_{species}', f"Wrong mapping for {species}"
    
    print("✓ Species mapping test passed")
    
    # Test 2: Dry-basis computation
    test_data = {
        'z_m': [0.0, 1.0, 2.0],
        'T_C': [100, 200, 300],
        'X_H2': [0.1, 0.2, 0.2],
        'X_CO': [0.1, 0.2, 0.2],
        'X_CO2': [0.1, 0.2, 0.2],
        'X_CH4': [0.1, 0.2, 0.2],
        'X_H2O': [0.4, 0.2, 0.2],
        'X_N2': [0.2, 0.0, 0.0]
    }
    
    df_test = pd.DataFrame(test_data)
    composition = compute_dry_basis_composition(df_test, species_mapping)
    
    # Debug: print the composition
    print(f"  DEBUG: Composition: {composition}")
    
    # Check that dry basis sums to 1.0 (excluding H2O)
    dry_sum = sum(composition[species] for species in expected_species if species != 'H2O')
    print(f"  DEBUG: Dry sum: {dry_sum}")
    # Allow for small numerical errors
    assert abs(dry_sum - 1.0) < 1e-3, f"Dry basis sum is {dry_sum}, should be 1.0"
    
    print("✓ Dry-basis computation test passed")
    print("All unit tests passed!")


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(description='Plot HP-POX 0D-1D reactor results')
    parser.add_argument('--case', type=str, help='Case ID to plot (e.g., A_case1_richter) or "all" for all cases')
    parser.add_argument('--save', action='store_true', default=True, help='Save figures (default: True)')
    parser.add_argument('--show', action='store_true', help='Display figures in windows')
    parser.add_argument('--force', action='store_true', help='Force regeneration of existing plots')
    parser.add_argument('--debug', action='store_true', help='Print debug information and run unit tests')
    parser.add_argument('--results-dir', type=str, default='hp_pox_results', 
                       help='Results directory (default: hp_pox_results)')
    
    args = parser.parse_args()
    
    # Run unit tests if debug mode
    if args.debug:
        run_unit_tests()
        print()
    
    results_dir = Path(args.results_dir)
    
    if not results_dir.exists():
        print(f"Error: Results directory {results_dir} does not exist")
        return 1
    
    if args.case == "all" or args.case is None:
        # Plot all cases that have axial_profiles.csv
        case_dirs = []
        for d in results_dir.iterdir():
            if d.is_dir() and (d / "axial_profiles.csv").exists():
                case_dirs.append(d)
        
        if not case_dirs:
            print(f"No case directories with axial_profiles.csv found in {results_dir}")
            return 1
        
        print(f"Found {len(case_dirs)} cases to process")
        
        for case_dir in sorted(case_dirs):
            plot_case(case_dir, save=args.save, show=args.show, force=args.force, debug=args.debug)
    
    elif args.case:
        # Plot specific case
        case_dir = results_dir / args.case
        plot_case(case_dir, save=args.save, show=args.show, force=args.force, debug=args.debug)
    else:
        print("Error: Please specify --case <case_id> or --case all")
        return 1
    
    print("Plotting complete!")
    return 0


if __name__ == "__main__":
    main()
