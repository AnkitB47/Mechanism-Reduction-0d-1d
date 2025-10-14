#!/usr/bin/env python3
"""
Baseline simulation tool for HP-POX model with robust PSR.
Generates fresh baseline results using the updated PSR robustness logic.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import json
from typing import List, Dict, Any

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hp_pox_physics.hp_pox_model import HPPOXModel
from tools.plot_axial_unified import plot_temperature, plot_species, plot_h2co, plot_species_vs_time, plot_velocity_pressure
from hp_pox_physics.config_cases import richter_case_1, richter_case_4


def simulate_baseline_cases(base_mech_path: str, out_root: Path, cases_root: Path, cases: List[str]):
    """
    Simulate baseline cases with robust PSR and generate full figure set.
    
    Args:
        base_mech_path: Path to base mechanism file
        out_root: Output root directory for baseline results
        cases_root: Root directory containing case configurations
        cases: List of case names to simulate
    """
    print("=" * 80)
    print("HP-POX BASELINE SIMULATION (GRI v2)")
    print("=" * 80)
    print(f"Base mechanism: {base_mech_path}")
    print(f"Output root: {out_root}")
    print(f"Cases root: {cases_root}")
    print(f"Cases: {', '.join(cases)}")
    print()

    # Initialize model with robust PSR settings
    model = HPPOXModel(mechanism_path=base_mech_path)
    print(f"Loaded mechanism: {len(model.thermo.species_names)} species, {len(model.thermo.gas.reaction_equations())} reactions")

    # Clean previous artifacts for specified cases
    print("Cleaning previous artifacts...")
    for case_id in cases:
        case_dir = out_root / case_id
        figs_dir = case_dir / 'figs'

        # Remove figs directory if it exists
        if figs_dir.exists():
            import shutil
            shutil.rmtree(figs_dir)
            print(f"  Removed figs directory: {figs_dir}")

        # Remove data files
        for filename in ['axial_profiles.csv', 'kpis.json', 'diagnostics.json']:
            file_path = case_dir / filename
            if file_path.exists():
                file_path.unlink()
                print(f"  Removed: {file_path}")

    print("Clean complete.")
    print()
    
    # Create output directories
    out_root.mkdir(parents=True, exist_ok=True)
    
    results_summary = []
    
    for case_id in cases:
        print(f"\n{'='*60}")
        print(f"SIMULATING CASE: {case_id}")
        print(f"{'='*60}")
        
        try:
            # Get case configuration using new dataclass system
            if 'case1' in case_id.lower():
                case_config = richter_case_1
            elif 'case4' in case_id.lower():
                case_config = richter_case_4
            else:
                raise ValueError(f"Unknown case: {case_id}")

            # Run the case with robust PSR
            case_out_dir = out_root / case_id
            case_out_dir.mkdir(parents=True, exist_ok=True)

            # Create inlet configuration (simplified for this validation)
            inlet_config = {
                'case1' if 'case1' in case_id.lower() else 'case4': {
                    'primary_steam': {
                        'mass_kg_per_h': case_config.primary_steam.mass_kg_per_s * 3600,
                        'T_C': case_config.primary_steam.T_K - 273.15,
                        'comp': {'H2O': 1.0}
                    },
                    'secondary_steam': {
                        'mass_kg_per_h': case_config.secondary_steam.mass_kg_per_s * 3600,
                        'T_C': case_config.secondary_steam.T_K - 273.15,
                        'comp': {'H2O': 1.0}
                    },
                    'oxygen': {
                        'mass_kg_per_h': case_config.oxygen.mass_kg_per_s * 3600,
                        'T_C': case_config.oxygen.T_K - 273.15,
                        'comp': {'O2': 1.0}
                    },
                    'nitrogen_optisox': {
                        'mass_kg_per_h': case_config.nitrogen.mass_kg_per_s * 3600,
                        'T_C': case_config.nitrogen.T_K - 273.15,
                        'comp': {'N2': 1.0}
                    },
                    'natural_gas': {
                        'mass_kg_per_h': case_config.natural_gas.mass_kg_per_s * 3600,
                        'T_C': case_config.natural_gas.T_K - 273.15,
                        'composition_volpct': case_config.natural_gas.composition_X
                    }
                }
            }

            results = model.run_case(case_id, inlet_config, case_out_dir)
            
            # Extract key metrics
            case_summary = {
                'case': case_id,
                'psr_converged': results['psr_diagnostics']['converged'],
                'psr_t_out': results['psr_diagnostics']['T_out'],
                'psr_branch': results['psr_diagnostics'].get('branch_status', 'UNKNOWN'),
                'ch4_slip_pct': results.get('sim_dry', {}).get('CH4', 0.0) * 100,
                'h2_co_ratio': results.get('h2_co_ratio', 0.0),
                't_out_avg': results.get('T_out_avg', 0.0)
            }
            results_summary.append(case_summary)
            
            # Print case summary
            print(f"\nCase {case_id} Results:")
            print(f"  PSR Converged: {case_summary['psr_converged']}")
            print(f"  PSR T_out: {case_summary['psr_t_out']:.1f}K")
            print(f"  PSR Branch: {case_summary['psr_branch']}")
            print(f"  CH4 Slip: {case_summary['ch4_slip_pct']:.2f}%")
            print(f"  H2/CO Ratio: {case_summary['h2_co_ratio']:.3f}")
            print(f"  T_out Avg: {case_summary['t_out_avg']:.1f}K")
            
            # Generate comprehensive figure set
            print(f"\nGenerating figures for {case_id}...")
            figs_dir = case_out_dir / 'figs'
            figs_dir.mkdir(exist_ok=True)
            
            # Load axial profiles
            axial_df = pd.read_csv(case_out_dir / 'axial_profiles.csv')
            
            # Generate all required figures using the correct plotting functions
            try:
                # Temperature vs z
                plot_temperature(axial_df, figs_dir)
                print(f"    PASS: T_vs_z.png")
            except Exception as e:
                print(f"    FAIL: T_vs_z.png failed: {e}")
            
            try:
                # Species vs z
                plot_species(axial_df, figs_dir)
                print(f"    PASS: species_vs_z.png")
            except Exception as e:
                print(f"    FAIL: species_vs_z.png failed: {e}")
            
            try:
                # Species vs time
                plot_species_vs_time(axial_df, figs_dir)
                print(f"    PASS: species_vs_time.png")
            except Exception as e:
                print(f"    FAIL: species_vs_time.png failed: {e}")
            
            try:
                # H2CO vs z
                plot_h2co(axial_df, figs_dir)
                print(f"    PASS: H2CO_vs_z.png")
            except Exception as e:
                print(f"    FAIL: H2CO_vs_z.png failed: {e}")
            
            try:
                # Velocity and pressure vs z
                plot_velocity_pressure(axial_df, figs_dir)
                print(f"    PASS: u_vs_z.png and u_p_vs_z.png")
            except Exception as e:
                print(f"    FAIL: u_vs_z.png and u_p_vs_z.png failed: {e}")
            
            print(f"  PASS: Case {case_id} completed successfully")
            
        except Exception as e:
            print(f"  FAIL: Case {case_id} failed: {e}")
            case_summary = {
                'case': case_id,
                'psr_converged': False,
                'psr_t_out': 0.0,
                'psr_branch': 'FAILED',
                'ch4_slip_pct': 0.0,
                'h2_co_ratio': 0.0,
                't_out_avg': 0.0,
                'error': str(e)
            }
            results_summary.append(case_summary)
            continue
    
    # Print final summary
    print(f"\n{'='*80}")
    print("BASELINE SIMULATION SUMMARY")
    print(f"{'='*80}")
    
    for summary in results_summary:
        status = "PASS" if summary['psr_converged'] else "FAIL"
        ch4_warning = " ⚠️" if summary['ch4_slip_pct'] > 10.0 else ""
        print(f"{summary['case']:20} | PSR: {status:4} | CH4: {summary['ch4_slip_pct']:5.2f}%{ch4_warning} | H2/CO: {summary['h2_co_ratio']:5.3f} | T_out: {summary['psr_t_out']:6.1f}K")
    
    # Save summary to file
    summary_file = out_root / 'baseline_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(results_summary, f, indent=2)
    
    print(f"\nBaseline simulation complete. Results saved to: {out_root}")
    print(f"Summary saved to: {summary_file}")


def main():
    """Main entry point for baseline simulation."""
    parser = argparse.ArgumentParser(description='Simulate HP-POX baseline cases with robust PSR')
    parser.add_argument('--base-mech', required=True, help='Path to base mechanism file')
    parser.add_argument('--out-root', required=True, help='Output root directory for baseline results')
    parser.add_argument('--cases-root', required=True, help='Root directory containing case configurations')
    parser.add_argument('--cases', nargs='+', required=True, help='Case names to simulate')
    parser.add_argument('--overwrite-figs', action='store_true', help='Overwrite existing figure files')
    parser.add_argument('--recompute-kpis', action='store_true', help='Recompute KPIs from PFR output')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    base_mech_path = Path(args.base_mech)
    out_root = Path(args.out_root)
    cases_root = Path(args.cases_root)
    cases = args.cases
    
    if not base_mech_path.exists():
        raise FileNotFoundError(f"Base mechanism not found: {base_mech_path}")
    
    if not cases_root.exists():
        raise FileNotFoundError(f"Cases root not found: {cases_root}")
    
    simulate_baseline_cases(str(base_mech_path), out_root, cases_root, cases)


if __name__ == '__main__':
    main()
