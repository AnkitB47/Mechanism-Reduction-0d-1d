#!/usr/bin/env python3
"""
PFR parameter sweep utility to test CH4 slip reduction through process parameter changes.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import json

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hp_pox_physics.hp_pox_model import HPPOXModel


def run_parameter_sweep():
    """Run PFR parameter sweep for CH4 slip reduction."""
    
    print("Running PFR parameter sweep for CH4 slip reduction...")
    
    # Load base model
    model = HPPOXModel("gri30.yaml", transport_model='mixture-averaged')
    
    cases = ['A_case1_richter', 'A_case4_richter', 'B_case1_134L', 'B_case4_134L']
    
    # Parameter perturbations
    perturbations = [
        {'name': 'baseline', 'o2_factor': 1.0, 'mdot_factor': 1.0, 'description': 'Original parameters'},
        {'name': 'o2_plus_2pct', 'o2_factor': 1.02, 'mdot_factor': 1.0, 'description': 'O2 +2% relative'},
        {'name': 'o2_plus_5pct', 'o2_factor': 1.05, 'mdot_factor': 1.0, 'description': 'O2 +5% relative'},
        {'name': 'mdot_minus_10pct', 'o2_factor': 1.0, 'mdot_factor': 0.9, 'description': 'Mass flow -10%'},
        {'name': 'mdot_minus_20pct', 'o2_factor': 1.0, 'mdot_factor': 0.8, 'description': 'Mass flow -20%'},
    ]
    
    results = []
    
    for case in cases:
        print(f"\nProcessing case: {case}")
        
        for pert in perturbations:
            print(f"  Testing: {pert['description']}")
            
            try:
                # Modify inlet streams based on perturbation
                inlet_streams = model.load_inlet_streams('config.yaml')
                
                # Apply O2 perturbation (increase oxygen flow rate)
                if pert['o2_factor'] != 1.0:
                    if 'oxygen' in inlet_streams:
                        inlet_streams['oxygen']['mass_kg_per_h'] *= pert['o2_factor']
                
                # Apply mass flow perturbation (affects residence time)
                if pert['mdot_factor'] != 1.0:
                    for stream_name in inlet_streams:
                        if stream_name != 'geometry':  # Don't modify geometry
                            inlet_streams[stream_name]['mass_kg_per_h'] *= pert['mdot_factor']
                
                # Run PFR with chemistry ON only
                results_dict = model.run_case(case, inlet_streams, Path('outputs/sweep'), chemistry=True)
                
                # Extract KPIs
                h2_out = results_dict['dry_composition'].get('H2', 0.0)
                co_out = results_dict['dry_composition'].get('CO', 0.0)
                ch4_out = results_dict['dry_composition'].get('CH4', 0.0)
                t_out_k = results_dict['outlet_state'].temperature
                h2co = h2_out / max(co_out, 1e-12)
                
                result_row = {
                    'case': case,
                    'perturbation': pert['name'],
                    'description': pert['description'],
                    'o2_factor': pert['o2_factor'],
                    'mdot_factor': pert['mdot_factor'],
                    'H2_out': h2_out,
                    'CO_out': co_out,
                    'CH4_out': ch4_out,
                    'CH4_slip_pct': ch4_out * 100.0,
                    'H2CO_ratio': h2co,
                    'T_out_K': t_out_k,
                    'success': True
                }
                
                print(f"    CH4 slip: {ch4_out*100:.2f}%, H2/CO: {h2co:.3f}")
                
            except Exception as e:
                print(f"    Failed: {e}")
                result_row = {
                    'case': case,
                    'perturbation': pert['name'],
                    'description': pert['description'],
                    'o2_factor': pert['o2_factor'],
                    'mdot_factor': pert['mdot_factor'],
                    'H2_out': 0.0,
                    'CO_out': 0.0,
                    'CH4_out': 0.0,
                    'CH4_slip_pct': 0.0,
                    'H2CO_ratio': 0.0,
                    'T_out_K': 0.0,
                    'success': False
                }
            
            results.append(result_row)
    
    # Create results DataFrame
    df_results = pd.DataFrame(results)
    
    # Save results
    output_file = Path('pfr_sweep_results.csv')
    df_results.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Analyze results
    print("\n" + "="*80)
    print("PFR PARAMETER SWEEP ANALYSIS")
    print("="*80)
    
    print("\nCH4 Slip Summary:")
    print("-" * 50)
    pivot_table = df_results.pivot_table(
        values='CH4_slip_pct', 
        index='case', 
        columns='perturbation', 
        aggfunc='mean'
    )
    print(pivot_table.round(2))
    
    print("\nH2/CO Ratio Summary:")
    print("-" * 50)
    h2co_pivot = df_results.pivot_table(
        values='H2CO_ratio', 
        index='case', 
        columns='perturbation', 
        aggfunc='mean'
    )
    print(h2co_pivot.round(3))
    
    # Find perturbations that achieve CH4 slip < 5%
    print("\nPerturbations achieving CH4 slip < 5%:")
    print("-" * 50)
    successful_perturbations = df_results[df_results['CH4_slip_pct'] < 5.0]
    
    if len(successful_perturbations) > 0:
        for _, row in successful_perturbations.iterrows():
            print(f"{row['case']} - {row['description']}: {row['CH4_slip_pct']:.2f}% CH4 slip, H2/CO = {row['H2CO_ratio']:.3f}")
    else:
        print("No perturbations achieved CH4 slip < 5%")
    
    # Find best perturbation for each case
    print("\nBest perturbation per case (lowest CH4 slip):")
    print("-" * 50)
    for case in cases:
        case_results = df_results[df_results['case'] == case]
        best_idx = case_results['CH4_slip_pct'].idxmin()
        best_row = case_results.loc[best_idx]
        print(f"{case}: {best_row['description']} - {best_row['CH4_slip_pct']:.2f}% CH4 slip")
    
    return df_results


if __name__ == "__main__":
    results = run_parameter_sweep()
