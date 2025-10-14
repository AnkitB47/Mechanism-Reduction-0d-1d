#!/usr/bin/env python3
"""
Adiabatic bound diagnostic tool.
Computes equilibrium adiabatic temperature for each case to establish physical bounds.
"""

import sys
import os
import json
import numpy as np
import cantera as ct
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hp_pox_physics.thermo import ThermoManager
from hp_pox_physics.hp_pox_model import HPPOXModel
from hp_pox_physics.config_cases import CaseConstants

def compute_adiabatic_flame_temperature(gas, P, T_in, Y_in):
    """
    Compute adiabatic flame temperature at constant pressure.
    
    Args:
        gas: Cantera gas object
        P: Pressure (Pa)
        T_in: Inlet temperature (K)
        Y_in: Inlet mass fractions
    
    Returns:
        T_ad: Adiabatic flame temperature (K)
        H_in: Inlet enthalpy (J/kg)
        H_out: Outlet enthalpy (J/kg)
    """
    # Set inlet state
    gas.TPY = T_in, P, Y_in
    H_in = gas.enthalpy_mass
    
    # Compute equilibrium at constant pressure and enthalpy
    gas.equilibrate('HP')
    T_ad = gas.T
    H_out = gas.enthalpy_mass
    
    # Verify enthalpy conservation
    H_error = abs(H_out - H_in) / abs(H_in) if abs(H_in) > 1e-9 else 0.0
    
    return T_ad, H_in, H_out, H_error

def assert_physical_bounds(T_out, T_ad, case_name):
    """
    Assert that outlet temperature is within physical bounds.
    
    Args:
        T_out: Actual outlet temperature (K)
        T_ad: Adiabatic flame temperature (K)
        case_name: Case name for error message
    
    Raises:
        AssertionError: If T_out > T_ad + 30K
    """
    max_allowed = T_ad + 30.0
    
    if T_out > max_allowed:
        raise AssertionError(
            f"PHYSICS_BOUND_VIOLATION in {case_name}: "
            f"T_out={T_out:.1f}K > T_ad+30K={max_allowed:.1f}K "
            f"(T_ad={T_ad:.1f}K)"
        )

def compute_elemental_balance(gas, Y_in, Y_out):
    """
    Compute elemental balance between inlet and outlet.
    
    Args:
        gas: Cantera gas object
        Y_in: Inlet mass fractions
        Y_out: Outlet mass fractions
    
    Returns:
        dict: Elemental balance errors for C, H, O, N
    """
    # Get molecular weights
    MW = gas.molecular_weights
    
    # Compute element mass fractions
    element_names = ['C', 'H', 'O', 'N']
    balance_errors = {}
    
    for element in element_names:
        if element in gas.element_names:
            idx = gas.element_names.index(element)
            
            # Compute mass flow of element
            mass_in = np.sum(Y_in * gas.elemental_mass_fraction(element))
            mass_out = np.sum(Y_out * gas.elemental_mass_fraction(element))
            
            # Relative error
            if mass_in > 1e-12:
                error = abs(mass_out - mass_in) / mass_in
            else:
                error = abs(mass_out - mass_in)
            
            balance_errors[element] = error
    
    return balance_errors

def compute_energy_residual(Q_chem, Q_applied, H_in, H_out):
    """
    Compute energy balance residual.
    
    Args:
        Q_chem: Chemical heat release (W)
        Q_applied: Applied heat (W) 
        H_in: Inlet enthalpy flow (W)
        H_out: Outlet enthalpy flow (W)
    
    Returns:
        Energy residual (dimensionless)
    """
    energy_in = Q_chem + Q_applied + H_in
    energy_out = H_out
    energy_diff = abs(energy_in - energy_out)
    
    # Normalize by total energy magnitude
    total_energy = abs(Q_chem) + abs(Q_applied) + abs(H_in) + abs(H_out)
    
    if total_energy > 1e-9:
        residual = energy_diff / total_energy
    else:
        residual = energy_diff
    
    return residual

def diagnose_case_adiabatic(case_name, mechanism_path, config_path):
    """
    Diagnose adiabatic bounds for a single case.
    
    Args:
        case_name: Case name
        mechanism_path: Path to mechanism file
        config_path: Path to config file
    
    Returns:
        dict: Diagnostic results
    """
    print(f"Diagnosing adiabatic bounds for {case_name}...")
    
    # Load mechanism
    thermo = ThermoManager(mechanism_path)
    gas = thermo.gas
    
    # Load case configuration
    model = HPPOXModel(mechanism_path)
    case_constants = CaseConstants()
    
    # Get case configuration
    case_config = case_constants.get_case_config(case_name)
    
    # Get inlet configuration (contains both case1 and case4)
    inlet_config = model.load_inlet_streams(config_path)
    
    # Create inlet state to get mixed composition
    inlet_state, m_dot_total = model.create_inlet_state(inlet_config, case_name)
    
    # Get inlet conditions
    T_in = inlet_state.temperature
    P = inlet_state.pressure
    Y_in = inlet_state.mass_fractions
    
    print(f"  Inlet: T={T_in:.1f}K, P={P/1000:.1f}kPa")
    print(f"  Composition: CH4={Y_in[thermo.species_names.index('CH4')]:.3f}, "
          f"O2={Y_in[thermo.species_names.index('O2')]:.3f}")
    
    # Compute adiabatic flame temperature
    T_ad, H_in, H_out, H_error = compute_adiabatic_flame_temperature(gas, P, T_in, Y_in)
    
    print(f"  Adiabatic flame: T_ad={T_ad:.1f}K")
    print(f"  Enthalpy error: {H_error:.2e}")
    
    # Compute elemental balance
    # For adiabatic case, outlet should be equilibrium composition
    gas.TPY = T_ad, P, Y_in
    gas.equilibrate('HP')
    Y_out_equil = gas.Y
    
    elemental_errors = compute_elemental_balance(gas, Y_in, Y_out_equil)
    
    print(f"  Elemental balance errors:")
    for element, error in elemental_errors.items():
        print(f"    {element}: {error:.2e}")
    
    # Store results
    results = {
        'case_name': case_name,
        'inlet_conditions': {
            'temperature': T_in,
            'pressure': P,
            'composition': dict(zip(thermo.species_names, Y_in))
        },
        'adiabatic_flame': {
            'temperature': T_ad,
            'enthalpy_in': H_in,
            'enthalpy_out': H_out,
            'enthalpy_error': H_error
        },
        'elemental_balance': elemental_errors,
        'physical_bounds': {
            'max_allowed_temp': T_ad + 30.0,
            'adiabatic_temp': T_ad
        }
    }
    
    return results

def main():
    """Main diagnostic function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Diagnose adiabatic bounds for HP-POX cases')
    parser.add_argument('--mechanism', required=True, help='Mechanism file path')
    parser.add_argument('--config', required=True, help='Config file path')
    parser.add_argument('--cases', nargs='+', required=True, help='Case names to diagnose')
    parser.add_argument('--output-dir', default='diagnostics', help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    results = {}
    
    for case_name in args.cases:
        try:
            case_results = diagnose_case_adiabatic(case_name, args.mechanism, args.config)
            results[case_name] = case_results
            
            # Save individual case results
            case_file = output_dir / f"T_ad_{case_name}.json"
            with open(case_file, 'w') as f:
                json.dump(case_results, f, indent=2)
            
            print(f"  ✓ Saved {case_file}")
            
        except Exception as e:
            print(f"  ✗ Error diagnosing {case_name}: {e}")
            results[case_name] = {'error': str(e)}
    
    # Save summary
    summary_file = output_dir / "adiabatic_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✓ Adiabatic bounds diagnostic complete")
    print(f"  Results saved to: {output_dir}")
    print(f"  Summary: {summary_file}")

if __name__ == "__main__":
    main()
