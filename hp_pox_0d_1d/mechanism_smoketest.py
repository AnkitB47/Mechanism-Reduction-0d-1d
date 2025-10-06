#!/usr/bin/env python3
"""
Mechanism Smoke Test for HP-POX 0D-1D Model

This script tests all available chemical mechanisms to ensure they load correctly
in Cantera and can perform basic equilibrium calculations under HP-POX conditions.

Tested mechanisms:
- SAN_DIEGO_MECH/c0_mech.yaml (read-only reference)
- SAN_DIEGO_MECH/c4_mech.yaml (read-only reference)  
- ARAMCOMECH3.0/aramco3.yaml (converted from CHEMKIN)
- USC_MECH_II/uscmech2.yaml (converted from CHEMKIN)
"""

import cantera as ct
import numpy as np
import sys
from pathlib import Path


def test_mechanism(mech_path: str, mech_name: str) -> dict:
    """
    Test a single mechanism for basic functionality.
    
    Args:
        mech_path: Path to the mechanism YAML file
        mech_name: Human-readable name for the mechanism
        
    Returns:
        Dictionary with test results
    """
    print(f"\n{'='*60}")
    print(f"Testing: {mech_name}")
    print(f"Path: {mech_path}")
    print(f"{'='*60}")
    
    results = {
        'name': mech_name,
        'path': mech_path,
        'load_success': False,
        'species_count': 0,
        'reaction_count': 0,
        'equilibrium_success': False,
        'equilibrium_composition': {},
        'warnings': []
    }
    
    try:
        # Test 1: Load mechanism
        print("1. Loading mechanism...")
        gas = ct.Solution(mech_path)
        results['load_success'] = True
        results['species_count'] = gas.n_species
        results['reaction_count'] = gas.n_reactions
        print(f"   ✓ Loaded successfully")
        print(f"   ✓ Species: {gas.n_species}")
        print(f"   ✓ Reactions: {gas.n_reactions}")
        
        # Test 2: Set HP-POX conditions and equilibrate
        print("2. Testing equilibrium calculation...")
        T = 1200.0  # K
        P = 50.0 * ct.one_atm  # 50 bar
        
        # Check what species are available and create appropriate composition
        available_species = gas.species_names
        
        if 'CH4' in available_species and 'O2' in available_species:
            X = 'CH4:1, O2:0.5, H2O:0.5, N2:0.01'  # HP-POX composition
        elif 'H2' in available_species and 'CO' in available_species:
            X = 'H2:1, CO:0.5, H2O:0.5, N2:0.01'
            print("   Note: Using H2/CO mixture (CH4 not available in this mechanism)")
        elif 'H2' in available_species and 'O2' in available_species:
            X = 'H2:1, O2:0.5, H2O:0.5, N2:0.01'
            print("   Note: Using H2/O2 mixture (CH4/CO not available in this mechanism)")
        else:
            # Use only available species
            species_list = []
            if 'H2' in available_species:
                species_list.append('H2:1')
            if 'O2' in available_species:
                species_list.append('O2:0.5')
            if 'H2O' in available_species:
                species_list.append('H2O:0.5')
            if 'N2' in available_species:
                species_list.append('N2:0.01')
            
            if not species_list:
                raise ValueError("No suitable species found for equilibrium test")
            
            X = ', '.join(species_list)
            print(f"   Note: Using limited species composition: {X}")
        
        gas.TPX = T, P, X
        gas.equilibrate('TP')
        
        results['equilibrium_success'] = True
        
        # Extract key species for HP-POX
        key_species = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'N2', 'O2']
        for species in key_species:
            if species in gas.species_names:
                results['equilibrium_composition'][species] = gas[species].X[0]
            else:
                results['equilibrium_composition'][species] = 0.0
        
        # Calculate H2/CO ratio
        h2_co_ratio = 0.0
        if results['equilibrium_composition']['H2'] > 0 and results['equilibrium_composition']['CO'] > 0:
            h2_co_ratio = results['equilibrium_composition']['H2'] / results['equilibrium_composition']['CO']
        
        print(f"   ✓ Equilibrium calculation successful")
        print(f"   ✓ Temperature: {gas.T:.1f} K")
        print(f"   ✓ Pressure: {gas.P/ct.one_atm:.1f} bar")
        print(f"   ✓ H2/CO ratio: {h2_co_ratio:.3f}")
        
        # Print key species composition
        print("   Key species composition (mole fractions):")
        for species in key_species:
            if species in results['equilibrium_composition']:
                print(f"     {species:4s}: {results['equilibrium_composition'][species]:.6f}")
        
        # Test 3: Check for numerical stability
        print("3. Checking numerical stability...")
        
        # Check for NaN or infinite values
        has_nan = False
        has_inf = False
        
        for i, name in enumerate(gas.species_names):
            if np.isnan(gas.X[i]) or np.isinf(gas.X[i]):
                if np.isnan(gas.X[i]):
                    has_nan = True
                if np.isinf(gas.X[i]):
                    has_inf = True
                results['warnings'].append(f"Species {name} has invalid mole fraction: {gas.X[i]}")
        
        if not has_nan and not has_inf:
            print("   ✓ No NaN or infinite values detected")
        else:
            if has_nan:
                print("   ⚠ NaN values detected in species composition")
            if has_inf:
                print("   ⚠ Infinite values detected in species composition")
        
        # Check temperature and pressure
        if np.isnan(gas.T) or np.isinf(gas.T) or gas.T <= 0:
            results['warnings'].append(f"Invalid temperature: {gas.T}")
            print(f"   ⚠ Invalid temperature: {gas.T}")
        else:
            print(f"   ✓ Temperature is valid: {gas.T:.1f} K")
            
        if np.isnan(gas.P) or np.isinf(gas.P) or gas.P <= 0:
            results['warnings'].append(f"Invalid pressure: {gas.P}")
            print(f"   ⚠ Invalid pressure: {gas.P}")
        else:
            print(f"   ✓ Pressure is valid: {gas.P/ct.one_atm:.1f} bar")
        
        print(f"\n✓ {mech_name} - ALL TESTS PASSED")
        
    except Exception as e:
        print(f"   ✗ Error: {e}")
        results['warnings'].append(f"Load/equilibrium error: {e}")
        print(f"\n✗ {mech_name} - TESTS FAILED")
    
    return results


def main():
    """Main function to run all mechanism tests."""
    print("HP-POX Mechanism Smoke Test")
    print("=" * 60)
    print("Testing all available chemical mechanisms for basic functionality")
    print("under HP-POX conditions (T=1200K, P=50bar, CH4+O2+H2O+N2)")
    
    # Define mechanisms to test
    mechanisms = [
        {
            'path': 'SAN_DIEGO_MECH/c0_mech.yaml',
            'name': 'San Diego Mech (Case 0)'
        },
        {
            'path': 'SAN_DIEGO_MECH/c4_mech.yaml', 
            'name': 'San Diego Mech (Case 4)'
        },
        {
            'path': 'ARAMCOMECH3.0/aramco3.yaml',
            'name': 'AramcoMech 3.0'
        },
        {
            'path': 'USC_MECH_II/uscmech2.yaml',
            'name': 'USC Mech II'
        }
    ]
    
    # Check if all mechanism files exist
    missing_files = []
    for mech in mechanisms:
        if not Path(mech['path']).exists():
            missing_files.append(mech['path'])
    
    if missing_files:
        print(f"\nError: Missing mechanism files:")
        for file in missing_files:
            print(f"  - {file}")
        print("\nPlease ensure all mechanism files are present before running tests.")
        return 1
    
    # Run tests
    all_results = []
    passed_count = 0
    
    for mech in mechanisms:
        result = test_mechanism(mech['path'], mech['name'])
        all_results.append(result)
        
        if result['load_success'] and result['equilibrium_success'] and len(result['warnings']) == 0:
            passed_count += 1
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total mechanisms tested: {len(mechanisms)}")
    print(f"Passed all tests: {passed_count}")
    print(f"Failed tests: {len(mechanisms) - passed_count}")
    
    print(f"\nDetailed Results:")
    for result in all_results:
        status = "PASS" if (result['load_success'] and result['equilibrium_success'] and len(result['warnings']) == 0) else "FAIL"
        print(f"  {result['name']:20s}: {status}")
        if result['load_success']:
            print(f"    Species: {result['species_count']:4d}, Reactions: {result['reaction_count']:4d}")
        if result['warnings']:
            for warning in result['warnings']:
                print(f"    Warning: {warning}")
    
    # Check if all tests passed
    if passed_count == len(mechanisms):
        print(f"\n✓ ALL MECHANISMS PASSED - Ready for HP-POX simulations!")
        return 0
    else:
        print(f"\n⚠ Some mechanisms failed - Check warnings above")
        return 1


if __name__ == "__main__":
    sys.exit(main())
