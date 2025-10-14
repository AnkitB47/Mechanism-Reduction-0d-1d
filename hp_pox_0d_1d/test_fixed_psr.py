#!/usr/bin/env python3
"""
Test the fixed PSR implementation with correct heat transfer.
"""

import os
import sys
sys.path.insert(0, '.')

from hp_pox_0d_1d.hp_pox_physics.thermo import ThermoManager, GasState
from hp_pox_0d_1d.hp_pox_physics.config_cases import CaseConstants
import numpy as np

def test_fixed_psr():
    print("=== TESTING FIXED PSR IMPLEMENTATION ===")
    
    # Load case configuration
    case_id = 'A_case1_richter'
    config = CaseConstants.get_case_config(case_id)
    
    print(f"Case: {case_id}")
    print(f"Q_burner: {config['burner_heat_loss']} W")
    print(f"Geometry: {config['geometry']}")
    
    # Load mechanism
    mechanism_path = 'hp_pox_0d_1d/gri30.yaml'
    thermo = ThermoManager(mechanism_path)
    
    # Create inlet state (simplified)
    T_in = 591.3  # K
    P = config["pressure"]
    
    # Simple CH4/O2/H2O mixture
    Y_in = np.zeros(thermo.n_species)
    Y_in[thermo.species_indices['CH4']] = 0.1
    Y_in[thermo.species_indices['O2']] = 0.2
    Y_in[thermo.species_indices['H2O']] = 0.7
    
    inlet_state = GasState(thermo.gas, T_in, P, Y_in)
    
    print(f"\\nInlet state:")
    print(f"  T_in: {T_in} K")
    print(f"  P_in: {P} Pa")
    print(f"  Y_CH4: {Y_in[thermo.species_indices['CH4']]:.3f}")
    
    # Test the PSR with fixed heat transfer
    from hp_pox_0d_1d.hp_pox_physics.psr import PSR
    
    # Create PSR instance
    volume_m3 = 0.016136
    psr = PSR(volume_m3, thermo)
    psr.case_id = case_id  # For diagnostics
    
    print(f"\\n=== PSR TEST WITH FIXED HEAT TRANSFER ===")
    
    # Test parameters
    m_dot = 0.165  # kg/s
    q_burner_kw = config["burner_heat_loss"] / 1000.0  # Convert W to kW
    
    print(f"Test parameters:")
    print(f"  Volume: {volume_m3} m^3")
    print(f"  Mass flow: {m_dot} kg/s")
    print(f"  Q_burner: {q_burner_kw} kW ({q_burner_kw*1000:.0f} W)")
    
    # Run PSR with the fixed implementation
    try:
        result = psr.solve_psr_continuation(
            inlet_state=inlet_state,
            m_dot=m_dot,
            q_burner_kw=q_burner_kw,
            P=P
        )
        
        # PSR returns (outlet_state, converged, diagnostics)
        outlet_state, converged, diagnostics = result
        
        print(f"\\n=== PSR RESULT ===")
        print(f"Converged: {converged}")
        print(f"T_out: {outlet_state.temperature:.1f} K")
        print(f"Branch: {diagnostics.get('branch_status', 'UNKNOWN')}")
        print(f"Q_chem: {diagnostics.get('q_chem_final', 0.0):.1f} W")
        print(f"Q_burner: {diagnostics.get('final_q_burner', 0.0):.1f} W")
        
        if converged:
            print("PASS: PSR CONVERGED SUCCESSFULLY!")
        else:
            print("FAIL: PSR FAILED TO CONVERGE")
            
    except Exception as e:
        print(f"ERROR: PSR ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    test_fixed_psr()
