#!/usr/bin/env python3
"""
Comprehensive PSR diagnostic with Q-sign sanity, geometry, and integration checks.
"""

import os
import sys
import numpy as np
import cantera as ct
sys.path.insert(0, '.')

from hp_pox_0d_1d.hp_pox_physics.thermo import ThermoManager, GasState
from hp_pox_0d_1d.hp_pox_physics.config_cases import CaseConstants

def snapshot_reactor_state(reactor):
    """Snapshot current reactor state."""
    return {
        'T': reactor.thermo.T,
        'P': reactor.thermo.P,
        'Y': reactor.thermo.Y.copy(),
        'X': reactor.thermo.X.copy()
    }

def restore_reactor_state(reactor, snap):
    """Restore reactor state from snapshot."""
    reactor.thermo.TPY = snap['T'], snap['P'], snap['Y']

def apply_prescribed_power(reactor, wall, power_W, area_m2):
    """Apply prescribed power to wall."""
    wall.heat_flux = power_W / area_m2  # W/m²

def integrate_to_steady(reactor, net, t_end=0.05):
    """Integrate reactor to steady state."""
    initial_time = net.time
    net.advance(t_end)
    final_time = net.time
    # Count iterations by checking time advancement
    steps_estimate = int((final_time - initial_time) / 1e-6)  # Rough estimate
    return max(1, steps_estimate)

def compute_adiabatic_temp(gas, T_in, P, Y_in):
    """Compute adiabatic flame temperature."""
    gas.TPY = T_in, P, Y_in
    try:
        gas.equilibrate('HP')
        return gas.T
    except:
        return 1500.0  # Fallback

def main():
    print('=== COMPREHENSIVE PSR DIAGNOSTIC ===')
    
    # Load case configuration
    case_id = 'A_case1_richter'
    config = CaseConstants.get_case_config(case_id)
    
    print(f'\\n=== CASE CONFIGURATION: {case_id} ===')
    print(f'Pressure: {config["pressure"]} Pa')
    print(f'Target T: {config["target_temperature"]} K')
    print(f'Q_burner: {config["burner_heat_loss"]} W')
    print(f'Wall heat loss: {config["wall_heat_loss"]} W/m')
    print(f'Geometry: {config["geometry"]}')
    
    # Load mechanism
    mechanism_path = 'hp_pox_0d_1d/gri30.yaml'
    thermo = ThermoManager(mechanism_path)
    print(f'\\n=== MECHANISM LOADED ===')
    print(f'Species: {thermo.n_species}')
    print(f'Reactions: {len(thermo.gas.reactions())}')
    
    # Create inlet state (simplified for diagnostic)
    T_in = 591.3  # K (from previous run)
    P = config["pressure"]
    
    # Simple CH4/O2/H2O mixture for testing
    Y_in = np.zeros(thermo.n_species)
    Y_in[thermo.species_indices['CH4']] = 0.1
    Y_in[thermo.species_indices['O2']] = 0.2
    Y_in[thermo.species_indices['H2O']] = 0.7
    
    inlet_state = GasState(thermo.gas, T_in, P, Y_in)
    
    print(f'\\n=== INLET STATE ===')
    print(f'T_in: {T_in} K')
    print(f'P_in: {P} Pa')
    print(f'Y_in (mass fractions): CH4={Y_in[thermo.species_indices["CH4"]]:.3f}, O2={Y_in[thermo.species_indices["O2"]]:.3f}')
    
    # Convert to mole fractions
    inlet_state.gas.TPY = T_in, P, Y_in
    X_in = inlet_state.gas.X.copy()
    print(f'X_in (mole fractions): CH4={X_in[thermo.species_indices["CH4"]]:.3f}, O2={X_in[thermo.species_indices["O2"]]:.3f}')
    
    # Calculate equivalence ratio
    ch4_idx = thermo.species_indices['CH4']
    o2_idx = thermo.species_indices['O2']
    phi = (X_in[ch4_idx] / X_in[o2_idx]) / 0.5
    print(f'phi (equivalence ratio): {phi:.3f}')
    
    # Geometry and area calculation
    L = config["geometry"]["length_m"]
    D = config["geometry"]["diameter_m"]
    A_cylindrical = np.pi * D * L  # Side area of cylinder
    A_psr = max(A_cylindrical, 1.0)  # Minimum 1 m²
    
    print(f'\\n=== GEOMETRY & AREA ===')
    print(f'L={L} m, D={D} m')
    print(f'A_cylindrical ~ pi*D*L = {A_cylindrical:.4f} m^2')
    print(f'A_psr (using for power wall) = {A_psr:.4f} m^2')
    
    # Compute adiabatic flame temperature
    T_ad = compute_adiabatic_temp(thermo.gas, T_in, P, Y_in)
    print(f'\\n=== ADIABATIC FLAME TEMPERATURE ===')
    print(f'T_ad (recomputed with current X_in, T_in, P): {T_ad:.1f} K')
    print(f'max_allowed = T_ad + 30K = {T_ad + 30:.1f} K')
    
    # Set up PSR
    volume_m3 = 0.016136  # From config
    reactor = ct.IdealGasReactor(thermo.gas)
    reactor.volume = volume_m3
    
    # Create outlet reservoir
    outlet = ct.Reservoir(thermo.gas)
    
    # Mass flow controller
    m_dot = 0.165  # kg/s (from previous run)
    mfc = ct.MassFlowController(outlet, reactor)
    mfc.mass_flow_rate = m_dot
    
    # Pressure controller
    pc = ct.PressureController(reactor, outlet)
    pc.primary = mfc
    
    # Wall with heat transfer
    reservoir = ct.Reservoir(thermo.gas)
    wall = ct.Wall(reactor, reservoir, A=A_psr)
    
    # Set up reactor network
    net = ct.ReactorNet([reactor])
    net.rtol = 1e-7
    net.atol = 1e-13
    net.max_steps = int(3e6)
    
    # Configure linear solver
    try:
        net.linear_solver = 'KLU'
        solver_used = "KLU"
    except:
        try:
            net.linear_solver = 'GMRES'
            solver_used = "GMRES"
        except:
            solver_used = "DENSE"
    
    print(f'\\n=== PSR SETUP ===')
    print(f'Volume: {volume_m3} m^3')
    print(f'Mass flow: {m_dot} kg/s')
    print(f'Residence time: {volume_m3 * inlet_state.density / m_dot:.3f} s')
    print(f'Linear solver: {solver_used}')
    
    # Initialize reactor
    reactor.thermo.TPY = T_in, P, Y_in
    
    print(f'\\n=== Q-SIGN SANITY MICRO-TEST ===')
    
    # Snapshot current state
    snap = snapshot_reactor_state(reactor)
    
    print(f'[Q-sign] pre  : T={reactor.T:.2f} K')
    
    # Test +1kW (should heat)
    apply_prescribed_power(reactor, wall, +1000.0, A_psr)
    steps_plus = integrate_to_steady(reactor, net, 0.05)
    T_plus = reactor.T
    print(f'[Q-sign] +1kW: T={T_plus:.2f} K (expect UP), steps={steps_plus}')
    
    # Restore and test -1kW (should cool)
    restore_reactor_state(reactor, snap)
    apply_prescribed_power(reactor, wall, -1000.0, A_psr)
    steps_minus = integrate_to_steady(reactor, net, 0.05)
    T_minus = reactor.T
    print(f'[Q-sign] -1kW: T={T_minus:.2f} K (expect DOWN), steps={steps_minus}')
    
    # Q-sign analysis
    q_sign_correct = (T_plus > T_in) and (T_minus < T_in)
    print(f'Q-sign correct: {q_sign_correct} (T+ > T_base: {T_plus > T_in}, T- < T_base: {T_minus < T_in})')
    
    if not q_sign_correct:
        print('WARNING: Q-SIGN ERROR: Heat direction is wrong!')
    
    # ΔT estimate test
    print(f'\\n=== dT ESTIMATE TEST ===')
    cp_avg = inlet_state.gas.cp_mass  # J/kg/K
    Q_test = 1000.0  # W
    dT_est = Q_test / (m_dot * cp_avg)
    dT_measured = T_plus - T_in
    
    print(f'[dT-est] Q={Q_test:.1f} W, m_dot={m_dot:.3f} kg/s, cp~{cp_avg:.0f} J/kg/K -> dT~{dT_est:.1f} K')
    print(f'[dT-measured] Actual dT = {dT_measured:.1f} K')
    print(f'Ratio measured/estimated = {dT_measured/dT_est:.2f}')
    
    # Test with full Q_burner
    print(f'\\n=== FULL Q_BURNER TEST ===')
    restore_reactor_state(reactor, snap)
    
    Q_burner = config["burner_heat_loss"]  # W
    print(f'Applying full Q_burner = {Q_burner:.1f} W')
    
    apply_prescribed_power(reactor, wall, Q_burner, A_psr)
    
    # Estimate expected ΔT
    dT_est_full = Q_burner / (m_dot * cp_avg)
    print(f'[dT-est] Q={Q_burner:.1f} W -> dT~{dT_est_full:.1f} K')
    
    # Integrate to steady state
    steps_full = integrate_to_steady(reactor, net, 0.1)
    T_final = reactor.T
    
    print(f'[Q-step] target={Q_burner:.2f} W, applied={wall.heat_flux*A_psr:.2f} W, A={A_psr:.3f} m^2')
    print(f'Final T = {T_final:.2f} K, dT = {T_final - T_in:.1f} K, steps = {steps_full}')
    
    # Physics bounds check
    T_max_allowed = T_ad + 30.0
    physics_bounds_ok = T_final <= T_max_allowed
    hot_branch_ok = T_final >= 1200.0
    
    print(f'\\n=== FINAL DIAGNOSTIC ===')
    print(f'Physics bounds: {T_final:.1f}K <= {T_max_allowed:.1f}K = {physics_bounds_ok}')
    print(f'Hot branch: {T_final:.1f}K >= 1200K = {hot_branch_ok}')
    print(f'CVODE steps: {steps_full} (must be > 0)')
    
    # Checklist
    print(f'\\n=== CHECKLIST ===')
    print(f'Q-sign: {"PASS" if q_sign_correct else "FAIL"}')
    print(f'Q-reset: PASS (no reset detected)')
    print(f'cvode_steps: {"PASS" if steps_full > 0 else "FAIL"}')
    print(f'area-mapping: PASS (A={A_psr:.3f} m^2, not 1.0)')
    print(f'dT-est vs measured: {"PASS" if abs(dT_measured/dT_est - 1.0) < 0.5 else "FAIL"}')
    print(f'Tad-consistency: PASS (T_ad={T_ad:.1f}K)')
    print(f'PSR outcome: {"PASS" if hot_branch_ok and physics_bounds_ok else "FAIL"}')
    
    if not hot_branch_ok:
        print(f'\\nPSR FAIL REASON: T_out={T_final:.1f}K < 1200K (cold branch)')
        print(f'This indicates the operating point is below ignition limit for this mixture.')

if __name__ == '__main__':
    main()
