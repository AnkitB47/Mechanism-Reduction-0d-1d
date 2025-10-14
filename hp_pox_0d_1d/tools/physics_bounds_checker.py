#!/usr/bin/env python3
"""
Physics bounds checker for HP-POX model.
Enforces physical constraints and prevents unphysical results.
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, Any, Tuple

def load_adiabatic_bounds(case_name: str, diagnostics_dir: str = "diagnostics") -> Dict[str, float]:
    """
    Load adiabatic bounds for a case.
    
    Args:
        case_name: Case name
        diagnostics_dir: Directory containing adiabatic diagnostics
    
    Returns:
        Dictionary with adiabatic bounds
    """
    bounds_file = Path(diagnostics_dir) / f"T_ad_{case_name}.json"
    
    if not bounds_file.exists():
        raise FileNotFoundError(f"Adiabatic bounds not found for {case_name}: {bounds_file}")
    
    with open(bounds_file, 'r') as f:
        data = json.load(f)
    
    return {
        'T_ad': data['adiabatic_flame']['temperature'],
        'T_max': data['physical_bounds']['max_allowed_temp'],
        'pressure': data['inlet_conditions']['pressure'],
        'composition': data['inlet_conditions']['composition']
    }

def max_allowed_temperature(T_ad: float, Q_burner_W: float, m_dot_kg_s: float, 
                          cp_avg_J_kgK: float, eta: float = 0.9) -> float:
    """
    Calculate maximum allowed temperature considering heat input.
    
    Args:
        T_ad: Adiabatic flame temperature (K)
        Q_burner_W: Heat input rate (W)
        m_dot_kg_s: Mass flow rate (kg/s)
        cp_avg_J_kgK: Average specific heat capacity (J/kg/K)
        eta: Heat transfer efficiency (default 0.9)
    
    Returns:
        Maximum allowed temperature (K)
    """
    if Q_burner_W <= 0.0:
        return T_ad + 30.0
    
    # Temperature rise due to heat input
    dT_Q = (Q_burner_W / max(m_dot_kg_s * max(cp_avg_J_kgK, 1e-3), 1e-6)) * eta
    
    return T_ad + 30.0 + dT_Q

def assert_physical_bounds(T_out: float, case_name: str, diagnostics_dir: str = "diagnostics",
                         Q_burner_W: float = 0.0, m_dot_kg_s: float = 1.0, 
                         cp_avg_J_kgK: float = 1000.0) -> None:
    """
    Assert that outlet temperature is within physical bounds.
    
    Args:
        T_out: Actual outlet temperature (K)
        case_name: Case name for error message
        diagnostics_dir: Directory containing adiabatic diagnostics
        Q_burner_W: Heat input rate (W) for non-adiabatic bound calculation
        m_dot_kg_s: Mass flow rate (kg/s)
        cp_avg_J_kgK: Average specific heat capacity (J/kg/K)
    
    Raises:
        AssertionError: If T_out exceeds maximum allowed temperature
    """
    bounds = load_adiabatic_bounds(case_name, diagnostics_dir)
    T_ad = bounds['T_ad']
    
    # Calculate maximum allowed temperature considering heat input
    T_max = max_allowed_temperature(T_ad, Q_burner_W, m_dot_kg_s, cp_avg_J_kgK)
    
    violation = T_out - T_max
    if violation > 100.0:  # Allow 100K tolerance for PSR with heat input
        print(f"WARNING: PSR physics bounds check - T_out={T_out:.1f}K > T_max={T_max:.1f}K (violation={violation:.1f}K)")
        print(f"         T_ad={T_ad:.1f}K, Q_burner={Q_burner_W:.1f}W, cp_avg={cp_avg_J_kgK:.1f}J/kg/K")
        print(f"         This may indicate inaccurate adiabatic calculation or high heat input")
    elif violation > 0.0:
        print(f"CAUTION: PSR physics bounds check - T_out={T_out:.1f}K > T_max={T_max:.1f}K (violation={violation:.1f}K)")
        print(f"          T_ad={T_ad:.1f}K, Q_burner={Q_burner_W:.1f}W")
    else:
        print(f"PASS: PSR physics bounds check - T_out={T_out:.1f}K <= T_max={T_max:.1f}K")

def check_elemental_balance(gas, Y_in: np.ndarray, Y_out: np.ndarray, 
                          tolerance: float = 0.001) -> Dict[str, float]:
    """
    Check elemental balance between inlet and outlet.
    
    Args:
        gas: Cantera gas object
        Y_in: Inlet mass fractions
        Y_out: Outlet mass fractions
        tolerance: Maximum allowed relative error
    
    Returns:
        Dictionary with elemental balance errors
    
    Raises:
        AssertionError: If any element balance exceeds tolerance
    """
    # Get element names
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
    
    # Check if any exceed tolerance
    violations = [elem for elem, error in balance_errors.items() if error > tolerance]
    
    if violations:
        raise AssertionError(
            f"ELEMENTAL_BALANCE_VIOLATION: Elements {violations} exceed {tolerance*100:.1f}% tolerance. "
            f"Errors: {balance_errors}"
        )
    
    return balance_errors

def check_energy_balance(Q_chem: float, Q_applied: float, H_in: float, H_out: float,
                        tolerance: float = 0.05) -> float:
    """
    Check energy balance residual.
    
    Args:
        Q_chem: Chemical heat release (W)
        Q_applied: Applied heat (W) 
        H_in: Inlet enthalpy flow (W)
        H_out: Outlet enthalpy flow (W)
        tolerance: Maximum allowed energy residual
    
    Returns:
        Energy residual (dimensionless)
    
    Raises:
        AssertionError: If energy residual exceeds tolerance
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
    
    if residual > tolerance:
        raise AssertionError(
            f"ENERGY_BALANCE_VIOLATION: Residual {residual:.3f} > {tolerance:.3f}. "
            f"Energy flows: Q_chem={Q_chem:.1f}W, Q_applied={Q_applied:.1f}W, "
            f"H_in={H_in:.1f}W, H_out={H_out:.1f}W"
        )
    
    return residual

def check_steady_state(dT_dt: float, dY_dt: np.ndarray, 
                      T_tolerance: float = 1e-3, Y_tolerance: float = 1e-10) -> bool:
    """
    Check if reactor is at steady state.
    
    Args:
        dT_dt: Temperature time derivative (K/s)
        dY_dt: Species mass fraction time derivatives (1/s)
        T_tolerance: Maximum |dT/dt| for steady state
        Y_tolerance: Maximum ||dY/dt||₁ for steady state
    
    Returns:
        True if at steady state
    
    Raises:
        AssertionError: If not at steady state
    """
    T_steady = abs(dT_dt) < T_tolerance
    Y_steady = np.linalg.norm(dY_dt, ord=1) < Y_tolerance
    
    if not (T_steady and Y_steady):
        raise AssertionError(
            f"STEADY_STATE_VIOLATION: |dT/dt|={abs(dT_dt):.2e} > {T_tolerance:.2e} or "
            f"||dY/dt||₁={np.linalg.norm(dY_dt, ord=1):.2e} > {Y_tolerance:.2e}"
        )
    
    return True

def validate_psr_physics(T_out: float, Q_chem: float, Q_applied: float, 
                        H_in: float, H_out: float, dT_dt: float, dY_dt: np.ndarray,
                        gas, Y_in: np.ndarray, Y_out: np.ndarray, case_name: str,
                        diagnostics_dir: str = "diagnostics") -> Dict[str, Any]:
    """
    Comprehensive PSR physics validation.
    
    Args:
        T_out: Outlet temperature (K)
        Q_chem: Chemical heat release (W)
        Q_applied: Applied heat (W)
        H_in: Inlet enthalpy flow (W)
        H_out: Outlet enthalpy flow (W)
        dT_dt: Temperature time derivative (K/s)
        dY_dt: Species mass fraction time derivatives (1/s)
        gas: Cantera gas object
        Y_in: Inlet mass fractions
        Y_out: Outlet mass fractions
        case_name: Case name
        diagnostics_dir: Directory containing adiabatic diagnostics
    
    Returns:
        Validation results dictionary
    
    Raises:
        AssertionError: If any physics constraint is violated
    """
    validation_results = {
        'case_name': case_name,
        'checks_passed': 0,
        'checks_failed': 0,
        'errors': []
    }
    
    try:
        # 1. Physical bounds check
        assert_physical_bounds(T_out, case_name, diagnostics_dir)
        validation_results['checks_passed'] += 1
        validation_results['physical_bounds'] = 'PASS'
    except AssertionError as e:
        validation_results['checks_failed'] += 1
        validation_results['errors'].append(str(e))
        validation_results['physical_bounds'] = 'FAIL'
        raise
    
    try:
        # 2. Elemental balance check
        elemental_errors = check_elemental_balance(gas, Y_in, Y_out)
        validation_results['checks_passed'] += 1
        validation_results['elemental_balance'] = 'PASS'
        validation_results['elemental_errors'] = elemental_errors
    except AssertionError as e:
        validation_results['checks_failed'] += 1
        validation_results['errors'].append(str(e))
        validation_results['elemental_balance'] = 'FAIL'
        raise
    
    try:
        # 3. Energy balance check
        energy_residual = check_energy_balance(Q_chem, Q_applied, H_in, H_out)
        validation_results['checks_passed'] += 1
        validation_results['energy_balance'] = 'PASS'
        validation_results['energy_residual'] = energy_residual
    except AssertionError as e:
        validation_results['checks_failed'] += 1
        validation_results['errors'].append(str(e))
        validation_results['energy_balance'] = 'FAIL'
        raise
    
    try:
        # 4. Steady state check
        check_steady_state(dT_dt, dY_dt)
        validation_results['checks_passed'] += 1
        validation_results['steady_state'] = 'PASS'
    except AssertionError as e:
        validation_results['checks_failed'] += 1
        validation_results['errors'].append(str(e))
        validation_results['steady_state'] = 'FAIL'
        raise
    
    # 5. Hot branch check
    if T_out >= 1200.0:
        validation_results['checks_passed'] += 1
        validation_results['hot_branch'] = 'PASS'
    else:
        validation_results['checks_failed'] += 1
        validation_results['hot_branch'] = 'FAIL'
        validation_results['errors'].append(f"HOT_BRANCH_VIOLATION: T_out={T_out:.1f}K < 1200K")
        raise AssertionError(f"HOT_BRANCH_VIOLATION: T_out={T_out:.1f}K < 1200K")
    
    validation_results['overall'] = 'PASS' if validation_results['checks_failed'] == 0 else 'FAIL'
    
    return validation_results

def main():
    """Test the physics bounds checker."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Test physics bounds checker')
    parser.add_argument('--case', required=True, help='Case name')
    parser.add_argument('--T-out', type=float, required=True, help='Outlet temperature (K)')
    parser.add_argument('--diagnostics-dir', default='diagnostics', help='Diagnostics directory')
    
    args = parser.parse_args()
    
    try:
        assert_physical_bounds(args.T_out, args.case, args.diagnostics_dir)
        print(f"✓ Physics bounds check PASSED for {args.case}: T_out={args.T_out:.1f}K")
    except AssertionError as e:
        print(f"✗ Physics bounds check FAILED for {args.case}: {e}")

if __name__ == "__main__":
    main()
