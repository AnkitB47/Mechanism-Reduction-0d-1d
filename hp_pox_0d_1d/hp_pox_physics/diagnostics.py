# hp_pox_physics/diagnostics.py
"""
Diagnostics and validation for HP-POX simulations.
Checks elemental balance, energy residual, and KPI matching.
"""

import cantera as ct
import numpy as np
from typing import Dict, Any, Tuple


def check_elemental_balance(outlet_state, inlet_state, m_dot_out: float, m_dot_in: float) -> Dict[str, float]:
    """Check elemental balance between inlet and outlet."""
    # Elemental conservation check (C, H, O, N atoms)
    elements = ['C', 'H', 'O', 'N']

    inlet_atoms = {}
    outlet_atoms = {}

    for elem in elements:
        # Handle both GasState and Solution objects
        if hasattr(inlet_state, 'gas'):
            inlet_atoms[elem] = inlet_state.gas.elemental_mass_fraction(elem) * m_dot_in
        else:
            inlet_atoms[elem] = inlet_state.elemental_mass_fraction(elem) * m_dot_in

        if hasattr(outlet_state, 'gas'):
            outlet_atoms[elem] = outlet_state.gas.elemental_mass_fraction(elem) * m_dot_out
        else:
            outlet_atoms[elem] = outlet_state.elemental_mass_fraction(elem) * m_dot_out

    # Calculate relative errors
    errors = {}
    for elem in elements:
        if inlet_atoms[elem] > 1e-10:  # Avoid division by zero
            errors[elem] = abs(inlet_atoms[elem] - outlet_atoms[elem]) / inlet_atoms[elem]
        else:
            errors[elem] = abs(inlet_atoms[elem] - outlet_atoms[elem])

    max_error = max(errors.values()) if errors else 0.0

    return {
        'elemental_errors': errors,
        'max_elemental_error': max_error,
        'elemental_balance_pass': max_error <= 0.001  # 0.1% relative error
    }


def check_energy_residual(
    inlet_state,
    outlet_state,
    m_dot_in: float,
    m_dot_out: float,
    q_burner: float,
    q_wall: float = 0.0
) -> Dict[str, float]:
    """Check energy balance residual."""
    # Enthalpy balance: H_in + Q_burner = H_out + Q_wall
    # Note: Q_burner positive means heat input to reactor
    # Q_wall positive means heat loss from reactor

    # Handle both GasState and Solution objects
    if hasattr(inlet_state, 'gas'):
        h_in = inlet_state.gas.enthalpy_mass * m_dot_in  # W
        h_out = outlet_state.gas.enthalpy_mass * m_dot_out  # W
    else:
        h_in = inlet_state.enthalpy_mass * m_dot_in  # W
        h_out = outlet_state.enthalpy_mass * m_dot_out  # W

    # Energy residual = |H_in + Q_burner - H_out - Q_wall| / max(|H_in|, |H_out|, |Q_burner|, |Q_wall|, 1e-9)
    residual = abs(h_in + q_burner - h_out - q_wall)
    denominator = max(abs(h_in), abs(h_out), abs(q_burner), abs(q_wall), 1e-9)

    return {
        'energy_residual': residual / denominator,
        'h_in': h_in,
        'h_out': h_out,
        'q_burner': q_burner,
        'q_wall': q_wall,
        'energy_balance_pass': (residual / denominator) <= 0.05  # 5% error
    }


def check_kpis(outlet_state, case_config) -> Dict[str, Any]:
    """Check KPIs against targets."""
    # Get dry outlet composition (exclude H2O)
    dry_species = {}
    wet_species = {}
    total_dry = 0.0

    # Handle both GasState and Solution objects
    if hasattr(outlet_state, 'gas'):
        gas_obj = outlet_state.gas
    else:
        gas_obj = outlet_state

    # Helper for safe species lookup
    def _idx(g, name):
        return g.species_index(name) if name in g.species_names else None

    # One-time lookup
    i_H2  = _idx(gas_obj, "H2")
    i_CO  = _idx(gas_obj, "CO")
    i_CO2 = _idx(gas_obj, "CO2")
    i_CH4 = _idx(gas_obj, "CH4")
    i_H2O = _idx(gas_obj, "H2O")
    i_O2  = _idx(gas_obj, "O2")
    i_N2  = _idx(gas_obj, "N2")

    # Safe reads (skip None)
    y_H2 = gas_obj.Y[i_H2] if i_H2 is not None else 0.0
    x_H2 = gas_obj.X[i_H2] if i_H2 is not None else 0.0
    y_CO = gas_obj.Y[i_CO] if i_CO is not None else 0.0
    x_CO = gas_obj.X[i_CO] if i_CO is not None else 0.0
    y_CO2 = gas_obj.Y[i_CO2] if i_CO2 is not None else 0.0
    x_CO2 = gas_obj.X[i_CO2] if i_CO2 is not None else 0.0
    y_CH4 = gas_obj.Y[i_CH4] if i_CH4 is not None else 0.0
    x_CH4 = gas_obj.X[i_CH4] if i_CH4 is not None else 0.0
    y_H2O = gas_obj.Y[i_H2O] if i_H2O is not None else 0.0
    x_H2O = gas_obj.X[i_H2O] if i_H2O is not None else 0.0
    y_O2 = gas_obj.Y[i_O2] if i_O2 is not None else 0.0
    x_O2 = gas_obj.X[i_O2] if i_O2 is not None else 0.0
    y_N2 = gas_obj.Y[i_N2] if i_N2 is not None else 0.0
    x_N2 = gas_obj.X[i_N2] if i_N2 is not None else 0.0

    # Calculate dry basis (exclude H2O)
    total_dry = y_H2 + y_CO + y_CO2 + y_CH4 + y_N2
    if total_dry > 0:
        dry_H2 = y_H2 / total_dry
        dry_CO = y_CO / total_dry
        dry_CO2 = y_CO2 / total_dry
        dry_CH4 = y_CH4 / total_dry
        dry_N2 = y_N2 / total_dry
    else:
        dry_H2 = dry_CO = dry_CO2 = dry_CH4 = dry_N2 = 0.0

    # Actual values (dry basis unless noted)
    actual_T = gas_obj.T - 273.15  # Convert to °C
    actual_H2 = dry_H2 * 100
    actual_N2 = dry_N2 * 100
    actual_CO = dry_CO * 100
    actual_CO2 = dry_CO2 * 100
    actual_CH4 = dry_CH4 * 100
    actual_H2O_wet = y_H2O * 100

    # Target values
    target_T = case_config.target_T_K - 273.15
    target_H2 = case_config.target_H2_pct
    target_N2 = case_config.target_N2_pct
    target_CO = case_config.target_CO_pct
    target_CO2 = case_config.target_CO2_pct
    target_CH4 = case_config.target_CH4_pct
    target_H2O_wet = case_config.target_H2O_wet_pct

    # Deviations
    deviations = {
        'T_C': abs(actual_T - target_T),
        'H2_pct': abs(actual_H2 - target_H2),
        'N2_pct': abs(actual_N2 - target_N2),
        'CO_pct': abs(actual_CO - target_CO),
        'CO2_pct': abs(actual_CO2 - target_CO2),
        'CH4_pct': abs(actual_CH4 - target_CH4),
        'H2O_wet_pct': abs(actual_H2O_wet - target_H2O_wet)
    }

    # H2/CO ratio
    h2_co_ratio = dry_H2 / dry_CO if dry_CO > 0 else 0.0

    # Check if all deviations are within tolerance
    tolerances = {
        'T_C': 20.0,  # ±20 K
        'H2_pct': 1.0,  # ±1.0%
        'N2_pct': 1.0,  # ±1.0%
        'CO_pct': 1.0,  # ±1.0%
        'CO2_pct': 1.0,  # ±1.0%
        'CH4_pct': 1.0,  # ±1.0%
        'H2O_wet_pct': 1.0  # ±1.0%
    }

    kpi_pass = all(dev <= tol for (dev, tol) in zip(deviations.values(), tolerances.values()))

    return {
        'kpi_pass': kpi_pass,
        'deviations': deviations,
        'actual_values': {
            'T_C': actual_T,
            'H2_pct': actual_H2,
            'N2_pct': actual_N2,
            'CO_pct': actual_CO,
            'CO2_pct': actual_CO2,
            'CH4_pct': actual_CH4,
            'H2O_wet_pct': actual_H2O_wet,
            'H2_CO_ratio': h2_co_ratio
        },
        'target_values': {
            'T_C': target_T,
            'H2_pct': target_H2,
            'N2_pct': target_N2,
            'CO_pct': target_CO,
            'CO2_pct': target_CO2,
            'CH4_pct': target_CH4,
            'H2O_wet_pct': target_H2O_wet
        }
    }


def run_comprehensive_diagnostics(
    inlet_state: ct.Solution,
    outlet_state: ct.Solution,
    m_dot_in: float,
    m_dot_out: float,
    q_burner: float,
    q_wall: float,
    case_config
) -> Dict[str, Any]:
    """Run all diagnostics and return comprehensive results."""
    # Elemental balance
    elem_balance = check_elemental_balance(outlet_state, inlet_state, m_dot_out, m_dot_in)

    # Energy balance
    energy_balance = check_energy_residual(inlet_state, outlet_state, m_dot_in, m_dot_out, q_burner, q_wall)

    # KPI check
    kpi_check = check_kpis(outlet_state, case_config)

    # Overall pass/fail
    overall_pass = elem_balance['elemental_balance_pass'] and energy_balance['energy_balance_pass'] and kpi_check['kpi_pass']

    return {
        'overall_pass': overall_pass,
        'elemental_balance': elem_balance,
        'energy_balance': energy_balance,
        'kpi_check': kpi_check
    }
