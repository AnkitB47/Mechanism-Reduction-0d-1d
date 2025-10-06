"""
Stream mixing with enthalpy-consistent temperature calculation.
Implements proper enthalpy-weighted mixing for multiple inlet streams.
"""

import numpy as np
from typing import List, Dict, Any, Tuple
from scipy.optimize import brentq
from .thermo import ThermoManager


def mix_streams(streams: List[Dict[str, Any]], thermo: ThermoManager) -> Tuple[float, np.ndarray, float]:
    """Mix multiple streams with enthalpy-consistent temperature calculation.
    
    Args:
        streams: List of stream dictionaries with keys:
            - m_dot: mass flow rate (kg/s)
            - T: temperature (K)
            - Y: mass fractions dict {species: fraction}
            - P: pressure (Pa)
        thermo: Thermodynamics manager
        
    Returns:
        Tuple of (T_mix, Y_mix, m_dot_total)
    """
    if not streams:
        raise ValueError("No streams provided for mixing")
    
    # Convert mass fractions to arrays
    n_species = thermo.n_species
    m_dot_total = 0.0
    Y_mix = np.zeros(n_species)
    h_mix = 0.0
    
    # Step 1: Calculate mass-weighted composition and enthalpy
    for stream in streams:
        m_dot = stream['m_dot']
        T = stream['T']
        P = stream['P']
        Y_dict = stream['Y']
        
        # Convert Y_dict to array
        Y = np.zeros(n_species)
        for species, fraction in Y_dict.items():
            if species in thermo.species_indices:
                idx = thermo.species_indices[species]
                Y[idx] = fraction
        
        # Normalize mass fractions
        Y = Y / np.sum(Y)
        
        # Set gas state and get enthalpy
        thermo.gas.TPY = T, P, Y
        h = thermo.gas.enthalpy_mass  # J/kg
        
        # Accumulate
        m_dot_total += m_dot
        Y_mix += m_dot * Y
        h_mix += m_dot * h
    
    # Step 2: Calculate mass-weighted average composition
    Y_mix = Y_mix / m_dot_total
    Y_mix = Y_mix / np.sum(Y_mix)  # Renormalize
    
    # Step 3: Calculate mass-weighted average enthalpy
    h_mix = h_mix / m_dot_total  # J/kg
    
    # Step 4: Find temperature that gives this enthalpy
    T_mix = _find_temperature_for_enthalpy(thermo, h_mix, Y_mix, streams[0]['P'])
    
    return T_mix, Y_mix, m_dot_total


def _find_temperature_for_enthalpy(thermo: ThermoManager, h_target: float, 
                                 Y: np.ndarray, P: float) -> float:
    """Find temperature that gives target enthalpy using secant method.
    
    Args:
        thermo: Thermodynamics manager
        h_target: Target specific enthalpy (J/kg)
        Y: Mass fractions
        P: Pressure (Pa)
        
    Returns:
        Temperature (K)
    """
    def enthalpy_residual(T):
        thermo.gas.TPY = T, P, Y
        h = thermo.gas.enthalpy_mass
        return h - h_target
    
    # Bracket the solution
    T_min = 300.0  # K
    T_max = 3000.0  # K
    
    # Check if solution is bracketed
    try:
        h_min = enthalpy_residual(T_min)
        h_max = enthalpy_residual(T_max)
        
        if h_min * h_max > 0:
            # Not bracketed, use linear interpolation as initial guess
            T_guess = T_min + (T_max - T_min) * (0 - h_min) / (h_max - h_min)
            T_guess = max(T_min, min(T_max, T_guess))
        else:
            T_guess = (T_min + T_max) / 2
        
        # Use Brent's method for robust root finding
        T_solution = brentq(enthalpy_residual, T_min, T_max, xtol=1e-6)
        return T_solution
        
    except Exception as e:
        print(f"Warning: Root finding failed ({e}), using linear interpolation")
        # Fallback to linear interpolation
        T_guess = T_min + (T_max - T_min) * (0 - h_min) / (h_max - h_min)
        return max(T_min, min(T_max, T_guess))


def test_mixing():
    """Unit test for stream mixing."""
    from .thermo import ThermoManager
    
    # Create test thermo manager
    thermo = ThermoManager("gri30.yaml")
    
    # Test case: all streams at same temperature
    streams = [
        {
            'm_dot': 0.1,
            'T': 1000.0,
            'P': 5e6,
            'Y': {'CH4': 0.5, 'O2': 0.5}
        },
        {
            'm_dot': 0.2,
            'T': 1000.0,
            'P': 5e6,
            'Y': {'H2O': 1.0}
        }
    ]
    
    T_mix, Y_mix, m_dot_total = mix_streams(streams, thermo)
    
    print(f"Test mixing:")
    print(f"  T_mix = {T_mix:.1f} K (expected: 1000.0 K)")
    print(f"  m_dot_total = {m_dot_total:.3f} kg/s (expected: 0.3 kg/s)")
    print(f"  Y_mix = {Y_mix}")
    
    # Verify temperature is correct
    assert abs(T_mix - 1000.0) < 1.0, f"Temperature mixing failed: {T_mix} != 1000.0"
    assert abs(m_dot_total - 0.3) < 1e-6, f"Mass flow mixing failed: {m_dot_total} != 0.3"
    
    print("✓ Mixing test passed")


if __name__ == "__main__":
    test_mixing()
