"""
Friction models for HP-POX reactor.

Implements Darcy-Weisbach equation with Blasius and Colebrook correlations.
"""

import numpy as np
from typing import Dict, Any


def f_blasius(Re: float) -> float:
    """Blasius: f = 0.316 * Re^(-0.25) for smooth pipe, 4e3 < Re < 1e5."""
    if Re < 2300:
        return 64.0 / Re  # Laminar
    else:
        return 0.316 / (Re**0.25)  # Turbulent Blasius


def f_haaland(Re: float, rel_roughness: float = 0.0) -> float:
    """
    Haaland explicit approximation for turbulent pipe:
      1/sqrt(f) = -1.8 * log10( (rel_roughness/3.7)^1.11 + 6.9/Re )
    """
    if Re < 2300:
        return 64.0 / Re  # Laminar
    else:
        # Haaland approximation
        term1 = (rel_roughness/3.7)**1.11
        term2 = 6.9/Re
        log_term = np.log10(term1 + term2)
        f = 1.0 / (-1.8 * log_term)**2
        return f


class Friction:
    """Friction model for HP-POX reactor."""
    
    def __init__(self, mech):
        """Initialize friction model.
        
        Args:
            mech: Mechanism loader
        """
        self.mech = mech
    
    def calculate_friction_factor(self, T: float, P: float, Y: np.ndarray, 
                                u: float, gas_state) -> float:
        """Calculate friction factor using Darcy-Weisbach.
        
        Args:
            T: Temperature (K)
            P: Pressure (Pa)
            Y: Mass fractions
            u: Velocity (m/s)
            gas_state: Gas state object
            
        Returns:
            Friction factor
        """
        # Calculate Reynolds number
        Re = self._calculate_reynolds_number(T, P, Y, u, gas_state)
        
        # Calculate friction factor
        if Re < 2300:
            # Laminar flow
            f = 64.0 / Re
        else:
            # Turbulent flow - use Blasius correlation
            f = 0.316 / (Re**0.25)
        
        return f
    
    def _calculate_reynolds_number(self, T: float, P: float, Y: np.ndarray, 
                                 u: float, gas_state) -> float:
        """Calculate Reynolds number.
        
        Args:
            T: Temperature (K)
            P: Pressure (Pa)
            Y: Mass fractions
            u: Velocity (m/s)
            gas_state: Gas state object
            
        Returns:
            Reynolds number
        """
        # Get properties
        rho = gas_state.density
        mu = gas_state.viscosity
        
        # PFR diameter
        diameter = 0.49  # m (from geometry config)
        
        # Reynolds number: Re = rho * u * D / mu
        Re = rho * u * diameter / mu
        
        return Re
    
    def calculate_pressure_drop(self, z: np.ndarray, T: np.ndarray, P: np.ndarray, 
                              Y: np.ndarray, u: np.ndarray, gas_states) -> float:
        """Calculate total pressure drop.
        
        Args:
            z: Axial positions (m)
            T: Temperature profile (K)
            P: Pressure profile (Pa)
            Y: Mass fractions profile
            u: Velocity profile (m/s)
            gas_states: Gas states at each position
            
        Returns:
            Total pressure drop (Pa)
        """
        total_pressure_drop = 0.0
        
        for i in range(len(z) - 1):
            dz = z[i+1] - z[i]
            f = self.calculate_friction_factor(T[i], P[i], Y[:, i], u[i], gas_states[i])
            
            # Pressure drop for this segment
            rho = gas_states[i].density
            dp = f * (rho * u[i]**2) / (2 * 0.49) * dz  # 0.49 m diameter
            total_pressure_drop += dp
        
        return total_pressure_drop
