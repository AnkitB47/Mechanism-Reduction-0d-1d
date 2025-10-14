"""
Thermodynamics and transport properties for HP-POX model.
Handles mechanism loading and gas state management.
"""

import cantera as ct
import numpy as np
from pathlib import Path
from typing import Dict, Any, Optional

# Enable KLU globally (if the build supports it); otherwise do nothing and let dense be used
try:
    import cantera.sundials as cs  # Cantera 3.x
    if hasattr(cs, "use_klu"):
        cs.use_klu(True)  # enable sparse/KLU globally if available
        print("  KLU sparse solver enabled globally")
except Exception:
    pass  # no KLU; keep default dense

# Threading & performance hygiene (no thread storms)
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# [ROBUST] ---- Common numerical guards ----
T_MIN = 250.0      # K  (prevents 0K setTemperature)
T_MAX = 4000.0     # K
Y_FLOOR = 1e-30
DT_MIN = 1e-12     # s  (prevents "tout too close to t0")
REL_DY_MAX = 0.10  # max 10% change/species per microstep
ABS_DY_MAX = 0.10
DT_GROWTH = 1.8
DT_SHRINK = 0.5

def yi_atol(Y):
    """Elementwise atols for batch/PSR chemistry vectors (floor avoids zero species stalling integration)"""
    return (1e-12 * np.maximum(1.0, Y)).clip(min=1e-18)

def clamp_state(T, Y, rho=None):
    """Hard physical clamps & NaN/Inf handling"""
    T_MIN, T_MAX = 250.0, 3500.0
    Y_SUM_TOL = 1e-12
    RHO_MIN = 1e-6
    
    T = float(np.clip(T, T_MIN, T_MAX))
    Y = np.clip(Y, 0.0, 1.0)
    s = Y.sum()
    if not np.isfinite(s) or s < Y_SUM_TOL:
        # fallback to inlet composition if normalization is impossible
        raise RuntimeError("Non-physical Y; abort chemistry for this cell")
    Y /= s
    
    if rho is not None:
        rho = max(float(rho), RHO_MIN)
        return T, Y, rho
    return T, Y

def is_finite_state(T, Y):
    """Global NaN/Inf barrier"""
    return np.isfinite(T) and np.all(np.isfinite(Y)) and 0.999 <= Y.sum() <= 1.001

def sanitize_state(T, P, Y):
    """[ROBUST] Clip T and Y; renormalize Y with comprehensive bounds checking"""
    # Temperature bounds: T ∈ [250, 4000] K
    if not np.isfinite(T) or T <= 0.0:
        T = T_MIN
    T = float(np.clip(T, T_MIN, T_MAX))
    
    # Pressure bounds: P > 1e-6 Pa
    P = float(max(P, 1e-6))
    
    # Species bounds: Yi ≥ 0, ∑Yi → 1 (renormalize)
    Y = np.asarray(Y, dtype=np.float64)  # Double precision
    Y[~np.isfinite(Y)] = 0.0
    Y = np.maximum(Y, Y_FLOOR)  # Species floor: Yi = max(Yi, 1e-30)
    Y /= Y.sum()  # Renormalize to ensure ∑Yi = 1
    
    return T, P, Y

def finite_or_backtrack(advance_fn, state_getter, dt_try):
    """
    [ROBUST] advance_fn(dt) must advance integrator by dt (positive).
    state_getter() -> (T,P,Y) AFTER advance.
    Returns (ok, dt_used).
    """
    dt = max(dt_try, DT_MIN)
    for _ in range(12):  # limited backtracks
        try:
            advance_fn(dt)
            T, P, Y = state_getter()
            if np.isfinite(T) and np.isfinite(P) and np.all(np.isfinite(Y)):
                return True, dt
        except Exception:
            pass
        dt *= DT_SHRINK
        if dt < DT_MIN:
            return False, dt
    return False, dt

def sanitize_TY(T, Y):
    """Sanitize temperature and mass fractions."""
    T = float(np.clip(T, T_MIN, T_MAX))
    Y = np.asarray(Y, float)
    Y[~np.isfinite(Y)] = 0.0
    Y[Y < 0] = 0.0
    s = Y.sum()
    if s <= 0:
        # keep last valid composition (caller provides)
        return T, None
    return T, (Y / s).clip(min=Y_FLOOR)


class ThermoManager:
    """Manages thermodynamic properties and mechanism loading."""
    
    def __init__(self, mechanism_path: str = "gri30.yaml", transport_model: str = "mixture-averaged"):
        """Initialize thermodynamics manager.
        
        Args:
            mechanism_path: Path to mechanism file (.yaml or .cti)
            transport_model: Transport model ('mixture-averaged' for laptop, 'multicomponent' for HPC)
        """
        self.mechanism_path = mechanism_path
        self.gas = self._load_mechanism(transport_model)
        self.species_names = self.gas.species_names
        self.n_species = len(self.species_names)
        self.molecular_weights = self.gas.molecular_weights
        
        # Create species index mapping
        self.species_indices = {name: i for i, name in enumerate(self.species_names)}
        
        # Transport fallback tracking
        self.transport_fallback_count = 0
        self.max_fallback_attempts = 4
        self.original_transport_model = self.gas.transport_model
        
        print(f"Loaded mechanism: {len(self.species_names)} species, {len(self.gas.reactions())} reactions")
        print(f"Transport model: {self.gas.transport_model}")
    
    def handle_transport_fallback(self, operation_name: str = "operation") -> bool:
        """Handle transport model fallback when integrator rejects >4 times in a row.
        
        Args:
            operation_name: Name of the operation for logging
            
        Returns:
            True if fallback was applied, False if already at fallback model
        """
        self.transport_fallback_count += 1
        
        if self.transport_fallback_count > self.max_fallback_attempts:
            if self.gas.transport_model == 'Multi':
                print(f"    Transport fallback: Switching from Multi to Mix for {operation_name}")
                self.gas.transport_model = 'Mix'
                return True
            else:
                print(f"    Warning: Already using fallback transport model for {operation_name}")
                return False
        else:
            return False
    
    def reset_transport_fallback(self):
        """Reset transport fallback counter and restore original model."""
        if self.transport_fallback_count > self.max_fallback_attempts:
            print(f"    Transport reset: Restoring {self.original_transport_model} model")
            self.gas.transport_model = self.original_transport_model
        self.transport_fallback_count = 0
    
    def _load_mechanism(self, transport_model: str = 'mixture-averaged') -> ct.Solution:
        """Load chemical mechanism with proper transport model (default: mixture-averaged for laptop)."""
        try:
            gas = ct.Solution(self.mechanism_path, transport_model=transport_model)
            print(f"  Loaded mechanism: {len(gas.species_names)} species, {len(gas.reactions())} reactions")
            print(f"  Transport model: {transport_model}")
            return gas
        except Exception as e:
            # Fallback to mixture-averaged if multicomponent fails
            if transport_model == 'multicomponent':
                print(f"  Warning: multicomponent transport failed, falling back to mixture-averaged")
                try:
                    gas = ct.Solution(self.mechanism_path, transport_model='mixture-averaged')
                    print(f"  Loaded mechanism: {len(gas.species_names)} species, {len(gas.reactions())} reactions")
                    print(f"  Transport model: mixture-averaged (fallback)")
                    return gas
                except Exception as e2:
                    raise RuntimeError(f"Failed to load mechanism {self.mechanism_path} with both transport models: {e2}")
            else:
                raise RuntimeError(f"Failed to load mechanism from {self.mechanism_path}: {e}")
    
    def create_gas_state(self, temperature: float, pressure: float, 
                        mass_fractions: np.ndarray) -> 'GasState':
        """Create a gas state object.
        
        Args:
            temperature: Temperature (K)
            pressure: Pressure (Pa)
            mass_fractions: Mass fractions array
            
        Returns:
            GasState object
        """
        return GasState(self.gas, temperature, pressure, mass_fractions)
    
    def get_species_index(self, species_name: str) -> int:
        """Get index of species by name."""
        return self.species_indices[species_name]


class GasState:
    """Represents a gas state with thermodynamic and transport properties."""
    
    def __init__(self, gas: ct.Solution, temperature: float, pressure: float, 
                 mass_fractions: np.ndarray):
        """Initialize gas state.
        
        Args:
            gas: Cantera gas object
            temperature: Temperature (K)
            pressure: Pressure (Pa)
            mass_fractions: Mass fractions array
        """
        self.gas = gas
        
        # [ROBUST] Sanitize and set state
        T, P, Y = sanitize_state(temperature, pressure, mass_fractions)
        self.gas.TPY = T, P, Y
        
        # Store sanitized values
        self.temperature = T
        self.pressure = P
        self.mass_fractions = Y.copy()
        
        # Get properties (always read from gas, no caching)
        self.density = self.gas.density
        self.mole_fractions = self.gas.X
        self.molecular_weights = self.gas.molecular_weights
        
        # Transport properties
        try:
            self.viscosity = self.gas.viscosity
            self.thermal_conductivity = self.gas.thermal_conductivity
        except:
            # Fallback approximations
            self.viscosity = 1.8e-5 * (T / 300.0)**0.7
            self.thermal_conductivity = 0.026 * (T / 300.0)**0.9
    
    @property
    def species_names(self):
        """Get species names list."""
        return list(self.gas.species_names)
    
    @property
    def Y(self):
        """[ROBUST] Get mass fractions array."""
        return self.gas.Y.copy()
    
    @property
    def T(self):
        """[ROBUST] Get temperature."""
        return float(self.gas.T)
    
    @property
    def P(self):
        """[ROBUST] Get pressure."""
        return float(self.gas.P)
    
    @property
    def cp_mass(self):
        """Get specific heat at constant pressure (J/kg/K)."""
        return float(self.gas.cp_mass)
    
    @property
    def cv_mass(self):
        """Get specific heat at constant volume (J/kg/K)."""
        return float(self.gas.cv_mass)
    
    @property
    def Rmix(self):
        """[ROBUST] Get mixture gas constant (J/kg/K)."""
        return ct.gas_constant / float(self.gas.mean_molecular_weight)
    
    def set_TPY(self, T, P, Y):
        """[ROBUST] Set T, P, Y with sanitization."""
        T, P, Y = sanitize_state(T, P, Y)
        self.gas.TPY = T, P, Y
        self.temperature = T
        self.pressure = P
        self.mass_fractions = Y.copy()
    
    def copy_state(self):
        """Copy current state."""
        return float(self.gas.T), float(self.gas.P), self.gas.Y.copy()
    
    def restore_state(self, T, P, Y):
        """Restore state with bounds checking."""
        self.gas.TPY = max(T, 250.0), max(P, 1.0), np.clip(Y, 0, 1)
        self.temperature = max(T, 250.0)
        self.pressure = max(P, 1.0)
        self.mass_fractions = np.clip(Y, 0, 1)
    
    def copy(self):
        """Create a copy of this gas state."""
        # Create new gas state with same properties
        # Use the original mechanism file path instead of gas.name
        new_gas = ct.Solution(self.gas.source)
        new_gas.TPY = self.gas.T, self.gas.P, self.gas.Y
        return GasState(new_gas, self.temperature, self.pressure, self.mass_fractions)
    
    def update_state(self, temperature: float, pressure: float, mass_fractions: np.ndarray):
        """Update gas state.
        
        Args:
            temperature: Temperature (K)
            pressure: Pressure (Pa)
            mass_fractions: Mass fractions array
        """
        # Use sanitized set_TPY
        self.set_TPY(temperature, pressure, mass_fractions)
        
        # Update properties (always read from gas)
        self.density = self.gas.density
        self.mole_fractions = self.gas.X
        
        # Transport properties
        try:
            self.viscosity = self.gas.viscosity
            self.thermal_conductivity = self.gas.thermal_conductivity
        except:
            # Fallback approximations
            self.viscosity = 1.8e-5 * (temperature / 300.0)**0.7
            self.thermal_conductivity = 0.026 * (temperature / 300.0)**0.9
    
    def get_net_production_rates(self) -> np.ndarray:
        """Get net production rates (mol/m³/s)."""
        return self.gas.net_production_rates
    
    def get_partial_molar_enthalpies(self) -> np.ndarray:
        """Get partial molar enthalpies (J/mol)."""
        return self.gas.partial_molar_enthalpies
    
    def get_mass_fraction(self, species_name: str) -> float:
        """Get mass fraction of species."""
        idx = self.gas.species_index(species_name)
        return self.mass_fractions[idx]
    
    def get_mole_fraction(self, species_name: str) -> float:
        """Get mole fraction of species."""
        idx = self.gas.species_index(species_name)
        return self.mole_fractions[idx]
    
    def get_dry_basis_composition(self) -> Dict[str, float]:
        """Get dry basis composition (excluding H2O only)."""
        dry_fractions = {}
        total_dry = 0.0
        
        # Calculate total excluding H2O
        for i, species in enumerate(self.gas.species_names):
            if species != 'H2O':
                dry_fractions[species] = self.mole_fractions[i]
                total_dry += self.mole_fractions[i]
        
        # Renormalize
        if total_dry > 0:
            for species in dry_fractions:
                dry_fractions[species] /= total_dry
        
        return dry_fractions
