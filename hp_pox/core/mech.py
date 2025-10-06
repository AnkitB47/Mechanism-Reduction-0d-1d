"""
Mechanism loader and caching for HP-POX model.

Handles Cantera mechanism loading, transport model setup,
and provides cached access to thermodynamic and transport properties.
"""

import cantera as ct
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

warnings.filterwarnings('ignore', category=UserWarning)


class MechanismLoader:
    """Cantera mechanism loader with caching."""
    
    def __init__(self, mechanism_path: str):
        """Initialize mechanism loader.
        
        Args:
            mechanism_path: Path to mechanism file (.yaml or .cti)
        """
        self.mechanism_path = Path(mechanism_path)
        self.gas = None
        self.species_names = []
        self.species_indices = {}
        self.molecular_weights = np.array([])
        
        self._load_mechanism()
        self._setup_transport()
        self._cache_properties()
    
    def _load_mechanism(self) -> None:
        """Load Cantera mechanism."""
        if not self.mechanism_path.exists():
            raise FileNotFoundError(f"Mechanism file not found: {self.mechanism_path}")
        
        try:
            self.gas = ct.Solution(str(self.mechanism_path))
        except Exception as e:
            raise RuntimeError(f"Failed to load mechanism {self.mechanism_path}: {e}")
    
    def _setup_transport(self) -> None:
        """Setup transport model with fallbacks."""
        try:
            self.gas.transport_model = 'mixture-averaged'
        except Exception:
            try:
                self.gas.transport_model = 'multicomponent'
            except Exception:
                print("Warning: Using default transport model")
    
    def _cache_properties(self) -> None:
        """Cache frequently used properties."""
        self.species_names = [s.name for s in self.gas.species()]
        self.species_indices = {name: i for i, name in enumerate(self.species_names)}
        self.molecular_weights = np.array([s.molecular_weight for s in self.gas.species()])
    
    def get_species_index(self, species_name: str) -> int:
        """Get species index by name."""
        if species_name not in self.species_indices:
            raise KeyError(f"Species '{species_name}' not found in mechanism")
        return self.species_indices[species_name]
    
    def get_molecular_weight(self, species_name: str) -> float:
        """Get molecular weight of species."""
        idx = self.get_species_index(species_name)
        return self.molecular_weights[idx]
    
    def create_gas_state(self, temperature: float, pressure: float, 
                        mass_fractions: np.ndarray) -> 'GasState':
        """Create gas state object."""
        return GasState(self.gas, temperature, pressure, mass_fractions)
    
    def get_element_names(self) -> List[str]:
        """Get list of element names."""
        return self.gas.element_names
    
    def get_atomic_weight(self, element: str) -> float:
        """Get atomic weight of element."""
        return self.gas.atomic_weight(element)


class GasState:
    """Gas state wrapper with cached properties."""
    
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
        self.temperature = temperature
        self.pressure = pressure
        self.mass_fractions = mass_fractions.copy()
        
        # Set gas state
        self.gas.TPY = temperature, pressure, mass_fractions
        
        # Cache properties
        self._density = None
        self._mole_fractions = None
        self._cp_mass = None
        self._viscosity = None
        self._thermal_conductivity = None
        self._net_production_rates = None
        self._partial_molar_enthalpies = None
    
    @property
    def density(self) -> float:
        """Get density (kg/m³)."""
        if self._density is None:
            self._density = self.gas.density
        return self._density
    
    @property
    def mole_fractions(self) -> np.ndarray:
        """Get mole fractions."""
        if self._mole_fractions is None:
            self._mole_fractions = self.gas.X.copy()
        return self._mole_fractions
    
    @property
    def cp_mass(self) -> float:
        """Get specific heat capacity (J/kg/K)."""
        if self._cp_mass is None:
            self._cp_mass = self.gas.cp_mass
        return self._cp_mass
    
    @property
    def viscosity(self) -> float:
        """Get viscosity (Pa·s)."""
        if self._viscosity is None:
            try:
                self._viscosity = self.gas.viscosity
            except NotImplementedError:
                # Fallback to simple correlation
                self._viscosity = 1.8e-5 * (self.temperature / 300.0)**0.7
        return self._viscosity
    
    @property
    def thermal_conductivity(self) -> float:
        """Get thermal conductivity (W/m/K)."""
        if self._thermal_conductivity is None:
            try:
                self._thermal_conductivity = self.gas.thermal_conductivity
            except NotImplementedError:
                # Fallback to simple correlation
                self._thermal_conductivity = 0.025 * (self.temperature / 300.0)**0.7
        return self._thermal_conductivity
    
    @property
    def net_production_rates(self) -> np.ndarray:
        """Get net production rates (kmol/m³/s)."""
        if self._net_production_rates is None:
            self._net_production_rates = self.gas.net_production_rates.copy()
        return self._net_production_rates
    
    @property
    def partial_molar_enthalpies(self) -> np.ndarray:
        """Get partial molar enthalpies (J/kmol)."""
        if self._partial_molar_enthalpies is None:
            self._partial_molar_enthalpies = self.gas.partial_molar_enthalpies.copy()
        return self._partial_molar_enthalpies
    
    def get_elemental_mass_fraction(self, element: str) -> float:
        """Get elemental mass fraction."""
        return self.gas.elemental_mass_fraction(element)
    
    def update_state(self, temperature: float, pressure: float, mass_fractions: np.ndarray) -> None:
        """Update gas state and clear cache."""
        self.temperature = temperature
        self.pressure = pressure
        self.mass_fractions = mass_fractions.copy()
        
        # Set gas state
        self.gas.TPY = temperature, pressure, mass_fractions
        
        # Clear cache
        self._density = None
        self._mole_fractions = None
        self._cp_mass = None
        self._viscosity = None
        self._thermal_conductivity = None
        self._net_production_rates = None
        self._partial_molar_enthalpies = None
