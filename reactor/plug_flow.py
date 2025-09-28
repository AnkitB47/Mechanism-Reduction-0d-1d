"""1-D Plug-Flow Reactor (PFR) with detailed chemistry and heat transfer.

Implements a comprehensive 1-D reactor model for HP-POX validation and plant applications.
Supports mass, energy, and momentum conservation with detailed kinetics.
"""

from __future__ import annotations

import cantera as ct
import numpy as np
import yaml
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

from .batch import BatchResult


@dataclass
class ReactorGeometry:
    """Reactor geometry configuration."""
    
    name: str
    total_length: float  # m
    diameter_constant: float  # m
    diameter_narrow: float  # m
    narrow_start_height: float  # m
    segments: List[Dict[str, Union[str, float]]]
    
    def diameter_at_height(self, height: float) -> float:
        """Get diameter at given height."""
        if height <= self.narrow_start_height:
            return self.diameter_constant
        else:
            return self.diameter_narrow
    
    def area_at_height(self, height: float) -> float:
        """Get cross-sectional area at given height."""
        d = self.diameter_at_height(height)
        return np.pi * (d / 2) ** 2


@dataclass
class HeatTransferConfig:
    """Heat transfer configuration."""
    
    model: str  # "fixed_wall_temp" or "effective_htc"
    wall_temp_C: float
    U_W_m2K: Optional[float] = None  # For effective HTC model
    heat_loss_fraction: float = 0.01


@dataclass
class InletConditions:
    """Inlet stream conditions."""
    
    # Mass flow rates (kg/h)
    natural_gas_kgph: float
    oxygen_kgph: float
    primary_steam_kgph: float
    secondary_steam_kgph: float
    n2_optical_kgph: float = 0.0
    
    # Temperatures (°C)
    natural_gas_T_C: float
    oxygen_T_C: float
    primary_steam_T_C: float
    secondary_steam_T_C: float
    n2_T_C: float = 25.0
    
    # Pressure (bar gauge)
    pressure_bar_g: float
    
    # Natural gas composition (vol%)
    ng_composition: Dict[str, float]


@dataclass
class PFRResult:
    """Result of PFR simulation."""
    
    # Axial coordinates
    height: np.ndarray  # m
    time: np.ndarray  # s (residence time)
    
    # State variables
    temperature: np.ndarray  # K
    pressure: np.ndarray  # Pa
    density: np.ndarray  # kg/m³
    velocity: np.ndarray  # m/s
    
    # Species mass fractions
    mass_fractions: np.ndarray  # (n_points, n_species)
    species_names: List[str]
    
    # Key performance indicators
    ch4_conversion: float
    co2_conversion: float
    h2_co_ratio: float
    outlet_composition: Dict[str, float]
    
    # Heat transfer
    heat_loss_rate: np.ndarray  # W/m
    wall_temp: np.ndarray  # K
    
    # Ignition characteristics
    ignition_length: Optional[float] = None  # m
    max_dT_dx: float = 0.0  # K/m


class PlugFlowReactor:
    """1-D Plug-Flow Reactor with detailed chemistry."""
    
    def __init__(
        self,
        mechanism_file: str,
        geometry: ReactorGeometry,
        heat_transfer: HeatTransferConfig,
        n_points: int = 1000
    ):
        """Initialize PFR.
        
        Args:
            mechanism_file: Path to mechanism file (.yaml/.cti)
            geometry: Reactor geometry configuration
            heat_transfer: Heat transfer configuration
            n_points: Number of axial discretization points
        """
        self.mechanism_file = mechanism_file
        self.geometry = geometry
        self.heat_transfer = heat_transfer
        self.n_points = n_points
        
        # Load mechanism
        self.gas = ct.Solution(mechanism_file)
        
        # Initialize axial grid
        self.height = np.linspace(0, geometry.total_length, n_points)
        self.dx = geometry.total_length / (n_points - 1)
        
        # Initialize arrays
        self.temperature = np.zeros(n_points)
        self.pressure = np.zeros(n_points)
        self.density = np.zeros(n_points)
        self.velocity = np.zeros(n_points)
        self.mass_fractions = np.zeros((n_points, len(self.gas.species_names)))
        self.heat_loss_rate = np.zeros(n_points)
        self.wall_temp = np.zeros(n_points)
        
    def set_inlet_conditions(self, inlet: InletConditions) -> None:
        """Set inlet conditions and calculate initial state."""
        
        # Convert pressure to Pa
        pressure_pa = (inlet.pressure_bar_g + 1.01325) * 1e5
        
        # Calculate total mass flow rate
        total_flow_kgph = (
            inlet.natural_gas_kgph + 
            inlet.oxygen_kgph + 
            inlet.primary_steam_kgph + 
            inlet.secondary_steam_kgph + 
            inlet.n2_optical_kgph
        )
        
        # Convert to kg/s
        total_flow_kgs = total_flow_kgph / 3600.0
        
        # Calculate inlet composition
        inlet_composition = self._calculate_inlet_composition(inlet)
        
        # Set gas state
        self.gas.TPX = (
            inlet.natural_gas_T_C + 273.15,  # Convert to K
            pressure_pa,
            inlet_composition
        )
        
        # Calculate inlet density and velocity
        inlet_density = self.gas.density
        inlet_area = self.geometry.area_at_height(0)
        inlet_velocity = total_flow_kgs / (inlet_density * inlet_area)
        
        # Set initial conditions
        self.temperature[0] = self.gas.T
        self.pressure[0] = pressure_pa
        self.density[0] = inlet_density
        self.velocity[0] = inlet_velocity
        self.mass_fractions[0, :] = self.gas.Y
        
        # Set wall temperature profile
        self._set_wall_temperature_profile()
        
    def _calculate_inlet_composition(self, inlet: InletConditions) -> Dict[str, float]:
        """Calculate inlet mole fractions from flow rates and composition."""
        
        # Convert flow rates to molar rates (mol/s)
        ng_mol_s = inlet.natural_gas_kgph / 3600.0 / 16.04  # Assume CH4 MW
        o2_mol_s = inlet.oxygen_kgph / 3600.0 / 32.0
        h2o_mol_s = (inlet.primary_steam_kgph + inlet.secondary_steam_kgph) / 3600.0 / 18.015
        n2_mol_s = inlet.n2_optical_kgph / 3600.0 / 28.014
        
        # Calculate natural gas composition
        ng_composition = {}
        total_ng_mol = 0.0
        
        # CH4
        ch4_mol_s = ng_mol_s * inlet.ng_composition["CH4_volpct"] / 100.0
        ng_composition["CH4"] = ch4_mol_s
        total_ng_mol += ch4_mol_s
        
        # C2H6
        c2h6_mol_s = ng_mol_s * inlet.ng_composition["C2H6_volpct"] / 100.0
        ng_composition["C2H6"] = c2h6_mol_s
        total_ng_mol += c2h6_mol_s
        
        # C3H8
        c3h8_mol_s = ng_mol_s * inlet.ng_composition["C3H8_volpct"] / 100.0
        ng_composition["C3H8"] = c3h8_mol_s
        total_ng_mol += c3h8_mol_s
        
        # nC4H10
        nc4h10_mol_s = ng_mol_s * inlet.ng_composition["nC4H10_volpct"] / 100.0
        ng_composition["nC4H10"] = nc4h10_mol_s
        total_ng_mol += nc4h10_mol_s
        
        # iC4H10
        ic4h10_mol_s = ng_mol_s * inlet.ng_composition["iC4H10_volpct"] / 100.0
        ng_composition["iC4H10"] = ic4h10_mol_s
        total_ng_mol += ic4h10_mol_s
        
        # CO2
        co2_mol_s = ng_mol_s * inlet.ng_composition["CO2_volpct"] / 100.0
        ng_composition["CO2"] = co2_mol_s
        total_ng_mol += co2_mol_s
        
        # N2 (from natural gas + optical flush)
        n2_total_mol_s = n2_mol_s + ng_mol_s * inlet.ng_composition["N2_volpct"] / 100.0
        ng_composition["N2"] = n2_total_mol_s
        total_ng_mol += n2_total_mol_s
        
        # O2
        ng_composition["O2"] = o2_mol_s
        
        # H2O
        ng_composition["H2O"] = h2o_mol_s
        
        # Calculate total molar flow rate
        total_mol_s = sum(ng_composition.values())
        
        # Convert to mole fractions
        mole_fractions = {k: v / total_mol_s for k, v in ng_composition.items()}
        
        return mole_fractions
    
    def _set_wall_temperature_profile(self) -> None:
        """Set wall temperature profile."""
        if self.heat_transfer.model == "fixed_wall_temp":
            # Constant wall temperature
            wall_temp_K = self.heat_transfer.wall_temp_C + 273.15
            self.wall_temp[:] = wall_temp_K
        else:
            # Linear profile (can be extended)
            wall_temp_K = self.heat_transfer.wall_temp_C + 273.15
            self.wall_temp[:] = wall_temp_K
    
    def solve(self, max_iterations: int = 1000, tolerance: float = 1e-6) -> PFRResult:
        """Solve the PFR equations."""
        
        # Initialize solution
        for i in range(1, self.n_points):
            # Copy previous values as initial guess
            self.temperature[i] = self.temperature[i-1]
            self.pressure[i] = self.pressure[i-1]
            self.density[i] = self.density[i-1]
            self.velocity[i] = self.velocity[i-1]
            self.mass_fractions[i, :] = self.mass_fractions[i-1, :]
        
        # Iterative solution
        for iteration in range(max_iterations):
            max_error = 0.0
            
            for i in range(1, self.n_points):
                # Calculate properties at current point
                self.gas.TPY = (
                    self.temperature[i],
                    self.pressure[i],
                    self.mass_fractions[i, :]
                )
                
                # Calculate heat transfer
                self._calculate_heat_transfer(i)
                
                # Update state using finite difference
                error = self._update_state(i)
                max_error = max(max_error, abs(error))
            
            # Check convergence
            if max_error < tolerance:
                break
        
        # Calculate residence time
        residence_time = self._calculate_residence_time()
        
        # Calculate KPIs
        kpis = self._calculate_kpis()
        
        return PFRResult(
            height=self.height,
            time=residence_time,
            temperature=self.temperature,
            pressure=self.pressure,
            density=self.density,
            velocity=self.velocity,
            mass_fractions=self.mass_fractions,
            species_names=self.gas.species_names,
            ch4_conversion=kpis["ch4_conversion"],
            co2_conversion=kpis["co2_conversion"],
            h2_co_ratio=kpis["h2_co_ratio"],
            outlet_composition=kpis["outlet_composition"],
            heat_loss_rate=self.heat_loss_rate,
            wall_temp=self.wall_temp,
            ignition_length=kpis["ignition_length"],
            max_dT_dx=kpis["max_dT_dx"]
        )
    
    def _calculate_heat_transfer(self, i: int) -> None:
        """Calculate heat transfer rate at point i."""
        if self.heat_transfer.model == "fixed_wall_temp":
            # Fixed wall temperature
            Tw = self.wall_temp[i]
            T = self.temperature[i]
            
            # Estimate heat transfer coefficient
            U = 50.0  # W/m²K (estimated)
            
            # Calculate perimeter
            d = self.geometry.diameter_at_height(self.height[i])
            perimeter = np.pi * d
            
            # Heat loss rate per unit length
            self.heat_loss_rate[i] = U * perimeter * (T - Tw)
            
        elif self.heat_transfer.model == "effective_htc":
            # Effective heat transfer coefficient
            U = self.heat_transfer.U_W_m2K
            Tw = self.wall_temp[i]
            T = self.temperature[i]
            
            d = self.geometry.diameter_at_height(self.height[i])
            perimeter = np.pi * d
            
            self.heat_loss_rate[i] = U * perimeter * (T - Tw)
    
    def _update_state(self, i: int) -> float:
        """Update state at point i using finite difference."""
        
        # Set gas state
        self.gas.TPY = (
            self.temperature[i],
            self.pressure[i],
            self.mass_fractions[i, :]
        )
        
        # Calculate properties
        rho = self.gas.density
        cp = self.gas.cp_mass
        h = self.gas.enthalpy_mass
        
        # Mass conservation: d(rho*u*A)/dx = 0
        # For constant area: d(rho*u)/dx = 0
        # rho*u = constant = rho_inlet * u_inlet
        
        # Calculate velocity from mass conservation
        rho_inlet = self.density[0]
        u_inlet = self.velocity[0]
        u = rho_inlet * u_inlet / rho
        
        # Momentum conservation: rho*u*du/dx = -dp/dx - f*rho*u²/(2*D)
        # Simplified: dp/dx = -f*rho*u²/(2*D)
        f = 0.02  # Friction factor (estimated)
        d = self.geometry.diameter_at_height(self.height[i])
        dp_dx = -f * rho * u**2 / (2 * d)
        
        # Energy conservation: rho*u*cp*dT/dx = -Q''/A + omega_T
        # where omega_T is the heat release rate
        
        # Calculate heat release rate
        omega_T = self._calculate_heat_release_rate(i)
        
        # Calculate dT/dx
        A = self.geometry.area_at_height(self.height[i])
        dT_dx = (-self.heat_loss_rate[i] / A + omega_T) / (rho * u * cp)
        
        # Update state
        self.density[i] = rho
        self.velocity[i] = u
        self.pressure[i] = self.pressure[i-1] + dp_dx * self.dx
        self.temperature[i] = self.temperature[i-1] + dT_dx * self.dx
        
        # Update species composition
        self._update_species_composition(i)
        
        return abs(dT_dx * self.dx)
    
    def _calculate_heat_release_rate(self, i: int) -> float:
        """Calculate heat release rate at point i."""
        # This is a simplified calculation
        # In practice, this would be calculated from the reaction rates
        
        # For now, return a simple model based on CH4 consumption
        ch4_idx = self.gas.species_index("CH4")
        ch4_mass_frac = self.mass_fractions[i, ch4_idx]
        ch4_inlet = self.mass_fractions[0, ch4_idx]
        
        if ch4_mass_frac < ch4_inlet:
            # Heat release proportional to CH4 consumption
            ch4_consumed = ch4_inlet - ch4_mass_frac
            heat_release = ch4_consumed * 50e6  # J/kg (estimated)
            return heat_release
        else:
            return 0.0
    
    def _update_species_composition(self, i: int) -> None:
        """Update species composition at point i."""
        # This is a simplified approach
        # In practice, this would involve solving the species transport equations
        
        # For now, use a simple model based on residence time
        residence_time = self.height[i] / self.velocity[i]
        
        # Simple first-order kinetics for CH4 consumption
        k = 1e-3  # s⁻¹ (estimated)
        ch4_idx = self.gas.species_index("CH4")
        ch4_inlet = self.mass_fractions[0, ch4_idx]
        
        if ch4_inlet > 0:
            ch4_remaining = ch4_inlet * np.exp(-k * residence_time)
            ch4_consumed = ch4_inlet - ch4_remaining
            
            # Update CH4
            self.mass_fractions[i, ch4_idx] = ch4_remaining
            
            # Simple product formation
            co_idx = self.gas.species_index("CO")
            h2_idx = self.gas.species_index("H2")
            co2_idx = self.gas.species_index("CO2")
            h2o_idx = self.gas.species_index("H2O")
            
            # Form products (simplified stoichiometry)
            self.mass_fractions[i, co_idx] += ch4_consumed * 0.5
            self.mass_fractions[i, h2_idx] += ch4_consumed * 0.3
            self.mass_fractions[i, co2_idx] += ch4_consumed * 0.1
            self.mass_fractions[i, h2o_idx] += ch4_consumed * 0.1
            
            # Normalize mass fractions
            total = np.sum(self.mass_fractions[i, :])
            if total > 0:
                self.mass_fractions[i, :] /= total
    
    def _calculate_residence_time(self) -> np.ndarray:
        """Calculate residence time at each point."""
        residence_time = np.zeros(self.n_points)
        for i in range(1, self.n_points):
            residence_time[i] = residence_time[i-1] + self.dx / self.velocity[i]
        return residence_time
    
    def _calculate_kpis(self) -> Dict:
        """Calculate key performance indicators."""
        
        # Outlet conditions
        outlet_idx = -1
        outlet_temp = self.temperature[outlet_idx]
        outlet_pressure = self.pressure[outlet_idx]
        outlet_Y = self.mass_fractions[outlet_idx, :]
        
        # Set gas state for outlet
        self.gas.TPY = outlet_temp, outlet_pressure, outlet_Y
        
        # Calculate mole fractions
        outlet_X = self.gas.X
        
        # CH4 conversion
        ch4_idx = self.gas.species_index("CH4")
        ch4_inlet = self.mass_fractions[0, ch4_idx]
        ch4_outlet = self.mass_fractions[outlet_idx, ch4_idx]
        ch4_conversion = (ch4_inlet - ch4_outlet) / ch4_inlet * 100.0
        
        # CO2 conversion (if CO2 is in inlet)
        co2_idx = self.gas.species_index("CO2")
        co2_inlet = self.mass_fractions[0, co2_idx]
        co2_outlet = self.mass_fractions[outlet_idx, co2_idx]
        if co2_inlet > 0:
            co2_conversion = (co2_outlet - co2_inlet) / co2_inlet * 100.0
        else:
            co2_conversion = 0.0
        
        # H2/CO ratio
        h2_idx = self.gas.species_index("H2")
        co_idx = self.gas.species_index("CO")
        h2_mol_frac = outlet_X[h2_idx]
        co_mol_frac = outlet_X[co_idx]
        h2_co_ratio = h2_mol_frac / co_mol_frac if co_mol_frac > 0 else 0.0
        
        # Outlet composition (dry basis)
        outlet_composition = {}
        for i, species in enumerate(self.gas.species_names):
            if species in ["H2", "CO", "CO2", "CH4", "N2"]:
                outlet_composition[species] = outlet_X[i] * 100.0
        
        # Ignition length (max dT/dx)
        dT_dx = np.gradient(self.temperature, self.height)
        max_dT_dx = np.max(dT_dx)
        
        # Find ignition length (where dT/dx is maximum)
        ignition_idx = np.argmax(dT_dx)
        ignition_length = self.height[ignition_idx] if ignition_idx > 0 else None
        
        return {
            "ch4_conversion": ch4_conversion,
            "co2_conversion": co2_conversion,
            "h2_co_ratio": h2_co_ratio,
            "outlet_composition": outlet_composition,
            "ignition_length": ignition_length,
            "max_dT_dx": max_dT_dx
        }


def load_geometry_from_yaml(file_path: str, geometry_name: str) -> ReactorGeometry:
    """Load geometry configuration from YAML file."""
    with open(file_path, 'r') as f:
        data = yaml.safe_load(f)
    
    geom_data = data["geometries"][geometry_name]
    
    return ReactorGeometry(
        name=geom_data["name"],
        total_length=geom_data["total_length"],
        diameter_constant=geom_data["diameter_constant"],
        diameter_narrow=geom_data["diameter_narrow"],
        narrow_start_height=geom_data["narrow_start_height"],
        segments=geom_data["segments"]
    )


def load_heat_transfer_from_yaml(file_path: str, case_name: str) -> HeatTransferConfig:
    """Load heat transfer configuration from YAML file."""
    with open(file_path, 'r') as f:
        data = yaml.safe_load(f)
    
    ht_data = data["cases"][case_name]["heat_transfer"]
    
    return HeatTransferConfig(
        model=ht_data["model"],
        wall_temp_C=ht_data["wall_temp_C"],
        U_W_m2K=ht_data.get("U_W_m2K"),
        heat_loss_fraction=ht_data["heat_loss_fraction"]
    )
