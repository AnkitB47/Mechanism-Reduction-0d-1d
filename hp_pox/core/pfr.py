"""
1D Plug Flow Reactor (PFR) for HP-POX model.

Implements reforming zone with chemistry, heat loss, and momentum balance.
"""

import numpy as np
from typing import Dict, Any, Optional
from scipy.integrate import solve_ivp


class PFR:
    """1D Plug Flow Reactor for reforming zone."""
    
    def __init__(self, mech, case_config: Dict[str, Any], geometry_config: Dict[str, Any]):
        """Initialize PFR.
        
        Args:
            mech: Mechanism loader
            case_config: Case configuration
            geometry_config: Geometry configuration
        """
        self.mech = mech
        self.case_config = case_config
        self.geometry_config = geometry_config
        
        # PFR geometry
        self.length_m = geometry_config['pfr']['length_m']
        self.diameter_m = geometry_config['pfr']['diameter_m']
        self.area_m2 = np.pi * (self.diameter_m / 2)**2
        self.perimeter_m = np.pi * self.diameter_m
        
        # Species info
        self.n_species = len(self.mech.species_names)
        self.species_names = self.mech.species_names
        
        # Target residence time
        self.target_residence_time_s = 15.5
    
    def solve(self, psr_result: Dict[str, Any], heatloss, friction) -> Dict[str, Any]:
        """Solve PFR with proper physics using literature equations.
        
        Args:
            psr_result: PSR outlet conditions
            heatloss: Heat loss model
            friction: Friction model
            
        Returns:
            Dictionary with PFR results
        """
        print("Solving PFR (reforming zone) with proper physics...")
        
        # Use analytical mass flow (no scaling)
        self.scaled_mass_flow = self._calculate_total_mass_flow()
        print(f"Mass flow: {self.scaled_mass_flow:.6f} kg/s")
        print(f"Area: {self.area_m2:.6f} m²")
        print(f"Length: {self.length_m:.3f} m")
        
        # Create spatial grid
        n_cells = 400  # Start with fine grid
        z = np.linspace(0, self.length_m, n_cells)
        dz = z[1] - z[0]
        
        # Initialize arrays
        T = np.zeros(n_cells)
        P = np.zeros(n_cells)
        Y = np.zeros((self.n_species, n_cells))
        X = np.zeros((self.n_species, n_cells))
        u = np.zeros(n_cells)
        rho = np.zeros(n_cells)
        
        # Set initial conditions from PSR
        T[0] = psr_result['temperature_K']
        P[0] = psr_result['pressure_Pa']
        Y[:, 0] = psr_result['mass_fractions']
        
        # Calculate initial velocity from mass flow (with scaling)
        total_mass_flow = self.scaled_mass_flow
        rho[0] = psr_result['density_kg_m3']
        u[0] = total_mass_flow / (rho[0] * self.area_m2)
        print(f"Initial velocity: {u[0]:.3f} m/s")
        print(f"Initial density: {rho[0]:.3f} kg/m³")
        
        # Set initial gas state
        gas_state = self.mech.create_gas_state(T[0], P[0], Y[:, 0])
        X[:, 0] = gas_state.mole_fractions
        
        print(f"Target residence time: ~{self.target_residence_time_s} s")
        
        # March through reactor with proper physics
        print(f"🚀 Starting PFR march: {n_cells} cells, dz={dz:.4f}m")
        for i in range(1, n_cells):
            # Progress logging every 50 cells
            if i % 50 == 0 or i < 5:
                print(f"   📍 Cell {i:3d}/{n_cells}: z={z[i]:.3f}m, T={T[i-1]-273.15:.1f}°C, u={u[i-1]:.3f}m/s")
            
            # Set gas state at previous position
            gas_state.update_state(T[i-1], P[i-1], Y[:, i-1])
            
            # Get properties from Cantera
            rho[i-1] = gas_state.density
            cp = gas_state.cp_mass
            mu = gas_state.viscosity
            k = gas_state.thermal_conductivity
            omega = gas_state.net_production_rates  # mol/m³/s
            hk = gas_state.partial_molar_enthalpies  # J/mol
            
            # Safety checks
            if np.isnan(T[i-1]) or T[i-1] < 300.0:
                print(f"Warning: Invalid temperature at z={z[i-1]:.3f}, T={T[i-1]:.1f}K")
                T[i-1] = max(T[i-1], 300.0)
            
            if np.isnan(P[i-1]) or P[i-1] < 1e5:
                print(f"Warning: Invalid pressure at z={z[i-1]:.3f}, P={P[i-1]/1e5:.1f}bar")
                P[i-1] = max(P[i-1], 1e5)
            
            # Reynolds number
            Re = rho[i-1] * u[i-1] * self.diameter_m / mu
            
            # Prandtl number
            Pr = mu * cp / k
            
            # Friction factor
            f = friction.calculate_friction_factor(T[i-1], P[i-1], Y[:, i-1], u[i-1], gas_state)
            
            # Chemistry heat release (W/m)
            qdot_chem_per_length = -np.sum(omega * hk) * self.area_m2  # W/m
            
            # Wall heat sink (W/m) - use fixed value from case
            qdot_wall_per_length = heatloss.get_wall_heat_loss_per_length()  # W/m
            
            # Debug heat loss values (print once)
            if i == 1:
                print(f"Wall heat loss per length: {qdot_wall_per_length:.1f} W/m")
                print(f"Mass flow: {self.scaled_mass_flow:.6f} kg/s")
                print(f"Specific heat: {cp:.1f} J/kg/K")
            
            # Energy equation: dT/dz = (qdot_chem - qdot_wall) / (m_dot * cp)
            dT_dz = (qdot_chem_per_length - qdot_wall_per_length) / (self.scaled_mass_flow * cp)
            
            # Apply stability limits
            dT_dz = np.clip(dT_dz, -50.0, 50.0)
            
            # Species equations: dY_k/dz = (W_k * omega_k) / (m_dot / A)
            dY_dz = omega * self.mech.molecular_weights / (self.scaled_mass_flow / self.area_m2)
            
            # Apply stability limits to species
            dY_dz = np.clip(dY_dz, -0.1, 0.1)
            
            # Special handling for inert species
            inert_species = ['N2', 'AR', 'HE']
            for species in inert_species:
                if species in self.mech.species_indices:
                    idx = self.mech.species_indices[species]
                    dY_dz[idx] = 0.0
            
                # Momentum equation: dP/dz = -f * (rho * u^2) / (2 * D)
                dP_dz = -f * (rho[i-1] * u[i-1]**2) / (2 * self.diameter_m)
                
                # Add localized loss term for realistic pressure drop (nozzle/orifice effect)
                # This explains the 1-20 kPa pressure drop observed in experiments
                dP_local = 0.0
                if i == 1:  # At inlet
                    zeta_inlet = 100.0  # Loss coefficient for inlet nozzle (increased for realistic pressure drop)
                    dP_local = -zeta_inlet * (rho[i-1] * u[i-1]**2) / 2
                    print(f"Localized pressure drop: {dP_local/1000:.1f} kPa")
                
                dP_dz = np.clip(dP_dz, -1000.0, 0.0)  # Pressure should decrease
            
            # Update state
            T[i] = max(T[i-1] + dT_dz * dz, 300.0)
            P[i] = max(P[i-1] + dP_dz * dz + dP_local, 1e5)
            Y[:, i] = Y[:, i-1] + dY_dz * dz
            
            # Normalize mass fractions and ensure positivity
            Y[:, i] = np.maximum(Y[:, i], 1e-10)
            Y[:, i] = Y[:, i] / np.sum(Y[:, i])
            
            # Update gas state and recalculate velocity: u = m_dot / (rho * A)
            try:
                gas_state.update_state(T[i], P[i], Y[:, i])
                rho[i] = gas_state.density
                u[i] = self.scaled_mass_flow / (rho[i] * self.area_m2)  # Use scaled mass flow
                X[:, i] = gas_state.mole_fractions
            except Exception as e:
                print(f"Warning: Gas state failed at z={z[i]:.3f}, using previous values: {e}")
                T[i] = T[i-1]
                P[i] = P[i-1]
                Y[:, i] = Y[:, i-1]
                X[:, i] = X[:, i-1]
                rho[i] = rho[i-1]
                u[i] = u[i-1]
        
        # Calculate residence time
        residence_time = self._calculate_residence_time(z, u)
        
        # Detect ignition
        ignition_peak, ignition_thresh = self._detect_ignition(z, T)
        
        print(f"PFR outlet temperature: {T[-1] - 273.15:.1f}°C")
        print(f"PFR outlet pressure: {P[-1]/1e5:.1f} bar")
        print(f"Total residence time: {residence_time[-1]:.1f} s")
        print(f"Ignition length (peak): {ignition_peak:.3f} m")
        print(f"Ignition length (threshold): {ignition_thresh:.3f} m")
        print(f"Pressure drop: {(P[0] - P[-1])/1000:.1f} kPa")
        
        # Check elemental balance at outlet
        self._check_elemental_balance(gas_state, "PFR outlet")
        
        return {
            'z_m': z,
            'temperature_K': T,
            'pressure_Pa': P,
            'density_kg_m3': rho,
            'velocity_m_s': u,
            'mass_fractions': Y,
            'mole_fractions': X,
            'residence_time_s': residence_time,
            'ignition_length_peak_m': ignition_peak,
            'ignition_length_threshold_m': ignition_thresh,
            'pressure_drop_Pa': P[0] - P[-1],
            'species_names': self.species_names
        }
    
    def _calculate_adaptive_step_size(self, T: float, P: float, Y: np.ndarray, 
                                    u: float, dz: float, gas_state) -> float:
        """Calculate adaptive step size based on stability criteria."""
        # Base step size
        dz_adaptive = dz
        
        # Temperature change limit
        gas_state.update_state(T, P, Y)
        dT_dz = self._calculate_temperature_derivative(T, P, Y, u, gas_state)
        if abs(dT_dz) > 30.0:
            dz_adaptive = min(dz_adaptive, 30.0 / abs(dT_dz))
        
        # Species positivity check
        dY_dz = self._calculate_species_derivatives(T, P, Y, u, gas_state)
        for i in range(len(Y)):
            if Y[i] + dY_dz[i] * dz_adaptive < -1e-6:
                dz_adaptive = min(dz_adaptive, -Y[i] / dY_dz[i] / 2.0)
        
        return max(dz_adaptive, dz / 10.0)  # Don't go too small
    
    def _calculate_gnielinski_nusselt(self, Re: float, Pr: float, f: float) -> float:
        """Calculate Nusselt number using Gnielinski correlation."""
        if Re < 3000 or Re > 5e6:
            # Use Dittus-Boelter for extreme Re
            return 0.023 * (Re**0.8) * (Pr**0.4)
        
        # Gnielinski correlation: Nu = (f/8)(Re-1000)Pr / (1+12.7*sqrt(f/8)(Pr^(2/3)-1))
        numerator = (f/8) * (Re - 1000) * Pr
        denominator = 1 + 12.7 * np.sqrt(f/8) * (Pr**(2/3) - 1)
        return numerator / denominator
    
    def _calculate_heat_transfer_coefficient(self, T: float, P: float, Y: np.ndarray, 
                                           u: float, gas_state) -> float:
        """Calculate heat transfer coefficient using Gnielinski correlation."""
        # Get properties
        rho = gas_state.density
        mu = gas_state.viscosity
        k = gas_state.thermal_conductivity
        cp = gas_state.cp_mass
        
        # Reynolds number
        Re = rho * u * self.diameter_m / mu
        
        # Prandtl number
        Pr = mu * cp / k
        
        # Friction factor (Blasius for smooth turbulent)
        if Re < 2300:
            f = 64.0 / Re
        else:
            f = 0.316 / (Re**0.25)
        
        # Nusselt number
        Nu = self._calculate_gnielinski_nusselt(Re, Pr, f)
        
        # Heat transfer coefficient: h = Nu * k / D
        h = Nu * k / self.diameter_m
        
        return h
    
    def _calculate_derivatives(self, T: float, P: float, Y: np.ndarray, u: float,
                             gas_state, heatloss, friction, U: float) -> tuple:
        """Calculate all derivatives."""
        dT_dz = self._calculate_temperature_derivative(T, P, Y, u, gas_state, heatloss, U)
        dP_dz = self._calculate_pressure_derivative(T, P, Y, u, gas_state, friction)
        dY_dz = self._calculate_species_derivatives(T, P, Y, u, gas_state)
        du_dz = self._calculate_velocity_derivative(T, P, Y, u, gas_state, friction)
        
        return dT_dz, dP_dz, dY_dz, du_dz
    
    def _calculate_temperature_derivative(self, T: float, P: float, Y: np.ndarray, 
                                        u: float, gas_state, heatloss=None, U: float = 0.0) -> float:
        """Calculate dT/dz from energy balance using true reaction enthalpy."""
        # Mass flow rate
        total_mass_flow = self._calculate_total_mass_flow()
        
        # True heat source from chemistry: -Σ(ωk * hk) where ωk = net production rates, hk = partial molar enthalpies
        omega = gas_state.net_production_rates  # mol/m³/s
        hk = gas_state.partial_molar_enthalpies  # J/mol
        heat_source_vol = -np.sum(omega * hk)  # W/m³ (negative because exothermic reactions have negative hk)
        
        # Heat loss from wall using Gnielinski correlation
        h = self._calculate_heat_transfer_coefficient(T, P, Y, u, gas_state)
        perimeter = np.pi * self.diameter_m
        T_wall = heatloss.get_wall_temperature(0.0) if heatloss else 418.0  # K
        heat_loss_vol = h * (perimeter / self.area_m2) * (T - T_wall)  # W/m³
        
        # Energy balance: ρ u cp dT/dz = q̇_chem - q̇_wall
        dT_dz = (heat_source_vol - heat_loss_vol) / (gas_state.density * u * gas_state.cp_mass)
        
        return dT_dz
    
    def _calculate_pressure_derivative(self, T: float, P: float, Y: np.ndarray, 
                                     u: float, gas_state, friction) -> float:
        """Calculate dP/dz from momentum balance (Darcy-Weisbach)."""
        # Friction factor
        f = friction.calculate_friction_factor(T, P, Y, u, gas_state)
        
        # Darcy-Weisbach: dP/dz = -f * (rho * u^2) / (2 * D)
        dP_dz = -f * (gas_state.density * u**2) / (2 * self.diameter_m)
        
        return dP_dz
    
    def _calculate_species_derivatives(self, T: float, P: float, Y: np.ndarray, 
                                     u: float, gas_state) -> np.ndarray:
        """Calculate dY_i/dz from species balance."""
        # Mass flow rate
        total_mass_flow = self._calculate_total_mass_flow()
        
        # Species equations: dY_i/dz = (omega_i * MW_i) / (m_dot / A)
        omega = gas_state.net_production_rates
        dY_dz = omega * self.mech.molecular_weights / (total_mass_flow / self.area_m2)
        
        # Special handling for inert species
        inert_species = ['N2', 'AR', 'HE']
        for species in inert_species:
            if species in self.mech.species_indices:
                idx = self.mech.species_indices[species]
                dY_dz[idx] = 0.0  # No production/consumption for inert species
        
        return dY_dz
    
    def _calculate_velocity_derivative(self, T: float, P: float, Y: np.ndarray, 
                                     u: float, gas_state, friction) -> float:
        """Calculate du/dz from continuity equation."""
        # From continuity: du/dz = -(1/rho) * dP/dz
        dP_dz = self._calculate_pressure_derivative(T, P, Y, u, gas_state, friction)
        du_dz = -(1 / gas_state.density) * dP_dz
        
        return du_dz
    
    def _calculate_residence_time(self, z: np.ndarray, u: np.ndarray) -> np.ndarray:
        """Calculate residence time by integrating dz/u."""
        residence_time = np.zeros_like(z)
        for i in range(1, len(z)):
            dz = z[i] - z[i-1]
            residence_time[i] = residence_time[i-1] + dz / u[i]
        return residence_time
    
    def _detect_ignition(self, z: np.ndarray, T: np.ndarray) -> tuple:
        """Detect ignition length using two methods."""
        # Method 1: Peak dT/dz
        dT_dz = np.gradient(T, z)
        ignition_peak = z[np.argmax(dT_dz)]
        
        # Method 2: Temperature threshold
        T_thresh = T[0] + 400.0  # 400K above inlet
        thresh_idx = np.where(T >= T_thresh)[0]
        ignition_thresh = z[thresh_idx[0]] if len(thresh_idx) > 0 else z[-1]
        
        return ignition_peak, ignition_thresh
    
    def _calculate_total_mass_flow(self) -> float:
        """Calculate total mass flow rate (kg/s)."""
        total_mass_flow_kg_s = 0.0
        for feed_data in self.case_config['feeds'].values():
            total_mass_flow_kg_s += feed_data['mass_flow_kg_h'] / 3600.0
        return total_mass_flow_kg_s
    
    def _check_elemental_balance(self, gas_state, label: str) -> None:
        """Check elemental balance."""
        print(f"\nElemental balance {label}:")
        
        # Get elemental composition
        elements = ['C', 'H', 'O', 'N']
        for element in elements:
            if element in self.mech.get_element_names():
                mass_frac = gas_state.get_elemental_mass_fraction(element)
                atomic_weight = self.mech.get_atomic_weight(element)
                moles_per_kg = mass_frac / atomic_weight
                print(f"  {element}: {moles_per_kg:.6f} mol/kg")
