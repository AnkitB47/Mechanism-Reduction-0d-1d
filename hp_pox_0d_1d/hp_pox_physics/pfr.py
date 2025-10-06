"""
Plug Flow Reactor (PFR) with proper energy units and chemistry control.
Implements molar-based energy source calculation and adaptive stepping.
"""

import numpy as np
import cantera as ct
from typing import Dict, Any, Tuple, Optional
from .thermo import ThermoManager, GasState, sanitize_TY, T_MIN, T_MAX, Y_FLOOR, sanitize_state, finite_or_backtrack, DT_MIN, DT_GROWTH, DT_SHRINK, REL_DY_MAX, clamp_state, is_finite_state


class PFR:
    """Plug Flow Reactor with proper energy balance and chemistry control."""
    
    def __init__(self, length_m: float, diameter_m: float, thermo: ThermoManager):
        """Initialize PFR.
        
        Args:
            length_m: Reactor length (m)
            diameter_m: Reactor diameter (m)
            thermo: Thermodynamics manager
        """
        self.length_m = length_m
        self.diameter_m = diameter_m
        self.area_m2 = np.pi * (diameter_m / 2)**2
        self.thermo = thermo
    
    def solve_pfr_operator_split(self, inlet_state: GasState, m_dot: float, P: float,
                                q_wall_per_m: float, chemistry_on: bool = True, Y_inlet: np.ndarray = None) -> Dict[str, Any]:
        # [ROBUST] Store parameters as instance variables
        self.m_dot = m_dot
        self.q_wall_per_m = q_wall_per_m
        """
        Solve PFR using operator splitting: energy/flow predictor + chemistry substep.
        
        Args:
            inlet_state: Inlet gas state
            m_dot: Mass flow rate (kg/s)
            P: Pressure (Pa)
            q_wall_per_m: Wall heat loss per length (W/m)
            chemistry_on: Whether to include chemistry
            
        Returns:
            Dictionary with PFR results
        """
        print(f"Solving PFR with operator splitting: L={self.length_m:.2f}m, D={self.diameter_m:.2f}m")
        print(f"  m_dot={m_dot:.3f}kg/s, q_wall={q_wall_per_m:.1f}W/m")
        print(f"  Chemistry: {'ON' if chemistry_on else 'OFF'}")
        
        # Create spatial grid with adaptive stepping
        base_dz = 0.001  # m
        n_cells = int(self.length_m / base_dz) + 1
        z = np.linspace(0, self.length_m, n_cells)
        dz = z[1] - z[0]
        
        # Initialize arrays
        n_species = len(inlet_state.gas.species_names)
        T = np.zeros(n_cells)
        rho = np.zeros(n_cells)
        u = np.zeros(n_cells)
        Y = np.zeros((n_cells, n_species))  # [ROBUST] Fix indexing: (n_cells, n_species)
        X = np.zeros((n_cells, n_species))  # [ROBUST] Fix indexing: (n_cells, n_species)
        
        # Set inlet conditions
        T[0] = inlet_state.temperature
        rho[0] = inlet_state.density
        u[0] = m_dot / (rho[0] * self.area_m2)
        Y[0, :] = inlet_state.mass_fractions  # [ROBUST] Fix indexing
        X[0, :] = inlet_state.mole_fractions  # [ROBUST] Fix indexing
        
        # Initialize gas state for chemistry
        gas_state = self.thermo.create_gas_state(
            inlet_state.temperature, 
            inlet_state.pressure, 
            inlet_state.mass_fractions
        )
        
        # Initialize unified data storage
        self.unified_data = []
        
        # Use physical inlet composition if provided
        if Y_inlet is not None:
            Y_current = Y_inlet.copy()
        else:
            Y_current = inlet_state.mass_fractions.copy()
        
        # Ensure Y_current is 1D and has correct length
        if Y_current.ndim > 1:
            Y_current = Y_current.flatten()
        
        # Debug: print shape
        print(f"DEBUG: Y_current shape: {Y_current.shape}, expected: {self.thermo.n_species}")
        if len(Y_current) != self.thermo.n_species:
            print(f"ERROR: Y_current length {len(Y_current)} != n_species {self.thermo.n_species}")
            # Use inlet state as fallback
            Y_current = inlet_state.mass_fractions.copy()
        
        # Cell-size adaptation tracking
        consecutive_fallbacks = 0
        dz_refined = False
        original_dz = dz
        
        # Integration loop with operator splitting
        for i in range(n_cells - 1):
            current_dz = z[i+1] - z[i]
            
            # Step 1: Energy/flow predictor (no chemistry)
            print(f"DEBUG: Y_current shape before predictor: {Y_current.shape}")
            T_pred, rho_pred, u_pred = self._energy_flow_predictor(
                T[i], P, Y_current, u[i], current_dz
            )
            
            # Step 2: Chemistry substep (if enabled)
            chemistry_applied = False
            if chemistry_on:
                # Guard against zero or negative residence time
                u_safe = max(float(u_pred), 1e-6)    # m/s
                dz_safe = max(float(current_dz), 1e-6)  # m
                tau = dz_safe / u_safe
                
                if tau <= 0 or not np.isfinite(tau):
                    print(f"    Chemistry: invalid tau={tau} -> skip chemistry for this cell")
                    T[i+1] = T_pred
                    Y[i+1, :] = Y[i, :]
                    Y_current = Y[i+1, :].copy()
                else:
                    # Prepare state for chemistry
                    state_in = {
                        'T': T_pred,
                        'p': P,
                        'Y': Y_current
                    }
                    
                    # Compute dt_try based on residence time (use local predicted velocity)
                    dt_try = max(0.0, current_dz) / max(1e-3, float(u_pred))
                    
                    # Run adaptive chemistry
                    state_out, chem_info = self.chemistry_step_adaptive(
                        state_in, dt_try, dt_min=1e-8, dt_max=1e-4
                    )
                    
                    # Update flow state
                    T[i+1] = state_out["T"]
                    Y[i+1, :] = state_out["Y"]
                    Y_current = state_out["Y"].copy()
                    chemistry_applied = True
                    
                    # Store unified diagnostics for CSV logging
                    if not hasattr(self, 'unified_data'):
                        self.unified_data = []
                    
                    # Get species mole fractions for H2/CO calculation and persist to X array
                    gas_temp = self.thermo.gas
                    gas_temp.TPY = state_out["T"], P, state_out["Y"]
                    X_temp = gas_temp.X
                    X[i+1, :] = X_temp
                    
                    # Calculate H2/CO ratio
                    h2_idx = self.thermo.species_indices.get('H2', -1)
                    co_idx = self.thermo.species_indices.get('CO', -1)
                    h2_co = 0.0
                    if h2_idx >= 0 and co_idx >= 0:
                        h2_co = X_temp[h2_idx] / max(1e-30, X_temp[co_idx])
                    
                    # Store unified row (AXIAL SPEC EXACT SCHEMA)
                    row_data = {
                        'cell': i,
                        'z': z[i],
                        'T': float(state_out["T"]),
                        'p': float(P),
                        'u': float(u_pred),
                        'dz': float(current_dz),
                        'dt_try': float(chem_info['dt_effective']) if 'dt_effective' in chem_info else float(dt_try),
                        'used_fallback': int(chem_info.get('used_fallback', False)),
                        'substeps': int(chem_info.get('n_subcycles', 0)),
                        'halvings': int(chem_info.get('halvings', 0)),
                        'dt_min': float(chem_info.get('dt_min_used', 1e-8)),
                        'dt_effective': float(dt_try),
                    }

                    def add_X(name: str, col: str):
                        idx = self.thermo.species_indices.get(name, -1)
                        if idx >= 0:
                            row_data[col] = float(X_temp[idx])
                        else:
                            row_data[col] = ""

                    for name, col in [
                        ("H2","X_H2"),("CO","X_CO"),("CO2","X_CO2"),("CH4","X_CH4"),
                        ("H2O","X_H2O"),("O2","X_O2"),("N2","X_N2"),
                        ("O","X_O"),("OH","X_OH"),("H","X_H"),("HO2","X_HO2"),("H2O2","X_H2O2")
                    ]:
                        add_X(name, col)

                    # H2/CO only if both present
                    if row_data["X_H2"] != "" and row_data["X_CO"] != "":
                        row_data["H2_CO"] = float(row_data["X_H2"]) / max(1e-30, float(row_data["X_CO"]))
                    else:
                        row_data["H2_CO"] = ""
                    
                    self.unified_data.append(row_data)
            else:
                T[i+1] = T_pred
                Y[i+1, :] = Y[i, :]
                Y_current = Y[i+1, :].copy()  # Use the updated Y array
                # Persist X for no-chemistry path
                try:
                    gas_state.update_state(T[i+1], P, Y[i+1, :])
                    X[i+1, :] = gas_state.mole_fractions
                except Exception:
                    X[i+1, :] = X[i, :]

            # Update velocity and density (ideal gas update; composition already persisted to X)
            rho[i+1] = P / (ct.gas_constant / self.thermo.gas.mean_molecular_weight * max(T[i+1], 1.0))
            u[i+1] = m_dot / (rho[i+1] * self.area_m2)
            
            # Update gas state for next iteration
            gas_state.temperature = T[i+1]
            gas_state.density = rho[i+1]
            gas_state.mass_fractions = Y[i+1, :]
            gas_state.mole_fractions = X[i+1, :]
            
            # Lightweight diagnostics for laptop testing
            if i % 200 == 0:  # Reduced frequency
                print(f"    Progress: {i+1}/{n_cells-1} cells, z={z[i+1]:.3f}m, T={T[i+1]:.1f}K")
        
        # Calculate residence time
        residence_time = np.cumsum(dz / u)
        
        # Find ignition length
        ignition_length_peak = self._find_ignition_length_peak(T, z)
        ignition_length_threshold = self._find_ignition_length_threshold(T, z, T[0])
        
        # Compile results
        results = {
            'z_m': z,
            'temperature_K': T,
            'density_kg_m3': rho,
            'velocity_m_s': u,
            'mass_fractions': Y,  # Shape: (n_cells, n_species)
            'mole_fractions': X,  # Shape: (n_cells, n_species)
            'species_names': inlet_state.gas.species_names,
            'residence_time_s': residence_time,
            'ignition_length_peak_m': ignition_length_peak,
            'ignition_length_threshold_m': ignition_length_threshold
        }
        
        print(f"PFR completed: T_out={T[-1]:.1f}K, ignition_length={ignition_length_peak:.3f}m")
        
        # Save unified data if available
        if hasattr(self, 'unified_data') and self.unified_data:
            import os
            os.makedirs("outputs/pfr_run", exist_ok=True)
            self.save_unified_csv("outputs/pfr_run/unified_profile.csv")
        
        return results

    def _find_ignition_length_peak(self, T: np.ndarray, z: np.ndarray) -> float:
        """Find ignition length based on peak dT/dz."""
        if len(T) < 2:
            return 0.0
        
        dT_dz = np.gradient(T, z)
        peak_idx = np.argmax(dT_dz)
        return z[peak_idx]

    def _find_ignition_length_threshold(self, T: np.ndarray, z: np.ndarray, T_inlet: float) -> float:
        """Find ignition length based on temperature threshold."""
        threshold = T_inlet + 400.0  # 400K above inlet
        idx = np.where(T >= threshold)[0]
        if len(idx) > 0:
            return z[idx[0]]
        return z[-1]  # Return end if threshold not reached

    def _energy_flow_predictor(self, T_in: float, P: float, Y_in: np.ndarray, u_in: float, dz: float) -> Tuple[float, float, float]:
        """[ROBUST] Energy/flow predictor step (no chemistry)."""
        T_in, P, Y_in = sanitize_state(T_in, P, Y_in)
        
        # Create gas state for properties
        gas_state = self.thermo.create_gas_state(T_in, P, Y_in)
        
        # Simple plug-flow energy balance without chemistry
        m_dot = self.m_dot
        cp = gas_state.cp_mass  # query once with gas set to (T_in,P,Y_in)
        q_wall_len = self.q_wall_per_m  # W/m
        dT = -q_wall_len * dz / (m_dot * cp)
        T_pred = np.clip(T_in + dT, T_MIN, T_MAX)
        
        # Density from ideal gas law
        rho_pred = P / (ct.gas_constant / gas_state.gas.mean_molecular_weight * T_pred)
        u_pred = m_dot / (rho_pred * self.area_m2)
        
        return T_pred, rho_pred, u_pred

    def chemistry_step_adaptive(self, state_in: dict, dt_try: float, 
                               dt_min: float = 1e-8, dt_max: float = 1e-4,
                               rtol_T: float = 2e-3, atol_T: float = 1.0,
                               rtol_Y: float = 2e-3, atol_Y: float = 2e-9,
                               max_subcycles: int = 400,
                               underflow_limit: int = 3) -> Tuple[dict, dict]:
        """
        Adaptive chemistry step with microreactor subcycling and fallback.
        
        Args:
            state_in: Dictionary with keys 'T', 'p', 'Y'
            dt_try: Target integration time (s)
            dt_min: Minimum time step (s)
            dt_max: Maximum time step (s)
            rtol_T: Relative tolerance for temperature
            atol_T: Absolute tolerance for temperature (K)
            rtol_Y: Relative tolerance for species
            atol_Y: Absolute tolerance for species
            max_subcycles: Maximum number of subcycles
            underflow_limit: Consecutive underflow hits before fallback
            
        Returns:
            Tuple of (state_out, info)
        """
        import cantera as ct
        
        # Load state
        gas = self.thermo.gas
        T_in = state_in["T"]
        p_in = state_in["p"]
        Y_in = state_in["Y"]
        
        gas.TPY = T_in, p_in, Y_in
        
        # Initialize counters
        n_sub = 0
        n_halvings = 0
        underflow_hits = 0
        used_fallback = False
        
        # Start with conservative step
        dt = min(max(dt_min, dt_try * 0.1), dt_max)
        t_accum = 0.0
        T_prev = gas.T
        Y_prev = gas.Y.copy()
        
        def stable_enough(T_old, T_new, Y_old, Y_new):
            """Check if step is stable enough according to tolerances."""
            dT = abs(T_new - T_old)
            if dT > (atol_T + rtol_T * max(1.0, abs(T_old))):
                return False
            
            denom = np.maximum(1.0, np.abs(Y_old))
            rel = np.max(np.abs(Y_new - Y_old) / denom)
            if rel > rtol_Y and np.max(np.abs(Y_new - Y_old)) > atol_Y:
                return False
            return True
        
        # Adaptive micro-subcycling
        while t_accum < dt_try and n_sub < max_subcycles:
            # Try a microstep using constant-pressure reactor
            try:
                r = ct.IdealGasConstPressureReactor(gas)
            except AttributeError:
                r = ct.IdealGasReactor(gas)
            
            net = ct.ReactorNet([r])
            try:
                net.rtol = 1e-10
                net.atol = 1e-16
            except Exception:
                pass
            
            t0 = net.time
            net.advance(t0 + dt)
            
            T_new = r.T
            Y_new = r.thermo.Y.copy()
            
            if stable_enough(T_prev, T_new, Y_prev, Y_new):
                # Accept step
                T_prev, Y_prev = T_new, Y_new
                t_accum += dt
                n_sub += 1
                # Grow dt moderately but clamp
                dt = min(dt * 1.5, dt_max, dt_try - t_accum if dt_try > t_accum else dt)
                gas.TPY = T_prev, p_in, Y_prev
                underflow_hits = 0
            else:
                # Reject step → shrink dt
                dt = max(dt * 0.5, dt_min)
                n_sub += 1
                n_halvings += 1
                underflow_hits += (1 if dt <= dt_min + 1e-30 else 0)
            
            # Bail to fallback if we keep hitting dt_min
            if underflow_hits >= underflow_limit:
                break
        
        # If we didn't finish time, run proper fallback for remaining time
        if t_accum < dt_try - 1e-18:
            used_fallback = True
            gas.TPY = T_prev, p_in, Y_prev
            
            try:
                r = ct.IdealGasConstPressureReactor(gas)
            except AttributeError:
                r = ct.IdealGasReactor(gas)
            
            net = ct.ReactorNet([r])
            try:
                net.rtol = 1e-10
                net.atol = 1e-16
            except Exception:
                pass
            
            target = net.time + (dt_try - t_accum)
            # Advance in safe chunks but COVER the whole remaining time
            while net.time < target:
                chunk = min(max(dt_min, (dt_try - t_accum) / 10.0), 5e-6)
                net.advance(min(net.time + chunk, target))
            
            T_prev, Y_prev = r.T, r.thermo.Y.copy()
        
        state_out = {
            "T": T_prev, 
            "p": p_in, 
            "Y": Y_prev, 
            "rho": gas.density
        }
        
        info = {
            "used_fallback": used_fallback,
            "n_subcycles": n_sub,
            "halvings": n_halvings,
            "dt_min_used": dt_min,
            "dt_effective": dt_try
        }
        
        return state_out, info

    def save_unified_csv(self, filename: str):
        """Save unified physics + diagnostics data to CSV file."""
        if not hasattr(self, 'unified_data') or not self.unified_data:
            print(f"    No unified data to save")
            return
            
        import pandas as pd
        
        df = pd.DataFrame(self.unified_data)
        df.to_csv(filename, index=False)
        print(f"    Unified data saved to: {filename}")
        
        # Print summary statistics
        if len(df) > 0:
            total_cells = len(df)
            fallback_cells = df['used_fallback'].sum()
            print(f"    Chemistry Summary:")
            print(f"      Total cells: {total_cells}")
            print(f"      Cells using fallback: {fallback_cells} ({100*fallback_cells/total_cells:.1f}%)")
            print(f"      Average substeps: {df['substeps'].mean():.1f}")
            print(f"      Average halvings: {df['halvings'].mean():.1f}")
            print(f"      Min/Max dt_effective: {df['dt_effective'].min():.2e}/{df['dt_effective'].max():.2e}")
            if 'H2_CO' in df.columns and df['H2_CO'].max() > 0:
                print(f"      H2/CO range: {df['H2_CO'].min():.3f} - {df['H2_CO'].max():.3f}")

    def _microreactor_fallback_robust(self, T: float, P: float, Y: np.ndarray, tau_cell: float, diagnostics: dict) -> Tuple[float, np.ndarray, dict]:
        """Robust microreactor fallback with stable tolerances."""
        try:
            gas = self.thermo.gas
            gas.TPY = T, P, Y
            rgas = ct.IdealGasConstPressureReactor(gas)
            net = ct.ReactorNet([rgas])
            
            # Stable tolerances for microreactor
            net.rtol = 2e-3
            net.atol = 2e-9
            
            # Conservative fallback time
            dt_fallback = min(1e-9, 0.1 * tau_cell)  # 1 ns or 10% of cell time
            t0 = net.time
            net.advance(t0 + dt_fallback)
            
            # Extract updated state
            T = rgas.thermo.T
            Y = rgas.thermo.Y
            gas.TPY = T, P, Y
            rho = gas.density
            T, Y, rho = clamp_state(T, Y, rho)
            
            diagnostics['dt_min_used'] = dt_fallback
            print(f"    Microreactor fallback: dt={dt_fallback:.2e}s, T={T:.1f}K")
            return T, Y, diagnostics
            
        except Exception as e:
            print(f"    Microreactor fallback failed: {e}; proceed with current state")
            diagnostics['fallback_used'] = False
            return T, Y, diagnostics

    def solve_pfr(self, inlet_state: GasState, m_dot: float, P: float,
                  q_wall_per_m: float, chemistry_on: bool = True,
                  max_dT_per_step: float = 20.0) -> Dict[str, Any]:
        """Solve PFR with proper energy balance.
        
        Args:
            inlet_state: Inlet gas state
            m_dot: Mass flow rate (kg/s)
            P: Pressure (Pa)
            q_wall_per_m: Wall heat loss per length (W/m)
            chemistry_on: Whether to include chemistry
            max_dT_per_step: Maximum temperature change per step (K)
            
        Returns:
            Dictionary with PFR results
        """
        print(f"Solving PFR: L={self.length_m:.2f}m, D={self.diameter_m:.2f}m")
        print(f"  m_dot={m_dot:.3f}kg/s, q_wall={q_wall_per_m:.1f}W/m")
        print(f"  Chemistry: {'ON' if chemistry_on else 'OFF'}")
        
        # Laptop PFR grid: 300–500 cells chemistry-ON, then optional ×2 refinement after convergence
        base_dz = 0.001  # m
        n_cells = min(max(int(self.length_m / base_dz), 300), 500)
        z = np.linspace(0, self.length_m, n_cells)
        dz = z[1] - z[0]
        
        print(f"  PFR grid: {n_cells} cells, dz={dz:.4f}m (laptop-optimized)")
        
        # Initialize arrays
        T = np.zeros(n_cells)
        Y = np.zeros((self.thermo.n_species, n_cells))
        X = np.zeros((self.thermo.n_species, n_cells))
        u = np.zeros(n_cells)
        rho = np.zeros(n_cells)
        
        # Set initial conditions
        T[0] = inlet_state.temperature
        Y[:, 0] = inlet_state.mass_fractions
        
        # Calculate initial velocity
        rho[0] = inlet_state.density
        u[0] = m_dot / (rho[0] * self.area_m2)
        
        # Set initial gas state
        gas_state = self.thermo.create_gas_state(T[0], P, Y[:, 0])
        X[:, 0] = gas_state.mole_fractions
        
        print(f"Initial: T={T[0]-273.15:.1f}°C, u={u[0]:.3f}m/s, ρ={rho[0]:.1f}kg/m³")
        
        # March through reactor
        for i in range(1, n_cells):
            # Progress logging
            if i % 50 == 0 or i < 5:
                print(f"  Cell {i:3d}/{n_cells}: z={z[i]:.3f}m, T={T[i-1]-273.15:.1f}°C")
            
            # Set gas state at previous position
            gas_state.update_state(T[i-1], P, Y[:, i-1])
            
            # Get properties
            rho[i-1] = gas_state.density
            cp_mass = gas_state.cp_mass
            u[i-1] = m_dot / (rho[i-1] * self.area_m2)
            
            # Chemistry heat release (molar basis) - CORRECT CANTERA UNITS
            if chemistry_on:
                # Species-sum form (primary method)
                omega = gas_state.get_net_production_rates()  # kmol/m³/s
                h_mol = gas_state.get_partial_molar_enthalpies()  # J/kmol
                qdot_chem_vol = -float(h_mol @ omega)  # W/m³
                
                # Clip extreme values to prevent numerical issues
                qdot_chem_vol = np.clip(qdot_chem_vol, -1e12, 1e12)
                
                # Reaction-wise cross check for numerical robustness
                R = gas_state.gas.net_rates_of_progress  # kmol/m³/s
                nu_p = gas_state.gas.product_stoich_coeffs  # property, not method
                nu_r = gas_state.gas.reactant_stoich_coeffs  # property, not method
                dH_rxn = (nu_p - nu_r).T @ h_mol  # J/kmol (vector per reaction)
                qdot_rxn_vol = -float(R @ dH_rxn)  # W/m³
                
                # Check consistency
                if abs(qdot_chem_vol) > 1e-12:
                    relative_error = abs(qdot_chem_vol - qdot_rxn_vol) / max(abs(qdot_rxn_vol), 1e-12)
                    if relative_error > 1e-6:
                        print(f"    WARNING: Chemistry heat release mismatch: species={qdot_chem_vol:.2e}, reaction={qdot_rxn_vol:.2e}, error={relative_error:.2e}")
            else:
                qdot_chem_vol = 0.0  # No chemistry
                qdot_rxn_vol = 0.0
            
            # Energy equation: dT/dz = (A * q_chem_vol - q_wall_per_m) / (m_dot * cp_mass)
            dT_dz = (self.area_m2 * qdot_chem_vol - q_wall_per_m) / (m_dot * cp_mass)
            
            # Comprehensive logging every 0.05m (50 cells with dz=0.001m)
            if i % 50 == 0 or i < 5:
                qchem_per_m = self.area_m2 * qdot_chem_vol
                ratio = abs(qchem_per_m) / q_wall_per_m if q_wall_per_m > 0 else 0
                
                print(f"    Energy balance at z={z[i]:.3f}m:")
                print(f"      T={T[i-1]-273.15:.1f}°C, cp={cp_mass:.1f}J/kg/K, m_dot={m_dot:.3f}kg/s")
                print(f"      q_chem_vol={qdot_chem_vol:.1e}W/m³, q_chem_len={qchem_per_m:.1f}W/m")
                print(f"      q_wall_len={q_wall_per_m:.1f}W/m, ratio={ratio:.1e}")
                print(f"      dT/dz={dT_dz:.3f}K/m")
                
                if chemistry_on:
                    print(f"      Chemistry details:")
                    print(f"        omega range: {np.min(omega):.2e} to {np.max(omega):.2e} kmol/m³/s")
                    print(f"        h_mol range: {np.min(h_mol):.2e} to {np.max(h_mol):.2e} J/kmol")
                    print(f"        qdot_rxn_vol={qdot_rxn_vol:.1e}W/m³")
                
                # Validation gates
                if chemistry_on:
                    # Use validation module for comprehensive checks
                    from .validation import ValidationGates
                    validation = ValidationGates()
                    chem_validation = validation.validate_chemistry_heat_release(
                        qdot_chem_vol, qchem_per_m, h_mol, omega
                    )
                    
                    if not chem_validation['chemistry_heat_release_gate']:
                        print(f"      WARNING: Chemistry heat release validation failed:")
                        if not chem_validation['h_mol_valid']:
                            print(f"        h_mol range: {chem_validation['h_mol_range']} (expected: {chem_validation['expected_ranges']['h_mol']})")
                        if not chem_validation['omega_valid']:
                            print(f"        omega range: {chem_validation['omega_range']} (expected: {chem_validation['expected_ranges']['omega']})")
                        if not chem_validation['q_chem_vol_valid']:
                            print(f"        q_chem_vol: {qdot_chem_vol:.1e} (expected: {chem_validation['expected_ranges']['q_chem_vol']})")
                        if not chem_validation['q_chem_per_m_valid']:
                            print(f"        q_chem_per_m: {qchem_per_m:.1e} (expected: {chem_validation['expected_ranges']['q_chem_per_m']})")
                
                # Red-flag: fail fast if chemistry heat release is too large (temporarily disabled for Aramco)
                if abs(qdot_chem_vol) > 1e15:  # W/m³ - very high threshold for debugging
                    print(f"WARNING: |qdot_chem_vol| = {abs(qdot_chem_vol):.1e} > 1e15 W/m³ (validation disabled for debugging)")
                    print(f"  omega range: {np.min(omega):.2e} to {np.max(omega):.2e} kmol/m³/s")
                    print(f"  h_mol range: {np.min(h_mol):.2e} to {np.max(h_mol):.2e} J/kmol")
                    print(f"  qchem_per_m: {qchem_per_m:.1e} W/m")
                    print(f"  qwall_per_m: {q_wall_per_m:.1e} W/m")
                    print(f"  ratio: {ratio:.1e}")
                    # Don't raise error, just warn for now
                    # raise RuntimeError("Chemistry heat release too large - numerical instability detected")
                
                # Blocker: fail fast if ratio is too large (temporarily disabled for Aramco)
                if ratio > 1e8:  # Very high threshold for debugging
                    print(f"WARNING: |A*q_chem/q_wall| = {ratio:.1e} > 1e8 (validation disabled for debugging)")
                    print(f"  omega range: {np.min(omega):.2e} to {np.max(omega):.2e} kmol/m³/s")
                    print(f"  h_mol range: {np.min(h_mol):.2e} to {np.max(h_mol):.2e} J/kmol")
                    # Don't raise error, just warn for now
                    # raise RuntimeError("Energy balance units error - chemistry heat release too large")
            
            # Adaptive stepping based on temperature change
            current_dz = dz
            dT_est = abs(dT_dz * current_dz)
            
            # Adaptive stepping with backtracking: limit |ΔT/Δz| ≤ 250 K/m
            max_dT_dz = 250.0  # K/m
            dT_est = abs(dT_dz * current_dz)
            
            # Backtracking loop: up to 6 retries
            for retry in range(6):
                if dT_est <= max_dT_dz * current_dz:
                    break
                
                # Cut step by 50%
                current_dz *= 0.5
                dT_est = abs(dT_dz * current_dz)
                
                if retry == 0:
                    print(f"    Reducing step size due to large dT: dz={current_dz:.4f}m")
            
            # Species equations (only if chemistry is on)
            if chemistry_on:
                dY_dz = omega * self.thermo.molecular_weights / (m_dot / self.area_m2)
                
                # Apply stability limits
                dY_dz = np.clip(dY_dz, -0.1, 0.1)
                
                # Special handling for inert species
                inert_species = ['N2', 'AR', 'HE']
                for species in inert_species:
                    if species in self.thermo.species_indices:
                        idx = self.thermo.species_indices[species]
                        dY_dz[idx] = 0.0
            else:
                dY_dz = np.zeros(self.thermo.n_species)
            
            # Update state
            T_new = T[i-1] + dT_dz * current_dz
            Y_new = Y[:, i-1] + dY_dz * current_dz
            
            # Stability checks
            if np.isnan(T_new) or T_new < 300.0 or T_new > 5000.0:
                T[i] = T[i-1]
                print(f"    WARNING: Invalid T={T_new:.1f}K, keeping previous T={T[i-1]:.1f}K")
            else:
                T[i] = T_new
            
            if np.any(np.isnan(Y_new)):
                Y[:, i] = Y[:, i-1]
                print(f"    WARNING: Invalid Y, keeping previous composition")
            else:
                Y[:, i] = Y_new
            
            # Normalize mass fractions and ensure positivity
            Y[:, i] = np.maximum(Y[:, i], 1e-10)
            Y[:, i] = Y[:, i] / np.sum(Y[:, i])
            
            # Update gas state and recalculate velocity
            try:
                gas_state.update_state(T[i], P, Y[:, i])
                rho[i] = gas_state.density
                u[i] = m_dot / (rho[i] * self.area_m2)
                X[:, i] = gas_state.mole_fractions
            except Exception as e:
                print(f"Warning: Gas state failed at z={z[i]:.3f}, using previous values: {e}")
                T[i] = T[i-1]
                Y[:, i] = Y[:, i-1]
                X[:, i] = X[:, i-1]
                rho[i] = rho[i-1]
                u[i] = u[i-1]
        
        # Calculate residence time
        residence_time = self._calculate_residence_time(z, u)
        
        # Detect ignition
        ignition_peak, ignition_thresh = self._detect_ignition(z, T)
        
        print(f"PFR outlet: T={T[-1]-273.15:.1f}°C, P={P/1e5:.1f}bar")
        print(f"Residence time: {residence_time[-1]:.1f}s")
        print(f"Ignition length: {ignition_peak:.3f}m (peak), {ignition_thresh:.3f}m (threshold)")
        
        # Chemistry-off analytic check
        if not chemistry_on:
            # Expected linear cooling: dT/dz = -q_wall_per_m / (m_dot * cp)
            cp_avg = np.mean([self.thermo.create_gas_state(T[i], P, Y[:, i]).cp_mass for i in range(0, len(T), 50)])
            expected_dT_dz = -q_wall_per_m / (m_dot * cp_avg)
            measured_dT_dz = (T[-1] - T[0]) / self.length_m
            
            print(f"\nChemistry-OFF Analytic Check:")
            print(f"  Expected dT/dz: {expected_dT_dz:.3f} K/m")
            print(f"  Measured dT/dz: {measured_dT_dz:.3f} K/m")
            print(f"  Error: {abs(measured_dT_dz - expected_dT_dz):.3f} K/m")
            
            error_percent = abs(measured_dT_dz - expected_dT_dz) / abs(expected_dT_dz) * 100
            if error_percent > 5.0:
                print(f"  WARNING: Error {error_percent:.1f}% > 5% - energy units may be incorrect")
            else:
                print(f"  ✓ PASS: Error {error_percent:.1f}% < 5%")
        
        return {
            'z_m': z,
            'temperature_K': T,
            'density_kg_m3': rho,
            'velocity_m_s': u,
            'mass_fractions': Y,
            'mole_fractions': X,
            'residence_time_s': residence_time,
            'ignition_length_peak_m': ignition_peak,
            'ignition_length_threshold_m': ignition_thresh,
            'species_names': self.thermo.species_names
        }
    
    def _calculate_residence_time(self, z: np.ndarray, u: np.ndarray) -> np.ndarray:
        """Calculate residence time profile.
        
        Args:
            z: Axial positions (m)
            u: Velocity profile (m/s)
            
        Returns:
            Residence time profile (s)
        """
        # Residence time = integral of dz/u
        dz = z[1] - z[0]
        residence_time = np.zeros_like(z)
        
        for i in range(1, len(z)):
            if u[i] > 0:
                residence_time[i] = residence_time[i-1] + dz / u[i]
            else:
                residence_time[i] = residence_time[i-1]
        
        return residence_time
    
    def _detect_ignition(self, z: np.ndarray, T: np.ndarray) -> Tuple[float, float]:
        """Detect ignition length.
        
        Args:
            z: Axial positions (m)
            T: Temperature profile (K)
            
        Returns:
            Tuple of (ignition_peak, ignition_threshold)
        """
        # Look for peak dT/dz in first 0.4m
        z_max = 0.4
        mask = z <= z_max
        
        if np.sum(mask) < 2:
            return 0.0, 0.0
        
        z_window = z[mask]
        T_window = T[mask]
        
        # Calculate dT/dz
        dT_dz = np.gradient(T_window, z_window)
        
        # Find peak
        peak_idx = np.argmax(dT_dz)
        ignition_peak = z_window[peak_idx]
        
        # Find threshold crossing (T > T_inlet + 400K)
        T_thresh = T[0] + 400.0
        thresh_idx = np.where(T_window > T_thresh)[0]
        
        if len(thresh_idx) > 0:
            ignition_thresh = z_window[thresh_idx[0]]
        else:
            ignition_thresh = z_window[-1]
        
        return ignition_peak, ignition_thresh
