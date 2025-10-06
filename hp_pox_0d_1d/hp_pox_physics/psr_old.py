"""
Perfectly Stirred Reactor (PSR) with fixed burner heat loss.
Uses fixed Q_burner values from Richter benchmark with robust numerical methods.
"""

import cantera as ct
import numpy as np
import time
from scipy.optimize import newton
from typing import Tuple, Optional, Dict, Any
from .thermo import ThermoManager, GasState


class Watchdog:
    """Watchdog for monitoring solver progress and preventing infinite loops."""
    
    def __init__(self, wall_time_limit: float = 180.0, max_steps: int = 1000):
        """Initialize watchdog.
        
        Args:
            wall_time_limit: Maximum wall time in seconds
            max_steps: Maximum number of steps before abort
        """
        self.wall_start = time.time()
        self.wall_time_limit = wall_time_limit
        self.max_steps = max_steps
        self.step_count = 0
        self.last_T = 0.0
        self.last_resid = float('inf')
        self.stagnation_hits = 0
        self.max_stagnation = 5
        
    def should_abort(self, msg: str = "") -> bool:
        """Check if solver should abort.
        
        Args:
            msg: Additional message for abort
            
        Returns:
            True if should abort, False otherwise
        """
        wall_time = time.time() - self.wall_start
        
        if wall_time > self.wall_time_limit:
            raise RuntimeError(f"WATCHDOG ABORT: Wall time limit exceeded ({wall_time:.1f}s > {self.wall_time_limit:.1f}s). {msg}")
        
        if self.step_count > self.max_steps:
            raise RuntimeError(f"WATCHDOG ABORT: Step limit exceeded ({self.step_count} > {self.max_steps}). {msg}")
        
        if self.stagnation_hits > self.max_stagnation:
            raise RuntimeError(f"WATCHDOG ABORT: Too many stagnation hits ({self.stagnation_hits} > {self.max_stagnation}). {msg}")
        
        return False
    
    def log_progress(self, T: float, residual: float, max_omega: float = 0.0):
        """Log solver progress.
        
        Args:
            T: Current temperature
            residual: Current residual
            max_omega: Maximum production rate magnitude
        """
        self.step_count += 1
        wall_time = time.time() - self.wall_start
        
        # Check for stagnation
        if abs(T - self.last_T) < 1.0 and abs(residual - self.last_resid) < 1e-6:
            self.stagnation_hits += 1
        else:
            self.stagnation_hits = 0
        
        self.last_T = T
        self.last_resid = residual
        
        # Print progress every 10 steps
        if self.step_count % 10 == 0:
            print(f"    Step {self.step_count}: T={T:.1f}K, residual={residual:.2e}, max|ω|={max_omega:.2e}, time={wall_time:.1f}s")
            import sys
            sys.stdout.flush()


class PSR:
    """Perfectly Stirred Reactor with constant pressure operation."""
    
    def __init__(self, volume_m3: float, thermo: ThermoManager):
        """Initialize PSR.
        
        Args:
            volume_m3: Reactor volume (m³)
            thermo: Thermodynamics manager
        """
        self.volume_m3 = volume_m3
        self.thermo = thermo
    
    def solve_psr_continuation(self, inlet_state: GasState, m_dot: float, 
                              P: float, q_burner_kw: float) -> Tuple[GasState, Dict[str, Any], bool]:
        """
        Solve PSR to steady state with fixed Q_burner heat loss using ignition-friendly continuation.
        
        Args:
            inlet_state: Inlet gas state
            m_dot: Mass flow rate (kg/s)
            P: Pressure (Pa)
            q_burner_kw: Fixed burner heat loss (kW, positive for cooling)
        
        Returns:
            Tuple of (outlet_state, diagnostics, converged)
        """
        print(f"  Solving PSR with Q_burner = {q_burner_kw:.2f} kW using ignition-friendly continuation")
        
        # Convert kW to W
        q_burner_w = q_burner_kw * 1000.0
        
        # Calculate residence time
        res_time_s = self.volume_m3 * inlet_state.density / m_dot
        dt_phys = max(0.2 * res_time_s, 0.001)  # Physical time step with minimum
        
        # A) Ignition-friendly seeding
        T_seed = max(1100.0, inlet_state.temperature + 250.0)
        T_seed = min(T_seed, 1400.0)  # Cap at 1400K
        print(f"  Ignition seeding: T_seed = {T_seed:.1f}K")
        
        # Create PSR reactor for seeding
        reactor = ct.IdealGasConstPressureReactor(self.thermo.gas)
        reactor.volume = self.volume_m3
        reactor.thermo.TPY = T_seed, P, inlet_state.mass_fractions
        
        # Auto-ignition phase with Q=0
        print(f"  Auto-ignition phase: Q=0, up to 10*τ = {10*res_time_s:.3f}s")
        sim = ct.ReactorNet([reactor])
        sim.rtol = 1e-8
        sim.atol = 1e-14
        sim.max_steps = int(2e6)
        
        # Set up vector atols
        n_species = len(reactor.thermo.species_names)
        atol_vector = np.full(n_species + 1, 1e-14)  # +1 for temperature
        atol_vector[0] = 1e-9  # Temperature atol
        
        # Set stricter atols for radicals
        radical_names = ['H', 'O', 'OH', 'HO2', 'CH3', 'CH2', 'CH', 'C2H5', 'C2H3', 'C2H', 'C3H7', 'C3H5', 'C3H3']
        for i, name in enumerate(reactor.thermo.species_names):
            if name in radical_names:
                atol_vector[i + 1] = 1e-18  # +1 because T is first
        
        try:
            sim.atol = atol_vector
        except:
            pass  # Use scalar atol if vector not supported
        
        # Run auto-ignition phase
        max_omega = 0.0
        max_qdot_chem = 0.0
        time_to_reactive = 0.0
        attempts = 0
        
        while attempts < 6:  # Max 6 attempts
            attempts += 1
            print(f"    Attempt {attempts}: T_seed = {T_seed:.1f}K")
            
            # Reset reactor state
            reactor.thermo.TPY = T_seed, P, inlet_state.mass_fractions
            sim.reinitialize()
            
            # Run auto-ignition
            t_max = 10 * res_time_s
            dt_step = 0.1 * res_time_s
            
            for t in np.arange(0, t_max, dt_step):
                try:
                    sim.advance(sim.time + dt_step)  # Advance from current time
                    
                    # Check reactivity
                    omega = reactor.thermo.net_production_rates
                    h_mol = reactor.thermo.partial_molar_enthalpies
                    qdot_chem_vol = -float(h_mol @ omega)
                    
                    max_omega = max(max_omega, float(np.max(np.abs(omega))))
                    max_qdot_chem = max(max_qdot_chem, abs(qdot_chem_vol))
                    
                    if max_omega > 1e-3 or abs(qdot_chem_vol) > 1e6:
                        time_to_reactive = t
                        print(f"    Reactive at t={t:.3f}s: max|ω|={max_omega:.2e}, |qdot|={abs(qdot_chem_vol):.2e}")
                        break
                        
                except Exception as e:
                    print(f"    Auto-ignition failed at t={t:.3f}s: {e}")
                    break
            
            # Check if reactive enough
            if max_omega > 1e-4 and max_qdot_chem > 1e5:
                print(f"    Ignition successful: max|ω|={max_omega:.2e}, |qdot|={max_qdot_chem:.2e}")
                break
            else:
                T_seed += 50.0  # Bump temperature
                if T_seed > 1400.0:
                    print(f"    Warning: Max T_seed reached, proceeding with current state")
                    break
        else:
            print(f"    Warning: Auto-ignition phase did not achieve reactivity")
        
        # B) Q_burner ramping with geometric increments
        print(f"  Q_burner ramping: 5% → 100% of target")
        q_increments = [0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00]
        q_burner_increments = [q * q_burner_w for q in q_increments]
        
        # Create wall for heat loss
        env = ct.Reservoir(self.thermo.gas)
        wall = ct.Wall(reactor, env)
        
        print(f"  Residence time: {res_time_s:.3f} s, dt_phys: {dt_phys:.3f} s")
        print(f"  Q_burner increments: {[f'{q/1000:.1f}' for q in q_burner_increments]} kW")
        
        # Initialize at Q=0
        reactor.thermo.TPY = 1100.0, P, inlet_state.mass_fractions
        wall.heat_transfer_coeff = 0.0
        
        converged = True
        final_residual = float('inf')
        
        # Initialize diagnostics
        diagnostics = {
            'T_out': 0.0,
            'Q_burner': q_burner_w,
            'Q_chem': 0.0,
            'residual': float('inf'),
            'converged': False,
            'steps': 0,
            'omega_max': 0.0,
            'omega_min': 0.0,
            'h_mol_max': 0.0,
            'h_mol_min': 0.0,
            'h_in_J_kg': inlet_state.gas.enthalpy_mass,
            'h_out_J_kg': 0.0,
            'enthalpy_difference_J_kg': 0.0,
            'n_increments': len(q_burner_increments),
            'final_heat_transfer_coeff': 0.0,
            'res_time_s': res_time_s,
            'dt_phys': dt_phys,
            'T_seed': T_seed,
            'ignition_attempts': attempts,
            'time_to_reactive': time_to_reactive,
            'max_omega_seeding': max_omega,
            'max_qdot_chem_seeding': max_qdot_chem
        }
        
        # Ramp Q_burner to target value
        for i, q_current in enumerate(q_burner_increments):
            print(f"  Increment {i+1}/{len(q_burner_increments)}: Q_burner = {q_current/1000:.2f} kW")
            
            # Set heat transfer coefficient (negative for cooling)
            # Use a more stable calculation for heat transfer coefficient
            T_diff = max(reactor.thermo.T - 300.0, 100.0)  # Minimum temperature difference
            wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
            
            # Bounded loop: advance in physical time steps
            residual_history = []
            best_residual = float('inf')
            stagnation_count = 0
            
            for k in range(50):  # Max 50 steps per increment
                # Calculate current residual
                h_in = inlet_state.gas.enthalpy_mass
                h_out = reactor.thermo.enthalpy_mass
                h_dot_in = h_in * m_dot
                h_dot_out = h_out * m_dot
                
                # Calculate chemistry heat release
                try:
                    omega = reactor.thermo.net_production_rates
                    h_mol = reactor.thermo.partial_molar_enthalpies
                    q_chem_vol = -float(h_mol @ omega)
                    q_chem_w = q_chem_vol * self.volume_m3
                except:
                    q_chem_w = 0.0
                
                residual = abs(q_chem_w + q_current + h_dot_out - h_dot_in) / max(abs(h_dot_in), 1.0)
                residual_history.append(residual)
                
                # Check for convergence
                if residual < 1e-3 and max_omega > 1e-4:  # Require both convergence and reactivity
                    print(f"    Converged at step {k+1}: residual = {residual:.2e}, max|ω| = {max_omega:.2e}")
                    break
                
                # Check for stagnation
                if len(residual_history) >= 3:
                    recent_improvement = (residual_history[-3] - residual) / residual_history[-3]
                    if recent_improvement < 0.01:  # <1% improvement
                        stagnation_count += 1
                        if stagnation_count >= 3:
                            print(f"    Stagnation detected at step {k+1}, halving increment")
                            q_current *= 0.5
                            wall.heat_transfer_coeff = -q_current / (reactor.volume * (reactor.thermo.T - 300.0))
                            sim.reinitialize()
                            stagnation_count = 0
                            continue
                    else:
                        stagnation_count = 0
                
                # Check for increasing residual
                if residual > best_residual * 1.1:  # 10% increase
                    print(f"    Residual increasing at step {k+1}, halving increment")
                    q_current *= 0.5
                    wall.heat_transfer_coeff = -q_current / (reactor.volume * (reactor.thermo.T - 300.0))
                    sim.reinitialize()
                    continue
                
                best_residual = min(best_residual, residual)
                
                # Advance in time with max horizon guard
                t_target = min(sim.time + dt_phys, sim.time + 2 * res_time_s)
                
                try:
                    sim.advance(t_target)
                    
                    # Enforce bounds
                    T = reactor.thermo.T
                    if T < 300.0 or T > 3500.0:
                        print(f"    Temperature out of bounds: T={T:.1f}K, rejecting step")
                        q_current *= 0.5
                        wall.heat_transfer_coeff = -q_current / (reactor.volume * (reactor.thermo.T - 300.0))
                        sim.reinitialize()
                        continue
                    
                    # Species floor and renormalization
                    Y = reactor.thermo.Y
                    Y = np.maximum(Y, 1e-30)
                    Y = Y / np.sum(Y)
                    reactor.thermo.Y = Y
                    
                    # Check density
                    rho = reactor.thermo.density
                    if rho <= 0.0:
                        print(f"    Invalid density: rho={rho:.2e}, rejecting step")
                        q_current *= 0.5
                        wall.heat_transfer_coeff = -q_current / (reactor.volume * (reactor.thermo.T - 300.0))
                        sim.reinitialize()
                        continue
                    
                except Exception as e:
                    print(f"    Integration failed at step {k+1}: {e}")
                    q_current *= 0.5
                    wall.heat_transfer_coeff = -q_current / (reactor.volume * (reactor.thermo.T - 300.0))
                    sim.reinitialize()
                    continue
            
            # Check if this increment converged
            if residual >= 1e-3:  # Further relaxed convergence criteria
                print(f"    Increment {i+1} did not converge: residual = {residual:.2e}")
                converged = False
                break
            
            final_residual = residual
        
        # Final state
        T_out = reactor.thermo.T
        Y_out = reactor.thermo.Y
        
        # Create outlet state
        outlet_state = self.thermo.create_gas_state(T_out, P, Y_out)
        
        # Calculate final chemistry heat release
        try:
            omega_out = reactor.thermo.net_production_rates
            h_mol_out = reactor.thermo.partial_molar_enthalpies
            q_chem_vol_out = -float(h_mol_out @ omega_out)
            q_chem_w_out = q_chem_vol_out * self.volume_m3
        except Exception as e:
            print(f"    Warning: Final chemistry calculation failed: {e}")
            q_chem_w_out = 0.0
            omega_out = np.zeros(len(self.thermo.species_names))
            h_mol_out = np.zeros(len(self.thermo.species_names))
        
            # Update diagnostics
            diagnostics.update({
                'T_out': T_out,
                'Q_chem': q_chem_w_out,
                'residual': final_residual,
                'converged': converged,
                'steps': sim.step_count if hasattr(sim, 'step_count') else 0,
                'omega_max': float(np.max(np.abs(omega_out))),
                'omega_min': float(np.min(np.abs(omega_out))),
                'h_mol_max': float(np.max(np.abs(h_mol_out))),
                'h_mol_min': float(np.min(np.abs(h_mol_out))),
                'h_out_J_kg': outlet_state.gas.enthalpy_mass,
                'enthalpy_difference_J_kg': outlet_state.gas.enthalpy_mass - inlet_state.gas.enthalpy_mass,
                'final_heat_transfer_coeff': wall.heat_transfer_coeff
            })
        
        # Print single line summary
        status = "CONVERGED" if converged else "FAILED"
        print(f"  PSR: Q_burner={q_burner_kw:.2f} kW | Q_chem={q_chem_w_out/1000:.2f} kW | Residual={final_residual:.2e} | T_out={T_out:.1f}K | {status}")
        
        return outlet_state, diagnostics, converged

    def solve_psr_with_fixed_heat_loss(self, inlet_state: GasState, m_dot: float, 
                                      P: float, q_burner_kw: float, t_end: float = 10.0) -> Tuple[GasState, float, bool, Dict[str, Any]]:
        """Solve PSR with fixed burner heat loss using CVODE BDF with continuation strategy.
        
        Args:
            inlet_state: Inlet gas state
            m_dot: Mass flow rate (kg/s)
            P: Pressure (Pa)
            q_burner_kw: Fixed burner heat loss (kW, positive for cooling)
            t_end: Integration time (s)
            
        Returns:
            Tuple of (outlet_state, Q_burner_W, converged, diagnostics)
        """
        print(f"Solving PSR with fixed Q_burner = {q_burner_kw:.2f} kW using CVODE BDF")
        
        # Convert kW to W
        q_burner_w = q_burner_kw * 1000.0
        
        # Create PSR reactor
        reactor = ct.IdealGasConstPressureReactor(self.thermo.gas)
        reactor.volume = self.volume_m3
        
        # Set inlet conditions
        reactor.thermo.TPY = inlet_state.temperature, P, inlet_state.mass_fractions
        
        # Create wall for heat loss
        env = ct.Reservoir(self.thermo.gas)
        wall = ct.Wall(reactor, env)
        
        # Continuation strategy: ramp Q_burner from 0 to target value
        n_increments = 10
        q_burner_increments = np.linspace(0.0, q_burner_w, n_increments)
        
        # Step 1: Solve at Q=0 from 1000K seed to steady state
        print("  Step 1: Solving at Q=0 from 1000K seed...")
        reactor.thermo.TPY = 1000.0, P, inlet_state.mass_fractions
        wall.heat_transfer_coeff = 0.0  # No heat loss initially
        
        sim = ct.ReactorNet([reactor])
        sim.rtol = 1e-9
        sim.atol = 1e-15
        sim.max_steps = int(2e6)
        
        # Advance to steady state with multiple attempts
        try:
            sim.advance(t_end)
            print(f"    Initial steady state: T={reactor.thermo.T:.1f}K")
        except Exception as e:
            print(f"    Warning: Initial steady state failed: {e}")
            # Try with more conservative settings
            sim.rtol = 1e-6
            sim.atol = 1e-12
            try:
                sim.advance(t_end)
                print(f"    Retry steady state: T={reactor.thermo.T:.1f}K")
            except Exception as e2:
                print(f"    Warning: Retry also failed: {e2}")
                # Continue with current state
        
        # Step 2: Ramp Q_burner to target value in increments
        print("  Step 2: Ramping Q_burner to target value...")
        for i, q_current in enumerate(q_burner_increments[1:], 1):
            # Set heat transfer coefficient (negative for cooling)
            # Use a more stable calculation for heat transfer coefficient
            T_diff = max(reactor.thermo.T - 300.0, 100.0)  # Minimum temperature difference
            wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
            
            # Advance to steady state
            try:
                sim.advance(t_end)
                print(f"    Increment {i}/{n_increments}: Q={q_current/1000:.1f}kW, T={reactor.thermo.T:.1f}K")
            except Exception as e:
                print(f"    Warning: Increment {i} failed: {e}")
                # Continue with current state
        
        # Final state
        T_out = reactor.thermo.T
        Y_out = reactor.thermo.Y
        
        # Create outlet state
        outlet_state = self.thermo.create_gas_state(T_out, P, Y_out)
        
        # Calculate final chemistry heat release
        try:
            omega_out = reactor.thermo.net_production_rates
            h_mol_out = reactor.thermo.partial_molar_enthalpies
            q_chem_vol_out = -float(h_mol_out @ omega_out)
            q_chem_w_out = q_chem_vol_out * self.volume_m3
        except Exception as e:
            print(f"    Warning: Final chemistry calculation failed: {e}")
            q_chem_w_out = 0.0
            omega_out = np.zeros(len(self.thermo.species_names))
            h_mol_out = np.zeros(len(self.thermo.species_names))
        
        # Calculate heat balance residual
        h_in = inlet_state.gas.enthalpy_mass
        h_out = outlet_state.gas.enthalpy_mass
        heat_balance_residual = abs((h_out - h_in) * m_dot + q_chem_w_out + q_burner_w) / max(abs(q_chem_w_out), abs(q_burner_w), 1.0)
        
        # Check convergence
        converged = heat_balance_residual < 1e-6 and abs(T_out - 1000.0) > 50.0  # Not stuck at initial temperature
        
        # Calculate PSR internal area for diagnostics
        aspect_ratio = 2.0
        D = (4 * self.volume_m3 / (np.pi * aspect_ratio))**(1/3)
        L = aspect_ratio * D
        A_psr = np.pi * D * L
        
        # Comprehensive diagnostics
        diagnostics = {
            'T_psr_out_K': T_out,
            'T_psr_out_C': T_out - 273.15,
            'Q_burner_W': q_burner_w,
            'Q_chem_W': q_chem_w_out,
            'heat_balance_residual': heat_balance_residual,
            'A_psr_m2': A_psr,
            'converged': converged,
            'omega_max': float(np.max(np.abs(omega_out))),
            'omega_min': float(np.min(np.abs(omega_out))),
            'h_mol_max': float(np.max(np.abs(h_mol_out))),
            'h_mol_min': float(np.min(np.abs(h_mol_out))),
            'h_in_J_kg': h_in,
            'h_out_J_kg': h_out,
            'enthalpy_difference_J_kg': h_out - h_in,
            'n_increments': n_increments,
            'final_heat_transfer_coeff': wall.heat_transfer_coeff
        }
        
        # Print single line summary
        status = "CONVERGED" if converged else "FAILED"
        print(f"PSR: Q_burner={q_burner_kw:.2f} kW | Q_chem={q_chem_w_out/1000:.2f} kW | Residual={heat_balance_residual:.2e} | T_out={T_out:.1f}K | {status}")
        
        return outlet_state, q_burner_w, converged, diagnostics
    
    def solve_psr_steady_state(self, inlet_state: GasState, m_dot: float, 
                              P: float, Q_burner: float = 0.0, t_end: float = 10.0) -> Tuple[GasState, bool]:
        """Solve PSR at steady state with given heat input.
        
        Args:
            inlet_state: Inlet gas state
            m_dot: Mass flow rate (kg/s)
            P: Pressure (Pa)
            Q_burner: Heat input rate (W, negative for cooling)
            t_end: Integration time (s)
            
        Returns:
            Tuple of (outlet_state, converged)
        """
        # Create constant pressure reactor
        reactor = ct.IdealGasConstPressureReactor(self.thermo.gas, energy='on')
        reactor.volume = self.volume_m3
        
        # Set initial state
        reactor.thermo.TPY = inlet_state.temperature, P, inlet_state.mass_fractions
        
        # Create mass flow controller
        inlet = ct.Reservoir(self.thermo.gas)
        inlet.thermo.TPY = inlet_state.temperature, P, inlet_state.mass_fractions
        
        mfc = ct.MassFlowController(inlet, reactor)
        mfc.mass_flow_rate = m_dot
        
        # Create wall for heat transfer
        wall = ct.Wall(reactor, reactor)
        wall.heat_transfer_coeff = 1.0  # W/m²/K
        wall.area = 1.0  # m²
        
        # Set heat flux
        if Q_burner != 0.0:
            wall.heat_flux = -Q_burner / wall.area  # Negative for cooling
        
        # Create reactor network
        sim = ct.ReactorNet([reactor])
        sim.rtol = 1e-8
        sim.atol = 1e-11
        
        # Integrate to steady state
        try:
            sim.advance_to_steady_state()
            converged = True
        except Exception as e:
            print(f"PSR steady state failed: {e}")
            sim.advance(t_end)
            converged = False
        
        # Get final state
        T_final = reactor.T
        Y_final = reactor.Y
        
        print(f"PSR outlet: T = {T_final-273.15:.1f}°C, Q = {Q_burner/1000:.1f}kW")
        
        # Create outlet state
        outlet_state = self.thermo.create_gas_state(T_final, P, Y_final)
        
        return outlet_state, converged
    
    def calculate_residence_time(self, m_dot: float, outlet_state: GasState) -> float:
        """Calculate residence time.
        
        Args:
            m_dot: Mass flow rate (kg/s)
            outlet_state: Outlet gas state
            
        Returns:
            Residence time (s)
        """
        # Residence time = V * rho / m_dot
        rho = outlet_state.density
        residence_time = self.volume_m3 * rho / m_dot
        return residence_time
