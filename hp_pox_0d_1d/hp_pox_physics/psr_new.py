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
        
        # Initialize watchdog
        watchdog = Watchdog(wall_time_limit=180.0, max_steps=2000)
        
        # Convert kW to W
        q_burner_w = q_burner_kw * 1000.0
        
        # Calculate residence time
        res_time_s = self.volume_m3 * inlet_state.density / m_dot
        dt_phys = max(0.2 * res_time_s, 0.001)  # Physical time step with minimum
        
        # A) Ignition-friendly seeding
        T_seed = max(1100.0, inlet_state.temperature + 250.0)
        T_seed = min(T_seed, 1400.0)  # Cap at 1400K
        print(f"  Ignition seeding: T_seed = {T_seed:.1f}K")
        
        # Auto-ignition phase with Q=0
        print(f"  Auto-ignition phase: Q=0, up to 10*τ = {10*res_time_s:.3f}s")
        
        max_omega = 0.0
        max_qdot_chem = 0.0
        time_to_reactive = 0.0
        attempts = 0
        reactive_seed = False
        
        while attempts < 6 and not reactive_seed:  # Max 6 attempts
            attempts += 1
            print(f"    Attempt {attempts}: T_seed = {T_seed:.1f}K")
            
            # Create fresh reactor for each attempt
            reactor = ct.IdealGasConstPressureReactor(self.thermo.gas)
            reactor.volume = self.volume_m3
            reactor.thermo.TPY = T_seed, P, inlet_state.mass_fractions
            
            # Set up CVODE/BDF solver
            sim = ct.ReactorNet([reactor])
            sim.rtol = 1e-8
            sim.atol = 1e-14
            sim.max_steps = int(2e6)
            
            # Set up vector atols
            n_species = len(reactor.thermo.species_names)
            atol_vector = np.full(n_species + 1, 1e-14)  # +1 for temperature
            atol_vector[0] = 1e-9  # Temperature atol
            
            # Set stricter atols for radicals
            radical_names = ['H', 'O', 'OH', 'HO2', 'H2O2', 'CH', 'CH2', 'CH3']
            for i, name in enumerate(reactor.thermo.species_names):
                if name in radical_names:
                    atol_vector[i + 1] = 1e-18  # +1 because T is first
            
            try:
                sim.atol = atol_vector
            except:
                pass  # Use scalar atol if vector not supported
            
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
                    
                    current_max_omega = float(np.max(np.abs(omega)))
                    current_max_qdot = abs(qdot_chem_vol)
                    
                    max_omega = max(max_omega, current_max_omega)
                    max_qdot_chem = max(max_qdot_chem, current_max_qdot)
                    
                    if current_max_omega > 1e-3 or current_max_qdot > 1e6:
                        time_to_reactive = t
                        print(f"    Reactive at t={t:.3f}s: max|ω|={current_max_omega:.2e}, |qdot|={current_max_qdot:.2e}")
                        reactive_seed = True
                        break
                        
                except Exception as e:
                    print(f"    Auto-ignition failed at t={t:.3f}s: {e}")
                    break
            
            # Check if reactive enough
            if max_omega > 1e-4 and max_qdot_chem > 1e5:
                print(f"    Ignition successful: max|ω|={max_omega:.2e}, |qdot|={max_qdot_chem:.2e}")
                reactive_seed = True
                break
            else:
                T_seed += 50.0  # Bump temperature
                if T_seed > 1400.0:
                    print(f"    Warning: Max T_seed reached, proceeding with current state")
                    break
        
        if not reactive_seed:
            raise RuntimeError("PSR: no reactive seed - auto-ignition phase failed to achieve reactivity")
        
        # B) Q_burner ramping with geometric increments (15-20 increments)
        print(f"  Q_burner ramping: 5% → 100% of target")
        q_increments = [0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.82, 0.88, 0.93, 0.97, 1.00]
        q_burner_increments = [q * q_burner_w for q in q_increments]
        
        # Create wall for heat loss
        env = ct.Reservoir(self.thermo.gas)
        wall = ct.Wall(reactor, env)
        
        print(f"  Residence time: {res_time_s:.3f} s, dt_phys: {dt_phys:.3f} s")
        print(f"  Q_burner increments: {[f'{q/1000:.1f}' for q in q_burner_increments]} kW")
        
        converged = True
        final_residual = float('inf')
        min_residual = float('inf')
        
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
            'max_qdot_chem_seeding': max_qdot_chem,
            'min_residual': float('inf'),
            'clips_applied': 0
        }
        
        # Ramp Q_burner to target value
        for i, q_current in enumerate(q_burner_increments):
            print(f"  Increment {i+1}/{len(q_burner_increments)}: Q_burner = {q_current/1000:.2f} kW")
            
            # Set heat transfer coefficient (negative for cooling)
            T_diff = max(reactor.thermo.T - 300.0, 100.0)  # Minimum temperature difference
            wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
            
            # Reinitialize solver for each increment
            sim.reinitialize()
            
            # Bounded inner loop: up to 50 chunks of dt_phys
            residual_history = []
            best_residual = float('inf')
            worsening_count = 0
            stagnation_count = 0
            
            for k in range(50):  # Max 50 steps per increment
                watchdog.should_abort(f"PSR increment {i+1}, step {k+1}")
                
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
                    current_max_omega = float(np.max(np.abs(omega)))
                except:
                    q_chem_w = 0.0
                    current_max_omega = 0.0
                
                residual = abs(q_chem_w + q_current + h_dot_out - h_dot_in) / max(abs(h_dot_in), 1.0)
                residual_history.append(residual)
                min_residual = min(min_residual, residual)
                
                # Log progress
                watchdog.log_progress(reactor.thermo.T, residual, current_max_omega)
                
                # Check for convergence
                if residual < 1e-3 and current_max_omega > 1e-4:
                    print(f"    Converged at step {k+1}: residual = {residual:.2e}, max|ω| = {current_max_omega:.2e}")
                    break
                
                # Check for worsening residual
                if residual > best_residual * 1.1:  # 10% increase
                    worsening_count += 1
                    if worsening_count >= 3:
                        print(f"    Residual worsening at step {k+1}, halving increment")
                        q_current *= 0.5
                        T_diff = max(reactor.thermo.T - 300.0, 100.0)
                        wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
                        sim.reinitialize()
                        worsening_count = 0
                        continue
                else:
                    worsening_count = 0
                
                # Check for stagnation
                if len(residual_history) >= 5:
                    recent_improvement = (residual_history[-5] - residual) / residual_history[-5]
                    if recent_improvement < 0.01:  # <1% improvement over 5 steps
                        stagnation_count += 1
                        if stagnation_count >= 1:  # Allow only 1 stagnation before halving
                            print(f"    Stagnation detected at step {k+1}, halving increment")
                            q_current *= 0.5
                            T_diff = max(reactor.thermo.T - 300.0, 100.0)
                            wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
                            sim.reinitialize()
                            stagnation_count = 0
                            continue
                    else:
                        stagnation_count = 0
                
                best_residual = min(best_residual, residual)
                
                # Advance in time with max horizon guard
                t_target = min(sim.time + dt_phys, sim.time + 2 * res_time_s)
                
                try:
                    sim.advance(t_target)
                    
                    # Enforce bounds
                    T = reactor.thermo.T
                    if T < 300.0 or T > 3500.0:
                        print(f"    Temperature out of bounds: T={T:.1f}K, clipping")
                        T = max(300.0, min(3500.0, T))
                        reactor.thermo.T = T
                        diagnostics['clips_applied'] += 1
                    
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
                        T_diff = max(reactor.thermo.T - 300.0, 100.0)
                        wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
                        sim.reinitialize()
                        continue
                    
                except Exception as e:
                    print(f"    Integration failed at step {k+1}: {e}")
                    q_current *= 0.5
                    T_diff = max(reactor.thermo.T - 300.0, 100.0)
                    wall.heat_transfer_coeff = -q_current / (reactor.volume * T_diff)
                    sim.reinitialize()
                    continue
            
            # Check if this increment converged
            if residual >= 1e-3 or current_max_omega <= 1e-4:
                print(f"    Increment {i+1} did not converge: residual = {residual:.2e}, max|ω| = {current_max_omega:.2e}")
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
            'final_heat_transfer_coeff': wall.heat_transfer_coeff,
            'min_residual': min_residual
        })
        
        # Print single line summary
        status = "CONVERGED" if converged else "FAILED"
        print(f"  PSR: Q_burner={q_burner_kw:.2f} kW | Q_chem={q_chem_w_out/1000:.2f} kW | Residual={final_residual:.2e} | T_out={T_out:.1f}K | {status}")
        
        return outlet_state, diagnostics, converged
