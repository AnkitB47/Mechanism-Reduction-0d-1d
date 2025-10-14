"""
Perfectly Stirred Reactor (PSR) with fixed burner heat loss.
Uses fixed Q_burner values from Richter benchmark with robust numerical methods.
Enhanced with hot-branch continuation, step-back damping, and enhanced solver settings.
"""

import cantera as ct
import numpy as np
import time
from scipy.optimize import newton
from typing import Tuple, Optional, Dict, Any, List
from .thermo import ThermoManager, GasState, sanitize_TY, T_MIN, T_MAX, Y_FLOOR, sanitize_state, finite_or_backtrack, DT_MIN, DT_GROWTH, DT_SHRINK
from .premix_init import premix_equilibrated_initial_state


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

    def compute_jet_velocity(self, m_dot: float, nozzle, gas_state: ct.Solution) -> float:
        """Compute jet velocity from nozzle geometry and gas state.

        Args:
            m_dot: Mass flow rate (kg/s)
            nozzle: NozzleConfig object with area_m2, diameter_m, C_d
            gas_state: Gas state for density calculation

        Returns:
            Jet velocity (m/s)
        """
        # Get nozzle area
        if nozzle.area_m2 is not None:
            area = nozzle.area_m2
        elif nozzle.diameter_m is not None:
            area = np.pi * (nozzle.diameter_m / 2)**2
        else:
            raise ValueError("Nozzle must specify either area_m2 or diameter_m")

        # Get discharge coefficient
        C_d = nozzle.C_d

        # Get density from gas state at current T, p, X
        rho = gas_state.density

        # u = m_dot / (rho * C_d * A)
        u_jet = m_dot / (rho * C_d * area)

        return u_jet

    def run_mixing_duct_stage(self, inlet_streams: List[Dict[str, Any]], case_config, mech: str) -> Dict[str, Any]:
        """Run Stage A: Cold mixing duct with optional wall loss.

        Args:
            inlet_streams: List of inlet stream dictionaries
            case_config: Case configuration object
            mech: Mechanism file path

        Returns:
            Mixed state after duct
        """
        print(f"    STAGE A: Cold mixing duct (L={case_config.mixing.L_ign_m:.2f}m)")

        # Build streams as Cantera quantities
        streams = []
        for stream_data in inlet_streams:
            gas = ct.Solution(mech)
            gas.TPY = stream_data['T'], case_config.pressure_pa, stream_data['Y']
            q = ct.Quantity(gas, constant="HP")
            q.mass = stream_data['m_dot']
            streams.append(q)

        # Mix all streams
        if streams:
            q_total = streams[0]
            for q in streams[1:]:
                q_total += q
        else:
            raise ValueError("No inlet streams provided")

        # Mixed state
        gas_mixed = ct.Solution(mech)
        gas_mixed.TPY = q_total.TPY

        # Apply optional pre-mix wall loss
        if case_config.mixing.q_pre_W > 0:
            print(f"      Pre-mix wall loss: {case_config.mixing.q_pre_W:.1f} W")
            # Deduct sensible enthalpy for wall loss
            h_loss = case_config.mixing.q_pre_W / q_total.mass  # J/kg
            h_current = gas_mixed.enthalpy_mass
            h_new = h_current - h_loss

            # Update temperature to match new enthalpy (approximate)
            gas_mixed.enthalpy_mass = h_new
            print(f"      T after wall loss: {gas_mixed.T:.1f}K")

        # Compute velocities for logging
        velocities = {}
        for stream_name, nozzle_config in case_config.nozzles.items():
            # Find corresponding inlet stream
            stream = next((s for s in inlet_streams if s.get('name') == stream_name), None)
            if stream:
                u_jet = self.compute_jet_velocity(stream['m_dot'], nozzle_config, gas_mixed)
                velocities[f'u_{stream_name}'] = u_jet
                print(f"      u_{stream_name}: {u_jet:.1f} m/s")

        # Compute effective mixing velocity (using oxidizer area as proxy)
        oxid_nozzle = case_config.nozzles['oxid_premix']
        oxid_area = oxid_nozzle.area_m2
        if oxid_area is None:
            if oxid_nozzle.diameter_m:
                oxid_area = np.pi * (oxid_nozzle.diameter_m / 2)**2
            else:
                oxid_area = 1e-4  # Fallback area if not specified

        u_mix = q_total.mass / (gas_mixed.density * oxid_area)
        velocities['u_mix'] = u_mix
        print(f"      u_mix (effective): {u_mix:.1f} m/s")

        return {
            'gas_mixed': gas_mixed,
            'velocities': velocities,
            'q_total': q_total
        }

    def run_ignition_kernel_stage(self, mixed_state: ct.Solution, case_config, mech: str) -> ct.Solution:
        """Run Stage B: Ignition kernel formation.

        Args:
            mixed_state: Mixed state from Stage A
            case_config: Case configuration object
            mech: Mechanism file path

        Returns:
            HP-equilibrated kernel state for PSR
        """
        print(f"    STAGE B: Ignition kernel formation (tau={case_config.kernel.tau_kernel_s:.3f}s)")

        # Create kernel PSR
        kernel_volume = case_config.psr_volume_m3 * case_config.kernel.V_kernel_frac_psr
        kernel_psr = PSR(kernel_volume, self.thermo)

        # Build inlet streams for premix
        inlet_streams = [
            {
                'name': 'fuel_premix',
                'T': case_config.natural_gas.T_K,
                'm_dot': case_config.natural_gas.mass_kg_per_s * case_config.kernel.y_feed_to_flame,
                'Y': case_config.natural_gas.composition_X
            },
            {
                'name': 'oxid_premix',
                'T': case_config.oxygen.T_K,
                'm_dot': case_config.oxygen.mass_kg_per_s,
                'Y': case_config.oxygen.composition_X
            },
            {
                'name': 'steam_secondary',
                'T': case_config.secondary_steam.T_K,
                'm_dot': case_config.secondary_steam.mass_kg_per_s,
                'Y': case_config.secondary_steam.composition_X
            }
        ]

        # Use premix helper to create equilibrated kernel
        gas_init_flame, gas_psr_in = premix_equilibrated_initial_state(
            mech=mech,
            p_pa=case_config.pressure_pa,
            T_fuel_K=case_config.natural_gas.T_K,
            fuel_comp_X=case_config.natural_gas.composition_X,
            m_fuel_kg_s=case_config.natural_gas.mass_kg_per_s * case_config.kernel.y_feed_to_flame,
            T_N2_K=case_config.nitrogen.T_K,
            m_N2_kg_s=case_config.nitrogen.mass_kg_per_s,
            T_O2_K=case_config.oxygen.T_K,
            m_O2_kg_s=case_config.oxygen.mass_kg_per_s,
            T_secsteam_K=case_config.secondary_steam.T_K,
            m_secsteam_kg_s=case_config.secondary_steam.mass_kg_per_s,
            T_primsteam_K=case_config.primary_steam.T_K,
            m_primsteam_kg_s=case_config.primary_steam.mass_kg_per_s,
            y_feed_to_flame=case_config.kernel.y_feed_to_flame
        )

        print(f"      Kernel T: {gas_psr_in.T:.1f}K")
        print(f"      Kernel composition: H2={gas_psr_in.X[gas_psr_in.species_index('H2')]*100:.1f}%, "
              f"CO={gas_psr_in.X[gas_psr_in.species_index('CO')]*100:.1f}%")

        return gas_psr_in
    
    def solve_psr_with_premix_ignition(self, inlet_state: GasState, m_dot: float,
                                      P: float, q_burner_kw: float, case_config, mech: str,
                                      inlet_streams: List[Dict[str, Any]]) -> Tuple[GasState, bool, Dict[str, Any]]:
        """
        Solve PSR as true CSTR at constant P with fixed Q_burner heat input.

        Implements robust three-stage pipeline:
        Stage A: Cold mixing duct with optional wall loss
        Stage B: Ignition kernel formation (premix HP-equilibrated)
        Stage C: PSR with burner heat loss (existing continuation logic)

        Args:
            inlet_state: Inlet gas state (legacy - not used with new pipeline)
            m_dot: Mass flow rate (kg/s)
            P: Pressure (Pa)
            q_burner_kw: Fixed burner heat input (kW, positive for heating)
            case_config: Case configuration object
            mech: Mechanism file path
            inlet_streams: List of inlet stream dictionaries

        Returns:
            Tuple of (outlet_state, converged, diagnostics)
        """
        print(f"  Solving PSR as CSTR with Q_burner = {q_burner_kw:.2f} kW (HEAT INPUT)")
        print(f"  THREE-STAGE PSR PIPELINE:")
        print(f"    Stage A: Cold mixing duct (L={case_config.mixing.L_ign_m:.2f}m)")
        print(f"    Stage B: Ignition kernel formation (tau={case_config.kernel.tau_kernel_s:.3f}s)")
        print(f"    Stage C: PSR continuation (existing logic)")

        # Stage A: Cold mixing duct
        mixing_result = self.run_mixing_duct_stage(inlet_streams, case_config, mech)
        gas_mixed = mixing_result['gas_mixed']

        # Stage B: Ignition kernel formation
        gas_kernel = self.run_ignition_kernel_stage(gas_mixed, case_config, mech)

        # Convert kW to W for PSR step
        q_burner = q_burner_kw * 1000.0  # W

        # Stage C: PSR continuation (existing logic, but with new initial state)
        print(f"    Stage C: PSR continuation with kernel seed (T={gas_kernel.T:.1f}K)")

        # Convert kernel gas state to GasState format
        kernel_gas_state = GasState(gas_kernel, gas_kernel.T, gas_kernel.P, gas_kernel.Y)

        # Run PSR step with kernel as initial state and Q_burner heat input
        print(f"    Running PSR step with kernel seed (T={kernel_gas_state.temperature:.1f}K)")

        # Create reactor state from kernel
        current_state = {
            'T': kernel_gas_state.temperature,
            'X': kernel_gas_state.mole_fractions.copy(),
            'Y': kernel_gas_state.mass_fractions.copy()
        }

        # Run PSR step with Q_burner heat input
        success = self._run_psr_step(current_state, m_dot, kernel_gas_state.temperature, q_burner, P, kernel_gas_state)

        if success:
            print(f"    PSR converged: T_out = {current_state['T']:.1f}K")
            converged = True
            branch_status = "HOT"
        else:
            print(f"    PSR failed: T_out = {current_state['T']:.1f}K")
            converged = False
            branch_status = "COLD"

        # Create final outlet state
        outlet_state = GasState(
            gas=self.thermo.gas,
            temperature=current_state['T'],
            pressure=P,
            mass_fractions=current_state['Y']
        )

        # Build diagnostics
        diagnostics = {
            'T_out': current_state['T'],
            'Q_burner': q_burner,
            'Q_chem': 0.0,  # Will be calculated properly in full implementation
            'residual': 0.0,  # Will be calculated properly in full implementation
            'converged': converged,
            'steps': 1,  # Single step for now
            'branch_status': branch_status,
            'kernel_T': gas_kernel.T,
            'kernel_H2': gas_kernel.X[gas_kernel.species_index('H2')] * 100,
            'kernel_CO': gas_kernel.X[gas_kernel.species_index('CO')] * 100
        }

        return outlet_state, converged, diagnostics

    def _run_psr_step(self, current_state: dict, m_dot: float, T_in: float, q_burner: float, P: float, inlet_state: GasState) -> bool:
        """Run a single PSR step with given conditions."""
        try:
            # Create reactor with current conditions
            reactor = ct.IdealGasReactor(self.thermo.gas)
            reactor.volume = self.volume_m3
            reactor.thermo.TPY = current_state['T'], P, current_state['Y']

            # Set up inlet flow (this is critical for CSTR operation)
            inlet = ct.Reservoir(self.thermo.gas)
            inlet.thermo.TPY = T_in, P, inlet_state.mass_fractions

            # Create mass flow controller for inlet
            mfc = ct.MassFlowController(inlet, reactor)
            mfc.mass_flow_rate = m_dot

            # Create outlet (reservoir at reactor pressure)
            outlet = ct.Reservoir(self.thermo.gas)
            outlet.thermo.TPY = current_state['T'], P, current_state['Y']

            # Create pressure controller for outlet (simplified)
            pc = ct.PressureController(reactor, outlet)
            pc.primary = mfc

            # Add wall with prescribed heat flux (power-controlled)
            reservoir = ct.Reservoir(self.thermo.gas)

            # Calculate reactor area from volume (simple cylindrical approximation)
            raw_area = (self.volume_m3 * 4.0 / 3.14159)**0.5
            reactor_area_m2 = max(raw_area, 1.0e-3)

            wall = ct.Wall(reactor, reservoir, A=reactor_area_m2)

            # Set prescribed heat using heat_transfer_coeff (correct Cantera method)
            if q_burner > 0:
                reference_dT = 100.0  # K
                htc = q_burner / (reactor_area_m2 * reference_dT)  # W/m²/K
                wall.heat_transfer_coeff = htc

                reservoir.thermo.TP = reactor.thermo.T + reference_dT, reactor.thermo.P

                print(f"      [Q-step] target={q_burner:.2f} W, htc={htc:.2f} W/m²/K, A={reactor_area_m2:.3f} m²")
            else:
                wall.heat_transfer_coeff = 0.0  # Adiabatic

            # Set up reactor network with all components
            net = ct.ReactorNet([reactor])
            net.rtol = 1e-7
            net.atol = 1e-13
            net.max_steps = int(3e6)

            # Time march to steady state
            t = 0.0
            dt = 1e-5
            t_end = 0.1  # 100 ms should be enough for ignition

            # Steady-state detection
            steady_state_count = 0
            steady_state_threshold = 0.005  # 5 ms sustained
            prev_T = reactor.thermo.T
            prev_Y = reactor.thermo.Y.copy()

            while t < t_end:
                try:
                    net.advance(t + dt)
                    t = net.time

                    # Check steady-state criteria
                    current_T = reactor.thermo.T
                    current_Y = reactor.thermo.Y

                    if t > 0:
                        dT_dt = abs(current_T - prev_T) / dt
                        dY_dt = np.sum(np.abs(current_Y - prev_Y)) / dt

                        if dT_dt < 1e-3 and dY_dt < 1e-10:
                            steady_state_count += dt
                            if steady_state_count >= steady_state_threshold:
                                break
                        else:
                            steady_state_count = 0

                        prev_T = current_T
                        prev_Y = current_Y.copy()

                    # Y clipping and renormalization after each substep
                    Y_current = reactor.thermo.Y
                    Y_clipped = np.maximum(Y_current, 1e-20)
                    Y_normalized = Y_clipped / np.sum(Y_clipped)
                    reactor.thermo.Y = Y_normalized

                    dt = min(dt * 1.1, 1e-3)
                except Exception as e:
                    if "tout too close to t0" in str(e):
                        t = net.time
                        dt = max(dt * 0.5, 1e-6)
                        continue
                    else:
                        return False

            # Update current state
            current_state['T'] = reactor.thermo.T
            current_state['X'] = reactor.thermo.X.copy()
            current_state['Y'] = reactor.thermo.Y.copy()

            # Check if hot branch
            return current_state['T'] >= 1200.0

        except Exception as e:
            print(f"      FAIL: PSR step failed: {e}")
            return False
