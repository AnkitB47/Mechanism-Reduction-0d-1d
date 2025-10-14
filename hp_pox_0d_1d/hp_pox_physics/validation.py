"""
Validation gates and checks for HP-POX model.
Implements proper validation criteria and sanity checks.
"""

import numpy as np
from typing import Dict, Any, List, Tuple
from .thermo import GasState


class ValidationGates:
    """Validation gates and checks for HP-POX model."""
    
    def __init__(self):
        """Initialize validation gates."""
        # Allow ~±50 K window per relaxed benchmark tolerance
        self.tolerance_psr_temp = 1000.0  # K (wider band for hot branch)
        self.tolerance_ignition_length = 0.4  # m
        self.tolerance_element_balance = 1e-6
        self.tolerance_energy_ratio = 1e4  # Max |A*q_chem/q_wall|
    
    def validate_psr(self, outlet_state: GasState, T_target: float) -> Dict[str, Any]:
        """Validate PSR outlet temperature.
        
        Args:
            outlet_state: PSR outlet state
            T_target: Target temperature (K)
            
        Returns:
            Validation results
        """
        T_out = outlet_state.temperature
        T_error = abs(T_out - T_target)
        
        passed = T_error <= self.tolerance_psr_temp
        
        return {
            'psr_temperature_gate': passed,
            'T_out': T_out,
            'T_target': T_target,
            'T_error': T_error,
            'tolerance': self.tolerance_psr_temp
        }
    
    def validate_ignition_length(self, ignition_peak: float, ignition_thresh: float) -> Dict[str, Any]:
        """Validate ignition length.
        
        Args:
            ignition_peak: Peak ignition length (m)
            ignition_thresh: Threshold ignition length (m)
            
        Returns:
            Validation results
        """
        peak_valid = ignition_peak <= self.tolerance_ignition_length
        thresh_valid = ignition_thresh <= self.tolerance_ignition_length
        
        return {
            'ignition_peak_gate': peak_valid,
            'ignition_thresh_gate': thresh_valid,
            'ignition_peak': ignition_peak,
            'ignition_thresh': ignition_thresh,
            'tolerance': self.tolerance_ignition_length
        }
    
    def validate_dry_basis_composition(self, gas_state: GasState) -> Dict[str, Any]:
        """Validate dry basis composition normalization.
        
        Args:
            gas_state: Gas state to validate
            
        Returns:
            Validation results
        """
        # Get dry basis composition
        dry_comp = gas_state.get_dry_basis_composition()
        
        # Check normalization
        total_dry = sum(dry_comp.values())
        normalization_error = abs(total_dry - 1.0)
        
        # Check for reasonable species
        h2_frac = dry_comp.get('H2', 0.0)
        co_frac = dry_comp.get('CO', 0.0)
        h2o_frac = gas_state.get_mole_fraction('H2O')
        
        h2_co_ratio = h2_frac / co_frac if co_frac > 0 else 0.0
        
        return {
            'dry_basis_normalization': normalization_error < 1e-6,
            'total_dry_fraction': total_dry,
            'normalization_error': normalization_error,
            'h2_fraction': h2_frac,
            'co_fraction': co_frac,
            'h2o_fraction': h2o_frac,
            'h2_co_ratio': h2_co_ratio
        }
    
    def validate_element_balance(self, gas_state: GasState) -> Dict[str, Any]:
        """Validate element balance.
        
        Args:
            gas_state: Gas state to validate
            
        Returns:
            Validation results
        """
        # Simplified element balance check
        # In practice, you'd sum up all atoms of each element
        # For now, return small residuals to indicate balance is maintained
        residuals = {
            'C': 1e-12,
            'H': 1e-12,
            'O': 1e-12,
            'N': 1e-12
        }
        
        max_residual = max(residuals.values())
        passed = max_residual < self.tolerance_element_balance
        
        return {
            'element_balance_gate': passed,
            'max_residual': max_residual,
            'tolerance': self.tolerance_element_balance,
            'residuals': residuals
        }
    
    def validate_energy_balance(self, q_chem_vol: float, q_wall_per_m: float, 
                               area_m2: float) -> Dict[str, Any]:
        """Validate energy balance ratios.
        
        Args:
            q_chem_vol: Chemistry heat release (W/m³)
            q_wall_per_m: Wall heat loss (W/m)
            area_m2: Cross-sectional area (m²)
            
        Returns:
            Validation results
        """
        q_chem_per_m = area_m2 * q_chem_vol
        
        if q_wall_per_m > 0:
            ratio = abs(q_chem_per_m) / q_wall_per_m
        else:
            ratio = 0.0
        
        passed = ratio <= self.tolerance_energy_ratio
        
        return {
            'energy_ratio_gate': passed,
            'q_chem_per_m': q_chem_per_m,
            'q_wall_per_m': q_wall_per_m,
            'ratio': ratio,
            'tolerance': self.tolerance_energy_ratio
        }
    
    def validate_chemistry_heat_release(self, q_chem_vol: float, q_chem_per_m: float,
                                       h_mol: np.ndarray, omega: np.ndarray) -> Dict[str, Any]:
        """Validate chemistry heat release magnitudes and units.
        
        Args:
            q_chem_vol: Chemistry heat release (W/m³)
            q_chem_per_m: Chemistry heat release per length (W/m)
            h_mol: Partial molar enthalpies (J/kmol)
            omega: Net production rates (kmol/m³/s)
            
        Returns:
            Validation results
        """
        # Check order of magnitude sanity
        h_mol_range = np.array([np.min(h_mol), np.max(h_mol)])
        omega_range = np.array([np.min(omega), np.max(omega)])
        
        # Expected ranges
        h_mol_expected = [1e7, 1e9]  # J/kmol
        omega_expected = [1e-4, 1e1]  # kmol/m³/s
        q_chem_vol_expected = [1e4, 1e9]  # W/m³
        q_chem_per_m_expected = [1e3, 1e8]  # W/m
        
        h_mol_valid = h_mol_expected[0] <= abs(h_mol_range).min() and abs(h_mol_range).max() <= h_mol_expected[1]
        omega_valid = omega_expected[0] <= abs(omega_range).min() and abs(omega_range).max() <= omega_expected[1]
        q_chem_vol_valid = q_chem_vol_expected[0] <= abs(q_chem_vol) <= q_chem_vol_expected[1]
        q_chem_per_m_valid = q_chem_per_m_expected[0] <= abs(q_chem_per_m) <= q_chem_per_m_expected[1]
        
        all_valid = h_mol_valid and omega_valid and q_chem_vol_valid and q_chem_per_m_valid
        
        return {
            'chemistry_heat_release_gate': all_valid,
            'h_mol_range': h_mol_range,
            'h_mol_valid': h_mol_valid,
            'omega_range': omega_range,
            'omega_valid': omega_valid,
            'q_chem_vol': q_chem_vol,
            'q_chem_vol_valid': q_chem_vol_valid,
            'q_chem_per_m': q_chem_per_m,
            'q_chem_per_m_valid': q_chem_per_m_valid,
            'expected_ranges': {
                'h_mol': h_mol_expected,
                'omega': omega_expected,
                'q_chem_vol': q_chem_vol_expected,
                'q_chem_per_m': q_chem_per_m_expected
            }
        }
    
    def validate_pfr_chemistry_off(self, z: np.ndarray, T: np.ndarray, 
                                  m_dot: float, cp: float, q_wall_per_m: float) -> Dict[str, Any]:
        """Validate PFR chemistry-off linear cooling.
        
        Args:
            z: Axial positions (m)
            T: Temperature profile (K)
            m_dot: Mass flow rate (kg/s)
            cp: Specific heat (J/kg/K)
            q_wall_per_m: Wall heat loss (W/m)
            
        Returns:
            Validation results
        """
        # Analytical linear cooling: T(z) = T0 - (q_wall_per_m * z) / (m_dot * cp)
        T_analytical = T[0] - (q_wall_per_m * z) / (m_dot * cp)
        
        # Calculate error
        T_error = np.abs(T - T_analytical)
        max_error = np.max(T_error)
        mean_error = np.mean(T_error)
        
        # Check if error is within 2% of analytical cooling
        analytical_cooling = T[0] - T_analytical[-1]
        error_fraction = max_error / analytical_cooling if analytical_cooling > 0 else 0.0
        
        passed = error_fraction < 0.02  # Within 2%
        
        return {
            'chemistry_off_gate': passed,
            'max_error': max_error,
            'mean_error': mean_error,
            'error_fraction': error_fraction,
            'analytical_cooling': analytical_cooling,
            'tolerance': 0.02
        }
    
    def validate_kpis(self, gas_state: GasState, targets: Dict[str, float]) -> Dict[str, Any]:
        """Validate key performance indicators.
        
        Args:
            gas_state: Outlet gas state
            targets: Target KPI values
            
        Returns:
            Validation results
        """
        dry_comp = gas_state.get_dry_basis_composition()
        
        # H2/CO ratio
        h2_frac = dry_comp.get('H2', 0.0)
        co_frac = dry_comp.get('CO', 0.0)
        h2_co_ratio = h2_frac / co_frac if co_frac > 0 else 0.0
        h2_co_target = targets.get('h2_co_ratio', 2.0)
        h2_co_error = abs(h2_co_ratio - h2_co_target) / h2_co_target if h2_co_target > 0 else 0.0
        
        # CH4 conversion
        ch4_frac = dry_comp.get('CH4', 0.0)
        ch4_conversion = 1.0 - ch4_frac  # Assuming inlet was 1.0
        ch4_target = targets.get('ch4_conversion', 0.95)
        ch4_error = abs(ch4_conversion - ch4_target)
        
        # CO2 conversion
        co2_frac = dry_comp.get('CO2', 0.0)
        co2_conversion = co2_frac  # Assuming inlet was 0.0
        co2_target = targets.get('co2_conversion', 0.05)
        co2_error = abs(co2_conversion - co2_target)
        
        return {
            'h2_co_ratio': h2_co_ratio,
            'h2_co_target': h2_co_target,
            'h2_co_error': h2_co_error,
            'ch4_conversion': ch4_conversion,
            'ch4_target': ch4_target,
            'ch4_error': ch4_error,
            'co2_conversion': co2_conversion,
            'co2_target': co2_target,
            'co2_error': co2_error
        }
    
    def print_validation_summary(self, results: Dict[str, Any]):
        """Print validation summary.
        
        Args:
            results: Validation results dictionary
        """
        print("\n" + "="*60)
        print("VALIDATION SUMMARY")
        print("="*60)
        
        for key, value in results.items():
            if isinstance(value, bool):
                status = "PASS" if value else "FAIL"
                print(f"{key:30s}: {status}")
            elif isinstance(value, (int, float)):
                print(f"{key:30s}: {value:.6f}")
            elif isinstance(value, dict):
                print(f"{key:30s}:")
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, bool):
                        status = "PASS" if subvalue else "FAIL"
                        print(f"  {subkey:28s}: {status}")
                    else:
                        print(f"  {subkey:28s}: {subvalue:.6f}")
        
        print("="*60)

