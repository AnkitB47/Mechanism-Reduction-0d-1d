"""
Key Performance Indicators (KPIs) for HP-POX model.

Calculates dry-basis compositions, conversions, and errors vs targets.
"""

import numpy as np
from typing import Dict, Any, List


class KPI:
    """KPI calculator for HP-POX model."""
    
    def __init__(self, case_config: Dict[str, Any]):
        """Initialize KPI calculator.
        
        Args:
            case_config: Case configuration
        """
        self.case_config = case_config
        self.targets = case_config['targets']
    
    def calculate(self, psr_result: Dict[str, Any], pfr_result: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate all KPIs.
        
        Args:
            psr_result: PSR results
            pfr_result: PFR results
            
        Returns:
            Dictionary with all KPIs
        """
        # Outlet conditions
        outlet_idx = -1
        T_out = pfr_result['temperature_K'][outlet_idx]
        P_out = pfr_result['pressure_Pa'][outlet_idx]
        Y_out = pfr_result['mass_fractions'][:, outlet_idx]
        X_out = pfr_result['mole_fractions'][:, outlet_idx]
        
        # Calculate KPIs
        h2_co_ratio = self._calculate_h2_co_ratio(X_out, pfr_result['species_names'])
        ch4_conversion = self._calculate_ch4_conversion(X_out, pfr_result['species_names'])
        co2_yield = self._calculate_co2_yield(X_out, pfr_result['species_names'])
        
        # Dry composition
        dry_composition = self._calculate_dry_composition(X_out, pfr_result['species_names'])
        wet_composition = self._calculate_wet_composition(X_out, pfr_result['species_names'])
        
        # Errors vs targets
        errors_pp = self._calculate_composition_errors(dry_composition)
        
        # Other KPIs
        ignition_peak = pfr_result['ignition_length_peak_m']
        ignition_thresh = pfr_result['ignition_length_threshold_m']
        pressure_drop = pfr_result['pressure_drop_Pa']
        residence_time_total = pfr_result['residence_time_s'][-1]
        residence_time_psr = psr_result['residence_time_s']
        
        # Sanity gates
        sanity_gates = self._check_sanity_gates(T_out - 273.15, h2_co_ratio, ignition_peak, 
                                              residence_time_total, pressure_drop, X_out, 
                                              pfr_result['species_names'])
        
        return {
            'h2_co_ratio': h2_co_ratio,
            'ch4_conversion_pct': ch4_conversion,
            'co2_yield_pct': co2_yield,
            'ignition_length_peak_m': ignition_peak,
            'ignition_length_threshold_m': ignition_thresh,
            'pressure_drop_Pa': pressure_drop,
            'residence_time_total_s': residence_time_total,
            'residence_time_psr_s': residence_time_psr,
            'outlet_temperature_C': T_out - 273.15,
            'outlet_pressure_bar': P_out / 1e5,
            'dry_composition_pct': dry_composition,
            'wet_composition_pct': wet_composition,
            'errors_pp': errors_pp,
            'targets': self.targets,
            'sanity_gates': sanity_gates
        }
    
    def _calculate_h2_co_ratio(self, X_out: np.ndarray, species_names: List[str]) -> float:
        """Calculate H2/CO ratio."""
        h2_idx = species_names.index('H2')
        co_idx = species_names.index('CO')
        
        h2_frac = X_out[h2_idx]
        co_frac = X_out[co_idx]
        
        return h2_frac / co_frac if co_frac > 0 else 0.0
    
    def _calculate_ch4_conversion(self, X_out: np.ndarray, species_names: List[str]) -> float:
        """Calculate CH4 conversion percentage."""
        # Get inlet CH4 fraction from natural gas composition
        ng_comp = self.case_config['feeds']['natural_gas']['composition']
        ch4_inlet_vol_pct = ng_comp['CH4']
        ch4_inlet_mol_frac = ch4_inlet_vol_pct / 100.0
        
        # Get outlet CH4 fraction
        ch4_idx = species_names.index('CH4')
        ch4_outlet_mol_frac = X_out[ch4_idx]
        
        # Calculate conversion
        conversion = (ch4_inlet_mol_frac - ch4_outlet_mol_frac) / ch4_inlet_mol_frac * 100.0
        
        return max(0.0, conversion)  # Don't allow negative conversion
    
    def _calculate_co2_yield(self, X_out: np.ndarray, species_names: List[str]) -> float:
        """Calculate CO2 yield percentage."""
        # Get inlet CO2 fraction from natural gas composition
        ng_comp = self.case_config['feeds']['natural_gas']['composition']
        co2_inlet_vol_pct = ng_comp['CO2']
        co2_inlet_mol_frac = co2_inlet_vol_pct / 100.0
        
        # Get outlet CO2 fraction
        co2_idx = species_names.index('CO2')
        co2_outlet_mol_frac = X_out[co2_idx]
        
        # Calculate yield (CO2 produced / total carbon)
        yield_pct = co2_outlet_mol_frac * 100.0
        
        return yield_pct
    
    def _calculate_dry_composition(self, X_out: np.ndarray, species_names: List[str]) -> Dict[str, float]:
        """Calculate dry composition (exclude H2O only, keep N2)."""
        # Find H2O index
        h2o_idx = species_names.index('H2O')
        
        # Calculate dry total (exclude H2O only)
        dry_total = np.sum(X_out) - X_out[h2o_idx]
        
        # Calculate dry composition
        dry_composition = {}
        target_species = ['H2', 'CO', 'CO2', 'CH4', 'N2']
        
        for species in target_species:
            if species in species_names:
                idx = species_names.index(species)
                dry_composition[species] = X_out[idx] / dry_total * 100.0
            else:
                dry_composition[species] = 0.0
        
        return dry_composition
    
    def _calculate_wet_composition(self, X_out: np.ndarray, species_names: List[str]) -> Dict[str, float]:
        """Calculate wet composition (all species)."""
        wet_composition = {}
        
        for i, species in enumerate(species_names):
            wet_composition[species] = X_out[i] * 100.0
        
        return wet_composition
    
    def _calculate_composition_errors(self, dry_composition: Dict[str, float]) -> Dict[str, float]:
        """Calculate errors vs targets in percentage points."""
        errors = {}
        targets = self.targets['dry_composition_pct']
        
        for species, target in targets.items():
            if species in dry_composition:
                actual = dry_composition[species]
                error = abs(actual - target)
                errors[species] = error
            else:
                errors[species] = 0.0
        
        return errors
    
    def _check_sanity_gates(self, T_out_C: float, h2_co_ratio: float, ignition_peak_m: float,
                           residence_time_s: float, pressure_drop_Pa: float, X_out: np.ndarray,
                           species_names: List[str]) -> Dict[str, Any]:
        """Check sanity gates and return results."""
        gates = {}
        
        # Temperature gate: |T_out - T_target| ≤ 20 K
        T_target = self.targets['outlet_temperature_C']
        T_error = abs(T_out_C - T_target)
        gates['temperature_gate'] = {
            'passed': bool(T_error <= 20.0),
            'value': float(T_out_C),
            'target': float(T_target),
            'error': float(T_error),
            'tolerance': 20.0
        }
        
        # Ignition length gate: 0.2 m ≤ L_ignition_peak ≤ 0.4 m
        gates['ignition_length_gate'] = {
            'passed': bool(0.2 <= ignition_peak_m <= 0.4),
            'value': float(ignition_peak_m),
            'target_range': [0.2, 0.4]
        }
        
        # Residence time gate: 14 s ≤ τ_total ≤ 17 s
        gates['residence_time_gate'] = {
            'passed': bool(14.0 <= residence_time_s <= 17.0),
            'value': float(residence_time_s),
            'target_range': [14.0, 17.0]
        }
        
        # H2/CO ratio gate: 1.8 ≤ (H2/CO)_dry ≤ 2.2
        gates['h2_co_ratio_gate'] = {
            'passed': bool(1.8 <= h2_co_ratio <= 2.2),
            'value': float(h2_co_ratio),
            'target_range': [1.8, 2.2]
        }
        
        # Pressure drop gate: Δp ∈ [1, 20] kPa
        pressure_drop_kPa = pressure_drop_Pa / 1000.0
        gates['pressure_drop_gate'] = {
            'passed': bool(1.0 <= pressure_drop_kPa <= 20.0),
            'value': float(pressure_drop_kPa),
            'target_range': [1.0, 20.0]
        }
        
        # O2 slip gate: X_O2,out < 10^-5
        o2_idx = species_names.index('O2') if 'O2' in species_names else -1
        o2_out = float(X_out[o2_idx]) if o2_idx >= 0 else 0.0
        gates['o2_slip_gate'] = {
            'passed': bool(o2_out < 1e-5),
            'value': o2_out,
            'target': 1e-5
        }
        
        # Overall gate status
        all_passed = all(gate['passed'] for gate in gates.values())
        gates['overall_status'] = 'PASSED' if all_passed else 'FAILED'
        
        return gates
