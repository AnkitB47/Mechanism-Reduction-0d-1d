"""Plasma-Assisted Surrogate Models.

Implements thermal and radical plasma surrogate models for enhanced combustion.
Supports both enthalpy rise and radical seeding approaches.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
import yaml

from reactor.plug_flow import PlugFlowReactor, ReactorGeometry, HeatTransferConfig, InletConditions


@dataclass
class PlasmaConfig:
    """Plasma configuration parameters."""
    
    model_type: str  # "thermal" or "radical"
    power_kW: float
    pressure_bar_g: float
    temperature_C: float
    equivalence_ratio: float
    injection_location: float  # m from inlet
    
    # Thermal model parameters
    enthalpy_rise_kJ_kg: Optional[float] = None
    temperature_rise_C: Optional[float] = None
    
    # Radical model parameters
    radical_fractions: Optional[Dict[str, float]] = None  # H, O, OH, etc.
    radical_temperature_C: Optional[float] = None


@dataclass
class PlasmaResult:
    """Result of plasma-assisted simulation."""
    
    # Base reactor result
    base_result: Dict
    
    # Plasma effects
    plasma_power_kW: float
    enthalpy_rise_kJ_kg: float
    temperature_rise_C: float
    radical_injection: Dict[str, float]
    
    # Enhanced KPIs
    enhanced_ch4_conversion: float
    enhanced_h2_co_ratio: float
    ignition_reduction_pct: float
    efficiency_improvement_pct: float


class ThermalPlasmaSurrogate:
    """Thermal plasma surrogate model.
    
    Adds enthalpy rise to the gas stream at the plasma injection point.
    """
    
    def __init__(self, config: PlasmaConfig):
        """Initialize thermal plasma surrogate."""
        self.config = config
        
    def apply_plasma_effect(
        self,
        reactor: PlugFlowReactor,
        result: Dict
    ) -> PlasmaResult:
        """Apply thermal plasma effect to reactor result.
        
        Args:
            reactor: PlugFlowReactor instance
            result: Base reactor result
            
        Returns:
            Enhanced result with plasma effects
        """
        
        # Find injection point
        injection_idx = self._find_injection_point(reactor)
        
        # Calculate enthalpy rise
        enthalpy_rise = self._calculate_enthalpy_rise(reactor, injection_idx)
        
        # Apply thermal effect
        enhanced_result = self._apply_thermal_effect(reactor, result, injection_idx, enthalpy_rise)
        
        # Calculate enhanced KPIs
        kpis = self._calculate_enhanced_kpis(reactor, enhanced_result, result)
        
        return PlasmaResult(
            base_result=result,
            plasma_power_kW=self.config.power_kW,
            enthalpy_rise_kJ_kg=enthalpy_rise,
            temperature_rise_C=enhanced_result["temperature_rise"],
            radical_injection={},
            enhanced_ch4_conversion=kpis["ch4_conversion"],
            enhanced_h2_co_ratio=kpis["h2_co_ratio"],
            ignition_reduction_pct=kpis["ignition_reduction"],
            efficiency_improvement_pct=kpis["efficiency_improvement"]
        )
    
    def _find_injection_point(self, reactor: PlugFlowReactor) -> int:
        """Find injection point index."""
        injection_height = self.config.injection_location
        return int(injection_height / reactor.dx)
    
    def _calculate_enthalpy_rise(self, reactor: PlugFlowReactor, injection_idx: int) -> float:
        """Calculate enthalpy rise from plasma power."""
        
        # Get gas properties at injection point
        reactor.gas.TPY = (
            reactor.temperature[injection_idx],
            reactor.pressure[injection_idx],
            reactor.mass_fractions[injection_idx, :]
        )
        
        # Calculate mass flow rate
        area = reactor.geometry.area_at_height(reactor.height[injection_idx])
        mass_flow_rate = reactor.density[injection_idx] * reactor.velocity[injection_idx] * area
        
        # Calculate enthalpy rise
        power_W = self.config.power_kW * 1000.0
        enthalpy_rise = power_W / mass_flow_rate  # J/kg
        
        return enthalpy_rise / 1000.0  # Convert to kJ/kg
    
    def _apply_thermal_effect(
        self,
        reactor: PlugFlowReactor,
        result: Dict,
        injection_idx: int,
        enthalpy_rise: float
    ) -> Dict:
        """Apply thermal effect to reactor result."""
        
        # Calculate temperature rise
        reactor.gas.TPY = (
            reactor.temperature[injection_idx],
            reactor.pressure[injection_idx],
            reactor.mass_fractions[injection_idx, :]
        )
        
        cp = reactor.gas.cp_mass  # J/kg/K
        temperature_rise = enthalpy_rise * 1000.0 / cp  # K
        
        # Create enhanced result
        enhanced_result = result.copy()
        enhanced_result["temperature_rise"] = temperature_rise
        
        # Apply temperature rise downstream
        for i in range(injection_idx, len(reactor.temperature)):
            enhanced_result["temperature"][i] += temperature_rise
        
        return enhanced_result
    
    def _calculate_enhanced_kpis(
        self,
        reactor: PlugFlowReactor,
        enhanced_result: Dict,
        base_result: Dict
    ) -> Dict:
        """Calculate enhanced KPIs."""
        
        # Recalculate composition with enhanced temperature
        # This is a simplified approach - in practice, would re-solve chemistry
        
        # Estimate CH4 conversion improvement
        base_conversion = base_result.get("ch4_conversion", 0.0)
        temp_rise = enhanced_result["temperature_rise"]
        
        # Simple Arrhenius-based improvement
        activation_energy = 150000.0  # J/mol
        gas_constant = 8.314  # J/mol/K
        base_temp = np.mean(reactor.temperature)
        
        rate_factor = np.exp(activation_energy * temp_rise / (gas_constant * base_temp * (base_temp + temp_rise)))
        enhanced_conversion = min(100.0, base_conversion * rate_factor)
        
        # Estimate H2/CO ratio improvement
        base_h2_co = base_result.get("h2_co_ratio", 1.0)
        enhanced_h2_co = base_h2_co * 1.1  # 10% improvement
        
        # Calculate ignition reduction
        base_ignition = base_result.get("ignition_length", 1.0)
        enhanced_ignition = base_ignition * 0.8  # 20% reduction
        ignition_reduction = (base_ignition - enhanced_ignition) / base_ignition * 100.0
        
        # Calculate efficiency improvement
        efficiency_improvement = 5.0  # 5% improvement (estimated)
        
        return {
            "ch4_conversion": enhanced_conversion,
            "h2_co_ratio": enhanced_h2_co,
            "ignition_reduction": ignition_reduction,
            "efficiency_improvement": efficiency_improvement
        }


class RadicalPlasmaSurrogate:
    """Radical plasma surrogate model.
    
    Injects radicals (H, O, OH) at the plasma injection point.
    """
    
    def __init__(self, config: PlasmaConfig):
        """Initialize radical plasma surrogate."""
        self.config = config
        
    def apply_plasma_effect(
        self,
        reactor: PlugFlowReactor,
        result: Dict
    ) -> PlasmaResult:
        """Apply radical plasma effect to reactor result.
        
        Args:
            reactor: PlugFlowReactor instance
            result: Base reactor result
            
        Returns:
            Enhanced result with plasma effects
        """
        
        # Find injection point
        injection_idx = self._find_injection_point(reactor)
        
        # Calculate radical injection rates
        radical_rates = self._calculate_radical_rates(reactor, injection_idx)
        
        # Apply radical effect
        enhanced_result = self._apply_radical_effect(reactor, result, injection_idx, radical_rates)
        
        # Calculate enhanced KPIs
        kpis = self._calculate_enhanced_kpis(reactor, enhanced_result, result)
        
        return PlasmaResult(
            base_result=result,
            plasma_power_kW=self.config.power_kW,
            enthalpy_rise_kJ_kg=0.0,
            temperature_rise_C=0.0,
            radical_injection=radical_rates,
            enhanced_ch4_conversion=kpis["ch4_conversion"],
            enhanced_h2_co_ratio=kpis["h2_co_ratio"],
            ignition_reduction_pct=kpis["ignition_reduction"],
            efficiency_improvement_pct=kpis["efficiency_improvement"]
        )
    
    def _find_injection_point(self, reactor: PlugFlowReactor) -> int:
        """Find injection point index."""
        injection_height = self.config.injection_location
        return int(injection_height / reactor.dx)
    
    def _calculate_radical_rates(
        self,
        reactor: PlugFlowReactor,
        injection_idx: int
    ) -> Dict[str, float]:
        """Calculate radical injection rates."""
        
        # Get gas properties at injection point
        reactor.gas.TPY = (
            reactor.temperature[injection_idx],
            reactor.pressure[injection_idx],
            reactor.mass_fractions[injection_idx, :]
        )
        
        # Calculate mass flow rate
        area = reactor.geometry.area_at_height(reactor.height[injection_idx])
        mass_flow_rate = reactor.density[injection_idx] * reactor.velocity[injection_idx] * area
        
        # Calculate radical injection rates based on power
        power_W = self.config.power_kW * 1000.0
        
        # Estimate radical production efficiency
        efficiency = 0.1  # 10% efficiency (estimated)
        radical_power = power_W * efficiency
        
        # Calculate radical mole fractions
        radical_fractions = self.config.radical_fractions or {
            "H": 0.4,
            "O": 0.3,
            "OH": 0.2,
            "O2": 0.1
        }
        
        # Convert to mass fractions
        radical_rates = {}
        for species, mole_frac in radical_fractions.items():
            if species in reactor.gas.species_names:
                species_idx = reactor.gas.species_index(species)
                mw = reactor.gas.molecular_weights[species_idx]
                mass_frac = mole_frac * mw / np.sum([radical_fractions[s] * reactor.gas.molecular_weights[reactor.gas.species_index(s)] 
                                                   for s in radical_fractions.keys() if s in reactor.gas.species_names])
                radical_rates[species] = mass_frac
        
        return radical_rates
    
    def _apply_radical_effect(
        self,
        reactor: PlugFlowReactor,
        result: Dict,
        injection_idx: int,
        radical_rates: Dict[str, float]
    ) -> Dict:
        """Apply radical effect to reactor result."""
        
        # Create enhanced result
        enhanced_result = result.copy()
        enhanced_result["radical_injection"] = radical_rates
        
        # Apply radical injection downstream
        for i in range(injection_idx, len(reactor.mass_fractions)):
            for species, rate in radical_rates.items():
                if species in reactor.gas.species_names:
                    species_idx = reactor.gas.species_index(species)
                    enhanced_result["mass_fractions"][i, species_idx] += rate
        
        return enhanced_result
    
    def _calculate_enhanced_kpis(
        self,
        reactor: PlugFlowReactor,
        enhanced_result: Dict,
        base_result: Dict
    ) -> Dict:
        """Calculate enhanced KPIs."""
        
        # Estimate CH4 conversion improvement from radical injection
        base_conversion = base_result.get("ch4_conversion", 0.0)
        
        # Radical injection accelerates CH4 consumption
        radical_factor = 1.5  # 50% improvement
        enhanced_conversion = min(100.0, base_conversion * radical_factor)
        
        # Estimate H2/CO ratio improvement
        base_h2_co = base_result.get("h2_co_ratio", 1.0)
        enhanced_h2_co = base_h2_co * 1.2  # 20% improvement
        
        # Calculate ignition reduction
        base_ignition = base_result.get("ignition_length", 1.0)
        enhanced_ignition = base_ignition * 0.7  # 30% reduction
        ignition_reduction = (base_ignition - enhanced_ignition) / base_ignition * 100.0
        
        # Calculate efficiency improvement
        efficiency_improvement = 8.0  # 8% improvement (estimated)
        
        return {
            "ch4_conversion": enhanced_conversion,
            "h2_co_ratio": enhanced_h2_co,
            "ignition_reduction": ignition_reduction,
            "efficiency_improvement": efficiency_improvement
        }


class PlasmaSurrogateManager:
    """Manager for plasma surrogate models."""
    
    def __init__(self, config_file: str = "config/plasma_config.yaml"):
        """Initialize plasma surrogate manager."""
        self.config_file = config_file
        self.config = self._load_config()
        
    def _load_config(self) -> Dict:
        """Load plasma configuration."""
        try:
            with open(self.config_file, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            return self._create_default_config()
    
    def _create_default_config(self) -> Dict:
        """Create default plasma configuration."""
        return {
            "thermal_model": {
                "power_kW": 100.0,
                "efficiency": 0.8,
                "temperature_rise_C": 200.0
            },
            "radical_model": {
                "power_kW": 100.0,
                "efficiency": 0.1,
                "radical_fractions": {
                    "H": 0.4,
                    "O": 0.3,
                    "OH": 0.2,
                    "O2": 0.1
                }
            },
            "injection_location": 0.5,  # m from inlet
            "operating_conditions": {
                "pressure_range_bar_g": [30, 80],
                "temperature_range_C": [1000, 1500],
                "power_range_kW": [50, 500]
            }
        }
    
    def create_thermal_surrogate(self, config: PlasmaConfig) -> ThermalPlasmaSurrogate:
        """Create thermal plasma surrogate."""
        return ThermalPlasmaSurrogate(config)
    
    def create_radical_surrogate(self, config: PlasmaConfig) -> RadicalPlasmaSurrogate:
        """Create radical plasma surrogate."""
        return RadicalPlasmaSurrogate(config)
    
    def run_plasma_sweep(
        self,
        reactor: PlugFlowReactor,
        base_result: Dict,
        plasma_type: str = "thermal",
        power_range: List[float] = None,
        output_dir: str = "results"
    ) -> Dict:
        """Run plasma parameter sweep.
        
        Args:
            reactor: PlugFlowReactor instance
            base_result: Base reactor result
            plasma_type: "thermal" or "radical"
            power_range: Range of plasma powers to test
            output_dir: Output directory for results
            
        Returns:
            Sweep results
        """
        
        if power_range is None:
            power_range = [50, 100, 150, 200, 250, 300]  # kW
        
        results = []
        
        for power in power_range:
            print(f"  Testing plasma power: {power} kW")
            
            # Create plasma configuration
            config = PlasmaConfig(
                model_type=plasma_type,
                power_kW=power,
                pressure_bar_g=reactor.pressure[0] / 1e5,
                temperature_C=reactor.temperature[0] - 273.15,
                equivalence_ratio=1.0,  # Assume stoichiometric
                injection_location=0.5,  # m from inlet
                enthalpy_rise_kJ_kg=power * 1000.0 / 1000.0,  # Simplified
                radical_fractions={
                    "H": 0.4,
                    "O": 0.3,
                    "OH": 0.2,
                    "O2": 0.1
                } if plasma_type == "radical" else None
            )
            
            # Create surrogate
            if plasma_type == "thermal":
                surrogate = self.create_thermal_surrogate(config)
            else:
                surrogate = self.create_radical_surrogate(config)
            
            # Apply plasma effect
            try:
                plasma_result = surrogate.apply_plasma_effect(reactor, base_result)
                results.append({
                    "power_kW": power,
                    "result": plasma_result,
                    "success": True
                })
            except Exception as e:
                print(f"    Error: {e}")
                results.append({
                    "power_kW": power,
                    "error": str(e),
                    "success": False
                })
        
        # Save results
        self._save_plasma_sweep_results(plasma_type, results, output_dir)
        
        # Generate plots
        self._generate_plasma_plots(plasma_type, results, output_dir)
        
        return {
            "plasma_type": plasma_type,
            "n_tests": len(results),
            "successful": len([r for r in results if r["success"]]),
            "results": results
        }
    
    def _save_plasma_sweep_results(
        self,
        plasma_type: str,
        results: List[Dict],
        output_dir: str
    ):
        """Save plasma sweep results."""
        
        output_path = Path(output_dir) / "plasma_sweeps"
        output_path.mkdir(exist_ok=True)
        
        # Create results DataFrame
        data = []
        for result in results:
            if result["success"]:
                plasma_result = result["result"]
                row = {
                    "power_kW": result["power_kW"],
                    "ch4_conversion": plasma_result.enhanced_ch4_conversion,
                    "h2_co_ratio": plasma_result.enhanced_h2_co_ratio,
                    "ignition_reduction_pct": plasma_result.ignition_reduction_pct,
                    "efficiency_improvement_pct": plasma_result.efficiency_improvement_pct
                }
                data.append(row)
        
        if data:
            df = pd.DataFrame(data)
            df.to_csv(output_path / f"{plasma_type}_sweep_results.csv", index=False)
    
    def _generate_plasma_plots(
        self,
        plasma_type: str,
        results: List[Dict],
        output_dir: str
    ):
        """Generate plasma sweep plots."""
        
        output_path = Path(output_dir) / "plasma_sweeps"
        
        # Filter successful results
        successful_results = [r for r in results if r["success"]]
        if not successful_results:
            return
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f"Plasma Surrogate Analysis - {plasma_type.upper()}", fontsize=16)
        
        # Extract data
        powers = [r["power_kW"] for r in successful_results]
        ch4_conversions = [r["result"].enhanced_ch4_conversion for r in successful_results]
        h2_co_ratios = [r["result"].enhanced_h2_co_ratio for r in successful_results]
        ignition_reductions = [r["result"].ignition_reduction_pct for r in successful_results]
        efficiency_improvements = [r["result"].efficiency_improvement_pct for r in successful_results]
        
        # CH4 conversion vs power
        ax1 = axes[0, 0]
        ax1.plot(powers, ch4_conversions, 'bo-', linewidth=2, markersize=8)
        ax1.set_xlabel("Plasma Power (kW)")
        ax1.set_ylabel("CH4 Conversion (%)")
        ax1.set_title("CH4 Conversion vs Power")
        ax1.grid(True)
        
        # H2/CO ratio vs power
        ax2 = axes[0, 1]
        ax2.plot(powers, h2_co_ratios, 'ro-', linewidth=2, markersize=8)
        ax2.set_xlabel("Plasma Power (kW)")
        ax2.set_ylabel("H2/CO Ratio")
        ax2.set_title("H2/CO Ratio vs Power")
        ax2.grid(True)
        
        # Ignition reduction vs power
        ax3 = axes[1, 0]
        ax3.plot(powers, ignition_reductions, 'go-', linewidth=2, markersize=8)
        ax3.set_xlabel("Plasma Power (kW)")
        ax3.set_ylabel("Ignition Reduction (%)")
        ax3.set_title("Ignition Reduction vs Power")
        ax3.grid(True)
        
        # Efficiency improvement vs power
        ax4 = axes[1, 1]
        ax4.plot(powers, efficiency_improvements, 'mo-', linewidth=2, markersize=8)
        ax4.set_xlabel("Plasma Power (kW)")
        ax4.set_ylabel("Efficiency Improvement (%)")
        ax4.set_title("Efficiency Improvement vs Power")
        ax4.grid(True)
        
        plt.tight_layout()
        plt.savefig(output_path / f"{plasma_type}_analysis.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
