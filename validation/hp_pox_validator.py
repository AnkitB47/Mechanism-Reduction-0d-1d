"""HP-POX Benchmark Validation System.

Validates 1-D reactor models against Richter et al. (Fuel 2015) benchmark data.
Supports Cases 1-4 with comprehensive KPI validation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import yaml

from reactor.plug_flow import PlugFlowReactor, ReactorGeometry, HeatTransferConfig, InletConditions
from reactor.batch import BatchResult


class HP_POXValidator:
    """HP-POX benchmark validation system."""
    
    def __init__(self, config_file: str = "config/hp_pox_cases.yaml"):
        """Initialize validator with HP-POX configuration."""
        self.config_file = config_file
        self.config = self._load_config()
        self.results = {}
        
    def _load_config(self) -> Dict:
        """Load HP-POX configuration."""
        with open(self.config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def validate_case(
        self,
        case_name: str,
        mechanism_file: str,
        n_points: int = 1000,
        output_dir: str = "results"
    ) -> Dict:
        """Validate a specific HP-POX case.
        
        Args:
            case_name: Case name (e.g., "case_1", "case_4")
            mechanism_file: Path to mechanism file
            n_points: Number of axial discretization points
            output_dir: Output directory for results
            
        Returns:
            Validation results dictionary
        """
        
        print(f"Validating {case_name}...")
        
        # Load case data
        case_data = self.config["cases"][case_name]
        
        # Create reactor geometry
        geometry = self._create_geometry()
        
        # Create heat transfer configuration
        heat_transfer = self._create_heat_transfer_config(case_name)
        
        # Create inlet conditions
        inlet = self._create_inlet_conditions(case_data)
        
        # Initialize reactor
        reactor = PlugFlowReactor(
            mechanism_file=mechanism_file,
            geometry=geometry,
            heat_transfer=heat_transfer,
            n_points=n_points
        )
        
        # Set inlet conditions
        reactor.set_inlet_conditions(inlet)
        
        # Solve reactor
        result = reactor.solve()
        
        # Validate against experimental data
        validation_results = self._validate_results(case_name, result, case_data)
        
        # Save results
        self._save_results(case_name, result, validation_results, output_dir)
        
        # Generate plots
        self._generate_plots(case_name, result, case_data, output_dir)
        
        self.results[case_name] = {
            "result": result,
            "validation": validation_results
        }
        
        return validation_results
    
    def _create_geometry(self) -> ReactorGeometry:
        """Create reactor geometry from config."""
        geom_data = self.config["geometry"]
        
        return ReactorGeometry(
            name=geom_data["name"],
            total_length=geom_data["total_length"],
            diameter_constant=geom_data["diameter_constant"],
            diameter_narrow=geom_data["diameter_narrow"],
            narrow_start_height=geom_data["narrow_start_height"],
            segments=geom_data["segments"]
        )
    
    def _create_heat_transfer_config(self, case_name: str) -> HeatTransferConfig:
        """Create heat transfer configuration from case data."""
        case_data = self.config["cases"][case_name]
        
        # Use wall temperature from measurements
        wall_temps = case_data["wall_temps"]
        avg_wall_temp = np.mean(list(wall_temps.values()))
        
        return HeatTransferConfig(
            model="fixed_wall_temp",
            wall_temp_C=avg_wall_temp,
            heat_loss_fraction=0.01
        )
    
    def _create_inlet_conditions(self, case_data: Dict) -> InletConditions:
        """Create inlet conditions from case data."""
        
        return InletConditions(
            natural_gas_kgph=case_data["natural_gas_kgph"],
            oxygen_kgph=case_data["oxygen_kgph"],
            primary_steam_kgph=case_data["primary_steam_kgph"],
            secondary_steam_kgph=case_data["secondary_steam_kgph"],
            n2_optical_kgph=case_data["n2_optical_kgph"],
            natural_gas_T_C=case_data["natural_gas_T_C"],
            oxygen_T_C=case_data["oxygen_T_C"],
            primary_steam_T_C=case_data["primary_steam_T_C"],
            secondary_steam_T_C=case_data["secondary_steam_T_C"],
            n2_T_C=case_data["n2_T_C"],
            pressure_bar_g=case_data["pressure_bar_g"],
            ng_composition=case_data["ng_composition"]
        )
    
    def _validate_results(
        self,
        case_name: str,
        result,
        case_data: Dict
    ) -> Dict:
        """Validate simulation results against experimental data."""
        
        validation = {
            "case": case_name,
            "passed": True,
            "errors": {},
            "warnings": {},
            "metrics": {}
        }
        
        # Temperature validation
        self._validate_temperature(result, case_data, validation)
        
        # Composition validation
        self._validate_composition(result, case_data, validation)
        
        # Flame length validation (if available)
        self._validate_flame_length(result, case_data, validation)
        
        # Pressure validation
        self._validate_pressure(result, case_data, validation)
        
        return validation
    
    def _validate_temperature(self, result, case_data: Dict, validation: Dict):
        """Validate temperature profiles."""
        
        # Get experimental temperature data
        exp_temps = case_data["axial_temps"]
        
        # Calculate simulated temperatures at measurement points
        sim_temps = {}
        for key, exp_temp in exp_temps.items():
            # Map measurement points to axial positions
            if "G1to4" in key:
                height = 0.5  # m
            elif "G5to8" in key:
                height = 1.0  # m
            elif "G9to11" in key:
                height = 1.5  # m
            elif "G12" in key:
                height = 2.0  # m
            elif "G13" in key:
                height = 2.5  # m
            elif "G14" in key:
                height = 3.0  # m
            else:
                continue
            
            # Interpolate simulated temperature
            sim_temp = np.interp(height, result.height, result.temperature)
            sim_temps[key] = sim_temp - 273.15  # Convert to °C
            
            # Calculate error
            error = abs(sim_temp - 273.15 - exp_temp)
            tolerance = self.config["validation"]["temperature_tolerance_C"]
            
            if error > tolerance:
                validation["errors"][f"temp_{key}"] = {
                    "expected": exp_temp,
                    "simulated": sim_temp - 273.15,
                    "error": error,
                    "tolerance": tolerance
                }
                validation["passed"] = False
        
        validation["metrics"]["temperature"] = sim_temps
    
    def _validate_composition(self, result, case_data: Dict, validation: Dict):
        """Validate outlet composition."""
        
        # Get expected composition
        exp_comp = case_data["expected_outlet"]
        
        # Get simulated composition
        sim_comp = result.outlet_composition
        
        # Validate each species
        for species, exp_volpct in exp_comp.items():
            if species in sim_comp:
                sim_volpct = sim_comp[species]
                error = abs(sim_volpct - exp_volpct)
                tolerance = self.config["validation"]["composition_tolerance_pct"]
                
                if error > tolerance:
                    validation["errors"][f"comp_{species}"] = {
                        "expected": exp_volpct,
                        "simulated": sim_volpct,
                        "error": error,
                        "tolerance": tolerance
                    }
                    validation["passed"] = False
            else:
                validation["warnings"][f"comp_{species}"] = "Species not found in simulation"
        
        validation["metrics"]["composition"] = sim_comp
    
    def _validate_flame_length(self, result, case_data: Dict, validation: Dict):
        """Validate flame length (if available)."""
        
        flame_data = case_data.get("flame", {})
        exp_length = flame_data.get("length_mm")
        
        if exp_length is not None:
            sim_length = result.ignition_length * 1000  # Convert to mm
            error = abs(sim_length - exp_length)
            tolerance = self.config["validation"]["flame_length_tolerance_mm"]
            
            if error > tolerance:
                validation["errors"]["flame_length"] = {
                    "expected": exp_length,
                    "simulated": sim_length,
                    "error": error,
                    "tolerance": tolerance
                }
                validation["passed"] = False
            
            validation["metrics"]["flame_length_mm"] = sim_length
    
    def _validate_pressure(self, result, case_data: Dict, validation: Dict):
        """Validate pressure drop."""
        
        # Calculate pressure drop
        inlet_pressure = result.pressure[0] / 1e5  # Convert to bar
        outlet_pressure = result.pressure[-1] / 1e5
        pressure_drop = inlet_pressure - outlet_pressure
        
        # Expected pressure drop (estimated)
        expected_drop = 0.1  # bar (estimated)
        error = abs(pressure_drop - expected_drop)
        tolerance = expected_drop * self.config["validation"]["pressure_drop_tolerance_pct"] / 100.0
        
        if error > tolerance:
            validation["warnings"]["pressure_drop"] = {
                "expected": expected_drop,
                "simulated": pressure_drop,
                "error": error,
                "tolerance": tolerance
            }
        
        validation["metrics"]["pressure_drop_bar"] = pressure_drop
    
    def _save_results(self, case_name: str, result, validation: Dict, output_dir: str):
        """Save validation results."""
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Save axial profiles
        profiles_df = pd.DataFrame({
            "height_m": result.height,
            "temperature_K": result.temperature,
            "pressure_Pa": result.pressure,
            "density_kg_m3": result.density,
            "velocity_m_s": result.velocity,
            "residence_time_s": result.time
        })
        
        # Add species mass fractions
        for i, species in enumerate(result.species_names):
            profiles_df[f"Y_{species}"] = result.mass_fractions[:, i]
        
        profiles_df.to_csv(output_path / f"{case_name}_profiles.csv", index=False)
        
        # Save validation summary
        validation_summary = {
            "case": case_name,
            "passed": validation["passed"],
            "ch4_conversion": result.ch4_conversion,
            "co2_conversion": result.co2_conversion,
            "h2_co_ratio": result.h2_co_ratio,
            "ignition_length_m": result.ignition_length,
            "max_dT_dx_K_m": result.max_dT_dx,
            "errors": validation["errors"],
            "warnings": validation["warnings"]
        }
        
        with open(output_path / f"{case_name}_validation.yaml", 'w') as f:
            yaml.dump(validation_summary, f, default_flow_style=False)
    
    def _generate_plots(self, case_name: str, result, case_data: Dict, output_dir: str):
        """Generate validation plots."""
        
        output_path = Path(output_dir)
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"HP-POX Validation - {case_name.upper()}", fontsize=16)
        
        # Temperature profile
        ax1 = axes[0, 0]
        ax1.plot(result.height, result.temperature - 273.15, 'b-', linewidth=2, label='Simulation')
        
        # Add experimental temperature points
        exp_temps = case_data["axial_temps"]
        heights = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        exp_values = [exp_temps["G1to4_mean_C"], exp_temps["G5to8_mean_C"], 
                     exp_temps["G9to11_mean_C"], exp_temps["G12_control_C"],
                     exp_temps["G13_center_C"], exp_temps["G14_wall_C"]]
        ax1.plot(heights, exp_values, 'ro', markersize=8, label='Experimental')
        
        ax1.set_xlabel('Height (m)')
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_title('Axial Temperature Profile')
        ax1.legend()
        ax1.grid(True)
        
        # Species profiles
        ax2 = axes[0, 1]
        key_species = ["CH4", "H2", "CO", "CO2", "H2O"]
        colors = ['b', 'r', 'g', 'm', 'c']
        
        for i, species in enumerate(key_species):
            if species in result.species_names:
                idx = result.species_names.index(species)
                ax2.plot(result.height, result.mass_fractions[:, idx] * 100, 
                        color=colors[i], linewidth=2, label=species)
        
        ax2.set_xlabel('Height (m)')
        ax2.set_ylabel('Mass Fraction (%)')
        ax2.set_title('Species Profiles')
        ax2.legend()
        ax2.grid(True)
        
        # Pressure profile
        ax3 = axes[1, 0]
        ax3.plot(result.height, result.pressure / 1e5, 'b-', linewidth=2)
        ax3.set_xlabel('Height (m)')
        ax3.set_ylabel('Pressure (bar)')
        ax3.set_title('Axial Pressure Profile')
        ax3.grid(True)
        
        # Velocity profile
        ax4 = axes[1, 1]
        ax4.plot(result.height, result.velocity, 'b-', linewidth=2)
        ax4.set_xlabel('Height (m)')
        ax4.set_ylabel('Velocity (m/s)')
        ax4.set_title('Axial Velocity Profile')
        ax4.grid(True)
        
        plt.tight_layout()
        plt.savefig(output_path / f"{case_name}_profiles.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Composition comparison plot
        self._plot_composition_comparison(case_name, result, case_data, output_path)
    
    def _plot_composition_comparison(self, case_name: str, result, case_data: Dict, output_path: Path):
        """Plot composition comparison between simulation and experiment."""
        
        exp_comp = case_data["expected_outlet"]
        sim_comp = result.outlet_composition
        
        # Prepare data for comparison
        species = list(exp_comp.keys())
        exp_values = [exp_comp[s] for s in species]
        sim_values = [sim_comp.get(s, 0) for s in species]
        
        # Create bar plot
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(species))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, exp_values, width, label='Experimental', alpha=0.8)
        bars2 = ax.bar(x + width/2, sim_values, width, label='Simulation', alpha=0.8)
        
        ax.set_xlabel('Species')
        ax.set_ylabel('Volume Fraction (%)')
        ax.set_title(f'Outlet Composition Comparison - {case_name.upper()}')
        ax.set_xticks(x)
        ax.set_xticklabels(species)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{height:.1f}', ha='center', va='bottom')
        
        for bar in bars2:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{height:.1f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(output_path / f"{case_name}_composition.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    def validate_all_cases(
        self,
        mechanism_file: str,
        cases: List[str] = None,
        output_dir: str = "results"
    ) -> Dict:
        """Validate all HP-POX cases.
        
        Args:
            mechanism_file: Path to mechanism file
            cases: List of cases to validate (default: ["case_1", "case_4"])
            output_dir: Output directory for results
            
        Returns:
            Summary of all validation results
        """
        
        if cases is None:
            cases = ["case_1", "case_4"]  # Focus on baseline and high-T cases
        
        summary = {
            "mechanism": mechanism_file,
            "cases_validated": cases,
            "overall_passed": True,
            "case_results": {}
        }
        
        for case in cases:
            try:
                validation = self.validate_case(case, mechanism_file, output_dir=output_dir)
                summary["case_results"][case] = validation
                
                if not validation["passed"]:
                    summary["overall_passed"] = False
                    
            except Exception as e:
                summary["case_results"][case] = {
                    "passed": False,
                    "error": str(e)
                }
                summary["overall_passed"] = False
        
        # Save summary
        summary_path = Path(output_dir) / "validation_summary.yaml"
        with open(summary_path, 'w') as f:
            yaml.dump(summary, f, default_flow_style=False)
        
        return summary
