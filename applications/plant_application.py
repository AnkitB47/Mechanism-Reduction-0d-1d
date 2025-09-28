"""Plant Application Framework.

Extends the HP-POX validation framework to industrial plant geometries A & B.
Supports operating envelope sweeps and plant-specific validation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import yaml
from dataclasses import dataclass

from reactor.plug_flow import PlugFlowReactor, ReactorGeometry, HeatTransferConfig, InletConditions
from validation.hp_pox_validator import HP_POXValidator


@dataclass
class OperatingEnvelope:
    """Operating envelope definition."""
    
    pressure_range_bar_g: List[float]
    temperature_range_C: List[float]
    flow_range_kgph: List[float]
    o2_c_ratio_range: List[float]
    diluent_fractions: Dict[str, float]  # N2, CO2, H2O


@dataclass
class PlantValidationData:
    """Plant validation data structure."""
    
    case_name: str
    pressure_bar_g: float
    temperature_C: float
    flow_rate_kgph: float
    o2_c_ratio: float
    
    # Measured plant data
    ch4_conversion_pct: float
    co2_conversion_pct: float
    h2_co_ratio: float
    outlet_composition: Dict[str, float]
    pressure_drop_bar: float
    ignition_length_m: float
    
    # Operating conditions
    ng_composition: Dict[str, float]
    inlet_temperatures: Dict[str, float]


class PlantApplication:
    """Plant application framework for geometries A & B."""
    
    def __init__(self, config_file: str = "config/plant_geometries.yaml"):
        """Initialize plant application with geometry configuration."""
        self.config_file = config_file
        self.config = self._load_config()
        self.results = {}
        
    def _load_config(self) -> Dict:
        """Load plant geometry configuration."""
        with open(self.config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def create_reactor(
        self,
        geometry_name: str,
        mechanism_file: str,
        n_points: int = 1000
    ) -> PlugFlowReactor:
        """Create reactor for plant geometry.
        
        Args:
            geometry_name: "plant_a" or "plant_b"
            mechanism_file: Path to mechanism file
            n_points: Number of axial discretization points
            
        Returns:
            Configured PlugFlowReactor
        """
        
        # Load geometry
        geometry = self._load_geometry(geometry_name)
        
        # Create heat transfer configuration
        heat_transfer = self._create_heat_transfer_config(geometry_name)
        
        # Create reactor
        reactor = PlugFlowReactor(
            mechanism_file=mechanism_file,
            geometry=geometry,
            heat_transfer=heat_transfer,
            n_points=n_points
        )
        
        return reactor
    
    def _load_geometry(self, geometry_name: str) -> ReactorGeometry:
        """Load geometry configuration."""
        geom_data = self.config["geometries"][geometry_name]
        
        return ReactorGeometry(
            name=geom_data["name"],
            total_length=geom_data["total_length"],
            diameter_constant=geom_data["diameter_constant"],
            diameter_narrow=geom_data["diameter_narrow"],
            narrow_start_height=geom_data["narrow_start_height"],
            segments=geom_data["segments"]
        )
    
    def _create_heat_transfer_config(self, geometry_name: str) -> HeatTransferConfig:
        """Create heat transfer configuration for plant geometry."""
        geom_data = self.config["geometries"][geometry_name]
        ht_data = geom_data["heat_transfer"]
        
        return HeatTransferConfig(
            model="fixed_wall_temp",
            wall_temp_C=ht_data["wall_temp_C"],
            heat_loss_fraction=ht_data["heat_loss_fraction"]
        )
    
    def run_operating_envelope_sweep(
        self,
        geometry_name: str,
        mechanism_file: str,
        envelope: OperatingEnvelope,
        n_points_per_dim: int = 5,
        output_dir: str = "results"
    ) -> Dict:
        """Run operating envelope sweep.
        
        Args:
            geometry_name: "plant_a" or "plant_b"
            mechanism_file: Path to mechanism file
            envelope: Operating envelope definition
            n_points_per_dim: Number of points per dimension
            output_dir: Output directory for results
            
        Returns:
            Sweep results dictionary
        """
        
        print(f"Running operating envelope sweep for {geometry_name}...")
        
        # Create reactor
        reactor = self.create_reactor(geometry_name, mechanism_file)
        
        # Generate parameter grid
        param_grid = self._generate_parameter_grid(envelope, n_points_per_dim)
        
        # Run simulations
        results = []
        for i, params in enumerate(param_grid):
            print(f"  Simulation {i+1}/{len(param_grid)}: {params}")
            
            try:
                # Create inlet conditions
                inlet = self._create_inlet_conditions(params)
                
                # Set inlet conditions
                reactor.set_inlet_conditions(inlet)
                
                # Solve reactor
                result = reactor.solve()
                
                # Store results
                results.append({
                    "parameters": params,
                    "result": result,
                    "kpis": self._calculate_kpis(result)
                })
                
            except Exception as e:
                print(f"    Error: {e}")
                results.append({
                    "parameters": params,
                    "error": str(e)
                })
        
        # Save results
        self._save_sweep_results(geometry_name, results, output_dir)
        
        # Generate response surfaces
        self._generate_response_surfaces(geometry_name, results, output_dir)
        
        return {
            "geometry": geometry_name,
            "envelope": envelope,
            "n_simulations": len(results),
            "successful": len([r for r in results if "error" not in r]),
            "results": results
        }
    
    def _generate_parameter_grid(self, envelope: OperatingEnvelope, n_points: int) -> List[Dict]:
        """Generate parameter grid for operating envelope sweep."""
        
        # Create parameter ranges
        pressure_range = np.linspace(
            envelope.pressure_range_bar_g[0],
            envelope.pressure_range_bar_g[1],
            n_points
        )
        
        temperature_range = np.linspace(
            envelope.temperature_range_C[0],
            envelope.temperature_range_C[1],
            n_points
        )
        
        flow_range = np.linspace(
            envelope.flow_range_kgph[0],
            envelope.flow_range_kgph[1],
            n_points
        )
        
        o2_c_range = np.linspace(
            envelope.o2_c_ratio_range[0],
            envelope.o2_c_ratio_range[1],
            n_points
        )
        
        # Generate grid
        param_grid = []
        for p in pressure_range:
            for T in temperature_range:
                for flow in flow_range:
                    for o2_c in o2_c_range:
                        param_grid.append({
                            "pressure_bar_g": p,
                            "temperature_C": T,
                            "flow_rate_kgph": flow,
                            "o2_c_ratio": o2_c,
                            "diluents": envelope.diluent_fractions
                        })
        
        return param_grid
    
    def _create_inlet_conditions(self, params: Dict) -> InletConditions:
        """Create inlet conditions from parameters."""
        
        # Calculate flow rates based on O2/C ratio
        flow_rate_kgph = params["flow_rate_kgph"]
        o2_c_ratio = params["o2_c_ratio"]
        
        # Assume natural gas is mostly CH4
        ng_flow = flow_rate_kgph * 0.7  # 70% natural gas
        o2_flow = ng_flow * o2_c_ratio * 4.0  # Stoichiometric O2 for CH4
        
        # Calculate diluent flows
        diluents = params["diluents"]
        h2o_flow = flow_rate_kgph * diluents.get("H2O", 0.1)
        n2_flow = flow_rate_kgph * diluents.get("N2", 0.05)
        
        # Use typical natural gas composition
        ng_composition = self.config["ng_compositions"]["typical"]
        
        return InletConditions(
            natural_gas_kgph=ng_flow,
            oxygen_kgph=o2_flow,
            primary_steam_kgph=h2o_flow * 0.7,
            secondary_steam_kgph=h2o_flow * 0.3,
            n2_optical_kgph=n2_flow,
            natural_gas_T_C=params["temperature_C"] - 50,  # Preheat
            oxygen_T_C=params["temperature_C"] - 100,
            primary_steam_T_C=params["temperature_C"] - 200,
            secondary_steam_T_C=params["temperature_C"] - 200,
            n2_T_C=25.0,
            pressure_bar_g=params["pressure_bar_g"],
            ng_composition=ng_composition
        )
    
    def _calculate_kpis(self, result) -> Dict:
        """Calculate key performance indicators."""
        
        return {
            "ch4_conversion": result.ch4_conversion,
            "co2_conversion": result.co2_conversion,
            "h2_co_ratio": result.h2_co_ratio,
            "ignition_length": result.ignition_length,
            "pressure_drop": (result.pressure[0] - result.pressure[-1]) / 1e5,
            "outlet_temp": result.temperature[-1] - 273.15,
            "max_temp": np.max(result.temperature) - 273.15
        }
    
    def _save_sweep_results(self, geometry_name: str, results: List[Dict], output_dir: str):
        """Save operating envelope sweep results."""
        
        output_path = Path(output_dir) / "plant_sweeps"
        output_path.mkdir(exist_ok=True)
        
        # Create results DataFrame
        data = []
        for result in results:
            if "error" not in result:
                row = result["parameters"].copy()
                row.update(result["kpis"])
                data.append(row)
        
        if data:
            df = pd.DataFrame(data)
            df.to_csv(output_path / f"{geometry_name}_sweep_results.csv", index=False)
        
        # Save summary
        summary = {
            "geometry": geometry_name,
            "total_simulations": len(results),
            "successful": len([r for r in results if "error" not in r]),
            "failed": len([r for r in results if "error" in r])
        }
        
        with open(output_path / f"{geometry_name}_summary.yaml", 'w') as f:
            yaml.dump(summary, f, default_flow_style=False)
    
    def _generate_response_surfaces(self, geometry_name: str, results: List[Dict], output_dir: str):
        """Generate response surface plots."""
        
        output_path = Path(output_dir) / "plant_sweeps"
        
        # Filter successful results
        successful_results = [r for r in results if "error" not in r]
        if not successful_results:
            return
        
        # Create DataFrame
        data = []
        for result in successful_results:
            row = result["parameters"].copy()
            row.update(result["kpis"])
            data.append(row)
        
        df = pd.DataFrame(data)
        
        # Generate response surface plots
        self._plot_response_surfaces(geometry_name, df, output_path)
    
    def _plot_response_surfaces(self, geometry_name: str, df: pd.DataFrame, output_path: Path):
        """Plot response surfaces for key parameters."""
        
        # CH4 conversion vs pressure and temperature
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f"Response Surfaces - {geometry_name.upper()}", fontsize=16)
        
        # CH4 conversion vs pressure and temperature
        ax1 = axes[0, 0]
        scatter = ax1.scatter(df["pressure_bar_g"], df["temperature_C"], 
                             c=df["ch4_conversion"], cmap="viridis", s=50)
        ax1.set_xlabel("Pressure (bar g)")
        ax1.set_ylabel("Temperature (°C)")
        ax1.set_title("CH4 Conversion (%)")
        plt.colorbar(scatter, ax=ax1)
        
        # H2/CO ratio vs pressure and temperature
        ax2 = axes[0, 1]
        scatter = ax2.scatter(df["pressure_bar_g"], df["temperature_C"], 
                             c=df["h2_co_ratio"], cmap="plasma", s=50)
        ax2.set_xlabel("Pressure (bar g)")
        ax2.set_ylabel("Temperature (°C)")
        ax2.set_title("H2/CO Ratio")
        plt.colorbar(scatter, ax=ax2)
        
        # Ignition length vs flow rate and O2/C ratio
        ax3 = axes[1, 0]
        scatter = ax3.scatter(df["flow_rate_kgph"], df["o2_c_ratio"], 
                             c=df["ignition_length"], cmap="coolwarm", s=50)
        ax3.set_xlabel("Flow Rate (kg/h)")
        ax3.set_ylabel("O2/C Ratio")
        ax3.set_title("Ignition Length (m)")
        plt.colorbar(scatter, ax=ax3)
        
        # Pressure drop vs flow rate and pressure
        ax4 = axes[1, 1]
        scatter = ax4.scatter(df["flow_rate_kgph"], df["pressure_bar_g"], 
                             c=df["pressure_drop"], cmap="Reds", s=50)
        ax4.set_xlabel("Flow Rate (kg/h)")
        ax4.set_ylabel("Pressure (bar g)")
        ax4.set_title("Pressure Drop (bar)")
        plt.colorbar(scatter, ax=ax4)
        
        plt.tight_layout()
        plt.savefig(output_path / f"{geometry_name}_response_surfaces.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    def validate_against_plant_data(
        self,
        geometry_name: str,
        mechanism_file: str,
        plant_data: List[PlantValidationData],
        output_dir: str = "results"
    ) -> Dict:
        """Validate against plant measurement data.
        
        Args:
            geometry_name: "plant_a" or "plant_b"
            mechanism_file: Path to mechanism file
            plant_data: List of plant validation data points
            output_dir: Output directory for results
            
        Returns:
            Validation results
        """
        
        print(f"Validating {geometry_name} against plant data...")
        
        # Create reactor
        reactor = self.create_reactor(geometry_name, mechanism_file)
        
        validation_results = []
        
        for data_point in plant_data:
            print(f"  Validating case: {data_point.case_name}")
            
            try:
                # Create inlet conditions
                inlet = InletConditions(
                    natural_gas_kgph=data_point.flow_rate_kgph * 0.7,
                    oxygen_kgph=data_point.flow_rate_kgph * 0.7 * data_point.o2_c_ratio * 4.0,
                    primary_steam_kgph=data_point.flow_rate_kgph * 0.1,
                    secondary_steam_kgph=data_point.flow_rate_kgph * 0.05,
                    n2_optical_kgph=data_point.flow_rate_kgph * 0.05,
                    natural_gas_T_C=data_point.inlet_temperatures.get("ng", data_point.temperature_C - 50),
                    oxygen_T_C=data_point.inlet_temperatures.get("o2", data_point.temperature_C - 100),
                    primary_steam_T_C=data_point.inlet_temperatures.get("steam", data_point.temperature_C - 200),
                    secondary_steam_T_C=data_point.inlet_temperatures.get("steam", data_point.temperature_C - 200),
                    n2_T_C=25.0,
                    pressure_bar_g=data_point.pressure_bar_g,
                    ng_composition=data_point.ng_composition
                )
                
                # Set inlet conditions
                reactor.set_inlet_conditions(inlet)
                
                # Solve reactor
                result = reactor.solve()
                
                # Calculate validation metrics
                validation = self._calculate_validation_metrics(data_point, result)
                validation_results.append(validation)
                
            except Exception as e:
                print(f"    Error: {e}")
                validation_results.append({
                    "case": data_point.case_name,
                    "error": str(e),
                    "passed": False
                })
        
        # Save validation results
        self._save_plant_validation_results(geometry_name, validation_results, output_dir)
        
        return {
            "geometry": geometry_name,
            "n_cases": len(plant_data),
            "results": validation_results
        }
    
    def _calculate_validation_metrics(
        self,
        plant_data: PlantValidationData,
        result
    ) -> Dict:
        """Calculate validation metrics for plant data."""
        
        validation = {
            "case": plant_data.case_name,
            "passed": True,
            "errors": {},
            "metrics": {}
        }
        
        # CH4 conversion validation
        ch4_error = abs(result.ch4_conversion - plant_data.ch4_conversion_pct)
        if ch4_error > 5.0:  # 5% tolerance
            validation["errors"]["ch4_conversion"] = {
                "expected": plant_data.ch4_conversion_pct,
                "simulated": result.ch4_conversion,
                "error": ch4_error
            }
            validation["passed"] = False
        
        # H2/CO ratio validation
        h2_co_error = abs(result.h2_co_ratio - plant_data.h2_co_ratio)
        if h2_co_error > 0.1:  # 0.1 tolerance
            validation["errors"]["h2_co_ratio"] = {
                "expected": plant_data.h2_co_ratio,
                "simulated": result.h2_co_ratio,
                "error": h2_co_error
            }
            validation["passed"] = False
        
        # Pressure drop validation
        pressure_drop = (result.pressure[0] - result.pressure[-1]) / 1e5
        pressure_error = abs(pressure_drop - plant_data.pressure_drop_bar)
        if pressure_error > 0.1:  # 0.1 bar tolerance
            validation["warnings"]["pressure_drop"] = {
                "expected": plant_data.pressure_drop_bar,
                "simulated": pressure_drop,
                "error": pressure_error
            }
        
        validation["metrics"] = {
            "ch4_conversion": result.ch4_conversion,
            "h2_co_ratio": result.h2_co_ratio,
            "pressure_drop": pressure_drop,
            "ignition_length": result.ignition_length
        }
        
        return validation
    
    def _save_plant_validation_results(
        self,
        geometry_name: str,
        validation_results: List[Dict],
        output_dir: str
    ):
        """Save plant validation results."""
        
        output_path = Path(output_dir) / "plant_validation"
        output_path.mkdir(exist_ok=True)
        
        # Save validation summary
        summary = {
            "geometry": geometry_name,
            "total_cases": len(validation_results),
            "passed": len([r for r in validation_results if r.get("passed", False)]),
            "failed": len([r for r in validation_results if not r.get("passed", True)]),
            "results": validation_results
        }
        
        with open(output_path / f"{geometry_name}_validation.yaml", 'w') as f:
            yaml.dump(summary, f, default_flow_style=False)
        
        # Save detailed results CSV
        data = []
        for result in validation_results:
            if "error" not in result:
                row = {"case": result["case"]}
                row.update(result["metrics"])
                data.append(row)
        
        if data:
            df = pd.DataFrame(data)
            df.to_csv(output_path / f"{geometry_name}_validation.csv", index=False)
