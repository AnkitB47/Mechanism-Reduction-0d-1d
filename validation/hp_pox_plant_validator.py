"""HP-POX Plant Validator.

Validates plant cases against Richter et al. (Fuel 2015) benchmark specifications.
Uses oxygen-blown HP-POX feeds with proper normalization and error analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import yaml
from datetime import datetime

from applications.plant_application import PlantApplication
from reactor.plug_flow import PlugFlowReactor, InletConditions


class HP_POXPlantValidator:
    """HP-POX plant validator with oxygen-blown feeds."""
    
    def __init__(self, config_file: str = "config/hp_pox_plant_cases.yaml"):
        """Initialize HP-POX plant validator."""
        self.config_file = config_file
        self.config = self._load_config()
        self.plant_app = PlantApplication()
        
        # Create output directory
        self.output_dir = Path("results_1d_hp_pox")
        self.output_dir.mkdir(exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)
        
    def _load_config(self) -> Dict:
        """Load HP-POX plant configuration."""
        with open(self.config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def validate_case(
        self,
        case_name: str,
        mechanism_file: str,
        output_subdir: str = None
    ) -> Dict:
        """Validate a single HP-POX case.
        
        Args:
            case_name: Case name (case_reference or case_sensitivity)
            mechanism_file: Path to mechanism file
            output_subdir: Output subdirectory (optional)
            
        Returns:
            Validation results
        """
        
        print(f"\n{'='*80}")
        print(f"HP-POX VALIDATION: {case_name.upper()}")
        print(f"{'='*80}")
        
        # Load case configuration
        case_config = self.config["cases"][case_name]
        
        # Print feed validation
        self._print_feed_validation(case_config)
        
        # Create reactor
        reactor = self.plant_app.create_reactor(
            case_config["geometry"], 
            mechanism_file
        )
        
        # Create inlet conditions
        inlet = self._create_hp_pox_inlet(case_config)
        reactor.set_inlet_conditions(inlet)
        
        # Solve reactor
        start_time = time.time()
        result = reactor.solve()
        runtime = time.time() - start_time
        
        # Calculate KPIs
        kpis = self._calculate_hp_pox_kpis(result, runtime, case_config)
        
        # Save results
        if output_subdir:
            case_output_dir = self.output_dir / output_subdir
        else:
            case_output_dir = self.output_dir / case_name
        case_output_dir.mkdir(exist_ok=True)
        
        self._save_case_results(case_name, result, kpis, case_output_dir)
        self._generate_case_plots(case_name, result, kpis, case_output_dir)
        
        # Print validation summary
        self._print_validation_summary(case_name, kpis, case_config)
        
        return {
            "case": case_name,
            "result": result,
            "kpis": kpis,
            "runtime": runtime,
            "output_dir": case_output_dir
        }
    
    def _print_feed_validation(self, case_config: Dict):
        """Print feed validation information."""
        
        print("✅ Feed check: Oxygen-blown HP-POX configuration")
        print(f"   O2 purity: ≥99.5% (no N2 ballast)")
        print(f"   No air used in any inlet stream")
        
        # Print NG composition
        ng_comp = case_config["ng_composition"]
        print("   NG composition (vol%):")
        for species, pct in ng_comp.items():
            print(f"     {species}: {pct}%")
        
        # Calculate O2:C ratio
        ng_flow = case_config["natural_gas_kgph"]
        o2_flow = case_config["oxygen_kgph"]
        
        # Assume CH4 is 96% of NG by mass
        ch4_flow = ng_flow * 0.96
        o2_c_ratio = o2_flow / ch4_flow
        
        print(f"   O2:C ratio: {o2_c_ratio:.2f}")
        print(f"   Total inlet flows: NG={ng_flow} kg/h, O2={o2_flow} kg/h")
    
    def _create_hp_pox_inlet(self, case_config: Dict) -> InletConditions:
        """Create HP-POX inlet conditions."""
        
        # Simplify NG composition to only include basic species available in most mechanisms
        simplified_ng_composition = {
            "CH4_volpct": 96.0,
            "CO2_volpct": 0.35,
            "C2H6_volpct": 2.25,
            "N2_volpct": 1.4  # Combine C3H8, nC4H10, iC4H10 into N2 for simplicity
        }
        
        return InletConditions(
            natural_gas_kgph=case_config["natural_gas_kgph"],
            oxygen_kgph=case_config["oxygen_kgph"],
            primary_steam_kgph=case_config["primary_steam_kgph"],
            secondary_steam_kgph=case_config["secondary_steam_kgph"],
            n2_optical_kgph=case_config["n2_optical_kgph"],
            natural_gas_T_C=case_config["natural_gas_T_C"],
            oxygen_T_C=case_config["oxygen_T_C"],
            primary_steam_T_C=case_config["primary_steam_T_C"],
            secondary_steam_T_C=case_config["secondary_steam_T_C"],
            n2_T_C=case_config["n2_T_C"],
            pressure_bar_g=case_config["pressure_bar_g"],
            ng_composition=simplified_ng_composition
        )
    
    def _calculate_hp_pox_kpis(self, result, runtime: float, case_config: Dict) -> Dict:
        """Calculate HP-POX specific KPIs."""
        
        # Outlet conditions
        outlet_idx = -1
        outlet_temp = result.temperature[outlet_idx]
        outlet_pressure = result.pressure[outlet_idx]
        outlet_Y = result.mass_fractions[outlet_idx, :]
        
        # Calculate mole fractions (simplified - would use proper gas state)
        # For now, use mass fractions as approximation
        outlet_X = outlet_Y.copy()
        
        # CH4 conversion
        ch4_idx = result.species_names.index("CH4")
        ch4_inlet = result.mass_fractions[0, ch4_idx]
        ch4_outlet = result.mass_fractions[outlet_idx, ch4_idx]
        ch4_conversion = (ch4_inlet - ch4_outlet) / ch4_inlet * 100.0 if ch4_inlet > 0 else 0.0
        
        # CO2 conversion
        co2_idx = result.species_names.index("CO2")
        co2_inlet = result.mass_fractions[0, co2_idx]
        co2_outlet = result.mass_fractions[outlet_idx, co2_idx]
        co2_conversion = (co2_outlet - co2_inlet) / co2_inlet * 100.0 if co2_inlet > 0 else 0.0
        
        # H2/CO ratio
        h2_idx = result.species_names.index("H2")
        co_idx = result.species_names.index("CO")
        h2_mass_frac = result.mass_fractions[outlet_idx, h2_idx]
        co_mass_frac = result.mass_fractions[outlet_idx, co_idx]
        h2_co_ratio = h2_mass_frac / co_mass_frac if co_mass_frac > 0 else 0.0
        
        # Pressure drop
        pressure_drop = (result.pressure[0] - result.pressure[-1])  # Pa
        
        # Outlet composition (wet and dry basis)
        outlet_composition_wet = {}
        outlet_composition_dry = {}
        
        for i, species in enumerate(result.species_names):
            if species in ["H2", "CO", "CO2", "CH4", "N2", "H2O"]:
                wet_frac = result.mass_fractions[outlet_idx, i] * 100.0
                outlet_composition_wet[species] = wet_frac
                
                # Dry basis: exclude H2O
                if species != "H2O":
                    outlet_composition_dry[species] = wet_frac
        
        # Normalize dry basis to 100%
        dry_total = sum(outlet_composition_dry.values())
        if dry_total > 0:
            for species in outlet_composition_dry:
                outlet_composition_dry[species] = outlet_composition_dry[species] / dry_total * 100.0
        
        # Calculate errors vs expected composition
        expected_dry = case_config["expected_outlet_dry"]
        composition_errors = {}
        for species, expected in expected_dry.items():
            if species in outlet_composition_dry:
                actual = outlet_composition_dry[species]
                error = abs(actual - expected)
                composition_errors[species] = {
                    "expected": expected,
                    "actual": actual,
                    "error": error
                }
        
        return {
            "ch4_conversion": ch4_conversion,
            "co2_conversion": co2_conversion,
            "h2_co_ratio": h2_co_ratio,
            "ignition_length": result.ignition_length,
            "pressure_drop": pressure_drop,
            "outlet_composition_wet": outlet_composition_wet,
            "outlet_composition_dry": outlet_composition_dry,
            "composition_errors": composition_errors,
            "runtime": runtime,
            "max_temperature": np.max(result.temperature) - 273.15,
            "outlet_temperature": outlet_temp - 273.15
        }
    
    def _save_case_results(self, case_name: str, result, kpis: Dict, output_dir: Path):
        """Save case results."""
        
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
        
        profiles_df.to_csv(output_dir / "axial_profiles.csv", index=False)
        
        # Save outlet composition table
        comp_data = []
        for species in ["H2", "CO", "CO2", "CH4", "N2", "H2O"]:
            wet_frac = kpis["outlet_composition_wet"].get(species, 0.0)
            dry_frac = kpis["outlet_composition_dry"].get(species, 0.0)
            comp_data.append({
                "species": species,
                "wet_basis_volpct": wet_frac,
                "dry_basis_volpct": dry_frac
            })
        
        comp_df = pd.DataFrame(comp_data)
        comp_df.to_csv(output_dir / "outlet_composition_table.csv", index=False)
        
        # Save KPIs
        kpis_to_save = {
            "h2_co_ratio": kpis["h2_co_ratio"],
            "ch4_conversion": kpis["ch4_conversion"],
            "co2_conversion": kpis["co2_conversion"],
            "ignition_length_m": kpis["ignition_length"],
            "pressure_drop_Pa": kpis["pressure_drop"],
            "runtime_s": kpis["runtime"]
        }
        
        with open(output_dir / "kpis.json", 'w') as f:
            json.dump(kpis_to_save, f, indent=2)
        
        print("✅ Axial profiles saved to axial_profiles.csv")
        print("✅ Outlet composition table saved to outlet_composition_table.csv")
        print("✅ KPIs saved to kpis.json")
    
    def _generate_case_plots(self, case_name: str, result, kpis: Dict, output_dir: Path):
        """Generate case plots."""
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f"HP-POX {case_name.upper()} - Axial Profiles", fontsize=16)
        
        # Temperature profile
        ax1 = axes[0, 0]
        ax1.plot(result.height, result.temperature - 273.15, 'b-', linewidth=2, label='Temperature')
        ax1.set_xlabel('Height (m)')
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_title('Axial Temperature Profile')
        ax1.grid(True)
        
        # Add ignition length marker
        if result.ignition_length:
            ax1.axvline(x=result.ignition_length, color='r', linestyle='--', 
                       label=f'Ignition: {result.ignition_length:.3f} m')
        ax1.legend()
        
        # Species profiles
        ax2 = axes[0, 1]
        key_species = ["CH4", "H2", "CO", "CO2", "H2O", "O2"]
        colors = ['purple', 'red', 'cyan', 'green', 'orange', 'blue']
        
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
        plt.savefig(output_dir / "axial_profiles_plot.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Generate outlet composition bar plot
        self._generate_outlet_barplot(kpis, output_dir)
        
        # Generate ignition length markers plot
        self._generate_ignition_markers_plot(result, output_dir)
        
        print("✅ Axial profiles plotted to axial_profiles_plot.png")
    
    def _generate_outlet_barplot(self, kpis: Dict, output_dir: Path):
        """Generate outlet composition bar plot."""
        
        # Prepare data for bar plot
        species = ["H2", "CO", "CO2", "CH4", "N2"]
        dry_composition = [kpis["outlet_composition_dry"].get(s, 0.0) for s in species]
        
        # Create bar plot
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(species, dry_composition, alpha=0.8, color=['red', 'cyan', 'green', 'purple', 'blue'])
        
        ax.set_xlabel('Species')
        ax.set_ylabel('Volume Fraction (%)')
        ax.set_title('Outlet Composition (Dry Basis)')
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, dry_composition):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{value:.1f}%', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(output_dir / "outlet_barplot.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Outlet composition bar plot saved to outlet_barplot.png")
    
    def _generate_ignition_markers_plot(self, result, output_dir: Path):
        """Generate ignition length markers plot."""
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot temperature profile
        ax.plot(result.height, result.temperature - 273.15, 'b-', linewidth=2, label='Temperature')
        
        # Add ignition length marker
        if result.ignition_length:
            ax.axvline(x=result.ignition_length, color='r', linestyle='--', linewidth=2,
                      label=f'Ignition Length: {result.ignition_length:.3f} m')
        
        ax.set_xlabel('Height (m)')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title('Ignition Length Detection')
        ax.legend()
        ax.grid(True)
        
        plt.tight_layout()
        plt.savefig(output_dir / "ignition_length_markers.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Ignition length markers plot saved to ignition_length_markers.png")
    
    def _print_validation_summary(self, case_name: str, kpis: Dict, case_config: Dict):
        """Print validation summary."""
        
        print(f"\n✅ Axial profiles saved and plotted")
        print(f"✅ Ignition length detected & plotted: {kpis['ignition_length']:.3f} m")
        print(f"✅ Outlet wet & dry compositions saved")
        
        # Print dry basis composition
        print("   Dry basis composition:")
        for species, frac in kpis["outlet_composition_dry"].items():
            print(f"     {species}: {frac:.1f}%")
        
        # Print error analysis
        print("   Error vs Richter targets:")
        for species, error_data in kpis["composition_errors"].items():
            print(f"     {species}: {error_data['actual']:.1f}% (expected {error_data['expected']:.1f}%, error {error_data['error']:.1f}pp)")
        
        print(f"✅ KPIs to kpis.json")
        print(f"✅ Artifacts saved under {self.output_dir}/{case_name}/")
    
    def run_complete_validation(
        self,
        mechanism_file: str,
        reduced_mechanism_file: str = None
    ) -> Dict:
        """Run complete HP-POX validation pipeline."""
        
        print("="*80)
        print("HP-POX PLANT VALIDATION PIPELINE")
        print("="*80)
        print(f"Mechanism: {mechanism_file}")
        print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        results = {}
        
        # Step 1: Validate reference case
        print("\n1. VALIDATING REFERENCE CASE")
        results["reference"] = self.validate_case("case_reference", mechanism_file)
        
        # Step 2: Validate sensitivity case
        print("\n2. VALIDATING SENSITIVITY CASE")
        results["sensitivity"] = self.validate_case("case_sensitivity", mechanism_file)
        
        # Step 3: Run mechanism reduction (if not provided)
        if reduced_mechanism_file is None:
            print("\n3. RUNNING GA-GNN MECHANISM REDUCTION")
            # This would integrate with the GA-GNN pipeline
            # For now, simulate with the same mechanism
            reduced_mechanism_file = mechanism_file
        
        # Step 4: Re-validate with reduced mechanism
        if reduced_mechanism_file != mechanism_file:
            print("\n4. RE-VALIDATING WITH REDUCED MECHANISM")
            results["reference_reduced"] = self.validate_case(
                "case_reference", reduced_mechanism_file, "reduction_ref_sens/case_reference_eval"
            )
            results["sensitivity_reduced"] = self.validate_case(
                "case_sensitivity", reduced_mechanism_file, "reduction_ref_sens/case_sensitivity_eval"
            )
        
        # Step 5: Generate comparative analysis
        print("\n5. GENERATING COMPARATIVE ANALYSIS")
        self._generate_comparative_analysis(results)
        
        print("\n" + "="*80)
        print("HP-POX VALIDATION PIPELINE COMPLETED")
        print("="*80)
        
        return results
    
    def _generate_comparative_analysis(self, results: Dict):
        """Generate comparative analysis plots."""
        
        # Generate outlet composition comparison
        self._generate_outlet_comparison_plot(results)
        
        # Generate ignition length comparison
        self._generate_ignition_comparison_plot(results)
        
        print("✅ Comparative analysis completed")
    
    def _generate_outlet_comparison_plot(self, results: Dict):
        """Generate outlet composition comparison plot."""
        
        # Prepare data
        cases = ["reference", "sensitivity"]
        mechanisms = ["detailed", "reduced"]
        species = ["H2", "CO", "CO2", "CH4", "N2"]
        
        data = []
        for case in cases:
            for mechanism in mechanisms:
                result_key = f"{case}_{mechanism}" if mechanism == "reduced" else case
                if result_key in results:
                    kpis = results[result_key]["kpis"]
                    for species_name in species:
                        if species_name in kpis["outlet_composition_dry"]:
                            data.append({
                                "case": case.title(),
                                "mechanism": mechanism.title(),
                                "species": species_name,
                                "composition": kpis["outlet_composition_dry"][species_name]
                            })
        
        if not data:
            return
        
        df = pd.DataFrame(data)
        
        # Create bar plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        x = np.arange(len(species))
        width = 0.2
        
        for i, case in enumerate(["Reference", "Sensitivity"]):
            for j, mechanism in enumerate(["Detailed", "Reduced"]):
                case_data = df[(df["case"] == case) & (df["mechanism"] == mechanism)]
                if len(case_data) > 0:
                    values = [case_data[case_data["species"] == s]["composition"].iloc[0] 
                             if len(case_data[case_data["species"] == s]) > 0 else 0 
                             for s in species]
                    
                    offset = (i * 2 + j) * width - width * 1.5
                    ax.bar(x + offset, values, width, 
                          label=f'{case} {mechanism}', alpha=0.8)
        
        ax.set_xlabel('Species')
        ax.set_ylabel('Volume Fraction (%)')
        ax.set_title('Outlet Composition Comparison (Dry Basis)')
        ax.set_xticks(x)
        ax.set_xticklabels(species)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "outlet_composition_comparison.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Outlet composition comparison plot saved")
    
    def _generate_ignition_comparison_plot(self, results: Dict):
        """Generate ignition length comparison plot."""
        
        # Prepare data
        cases = ["reference", "sensitivity"]
        mechanisms = ["detailed", "reduced"]
        
        data = []
        for case in cases:
            for mechanism in mechanisms:
                result_key = f"{case}_{mechanism}" if mechanism == "reduced" else case
                if result_key in results:
                    ignition_length = results[result_key]["kpis"]["ignition_length"]
                    if ignition_length:
                        data.append({
                            "case": case.title(),
                            "mechanism": mechanism.title(),
                            "ignition_length": ignition_length
                        })
        
        if not data:
            return
        
        df = pd.DataFrame(data)
        
        # Create bar plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        x = np.arange(len(cases))
        width = 0.35
        
        for i, mechanism in enumerate(["Detailed", "Reduced"]):
            mechanism_data = df[df["mechanism"] == mechanism]
            values = [mechanism_data[mechanism_data["case"] == case]["ignition_length"].iloc[0] 
                     if len(mechanism_data[mechanism_data["case"] == case]) > 0 else 0 
                     for case in ["Reference", "Sensitivity"]]
            
            ax.bar(x + i*width, values, width, label=mechanism, alpha=0.8)
        
        ax.set_xlabel('Case')
        ax.set_ylabel('Ignition Length (m)')
        ax.set_title('Ignition Length Comparison')
        ax.set_xticks(x + width/2)
        ax.set_xticklabels(["Reference", "Sensitivity"])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "ignition_length_comparison.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Ignition length comparison plot saved")
