"""
Data I/O module for HP-POX model.

Handles loading/saving configuration files, printing feed headers,
and managing output data structures.
"""

import yaml
import json
import csv
from pathlib import Path
from typing import Dict, Any, Optional
import numpy as np


class DataIO:
    """Data input/output operations."""
    
    def __init__(self):
        """Initialize data I/O handler."""
        self.config_dir = Path(__file__).parent.parent / 'config'
        self.mechanisms_dir = Path(__file__).parent.parent / 'mechanisms'
    
    def load_case_config(self, case: str) -> Dict[str, Any]:
        """Load case-specific configuration."""
        config_file = self.config_dir / f'{case}.yaml'
        
        if not config_file.exists():
            # Create default case configurations
            self._create_default_case_configs()
        
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def load_geometry_config(self, geometry_path: Optional[str] = None) -> Dict[str, Any]:
        """Load geometry configuration."""
        if geometry_path:
            config_file = Path(geometry_path)
        else:
            config_file = self.config_dir / 'geometry.yaml'
        
        if not config_file.exists():
            # Create default geometry configuration
            self._create_default_geometry_config()
        
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def load_heatloss_config(self, heatloss_path: Optional[str] = None) -> Dict[str, Any]:
        """Load heat loss configuration."""
        if heatloss_path:
            config_file = Path(heatloss_path)
        else:
            # Use default heat loss model
            return {'model': 'constant', 'wall_temperature_K': 418.0}
        
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def _create_default_case_configs(self) -> None:
        """Create default case configuration files."""
        # Case 1 configuration
        case1_config = {
            'name': 'Case 1 (50 bar(g), target 1200°C)',
            'pressure_gauge_bar': 50.0,
            'target_temperature_C': 1200.0,
            'feeds': {
                'primary_steam': {
                    'mass_flow_kg_h': 78.02,
                    'temperature_C': 289.6,
                    'composition': {'H2O': 1.0}
                },
                'secondary_steam': {
                    'mass_flow_kg_h': 23.15,
                    'temperature_C': 240.2,
                    'composition': {'H2O': 1.0}
                },
                'oxygen': {
                    'mass_flow_kg_h': 265.37,
                    'temperature_C': 240.2,
                    'composition': {'O2': 1.0}
                },
                'natural_gas': {
                    'mass_flow_kg_h': 224.07,
                    'temperature_C': 359.1,
                    'composition': {
                        'CH4': 95.97,
                        'CO2': 0.34,
                        'C2H6': 2.26,
                        'C3H8': 0.46,
                        'nC4H10': 0.05,
                        'iC4H10': 0.06,
                        'N2': 0.86
                    }
                }
            },
            'heat_loss': {
                'burner_cooling_kW': 34.29,
                'reactor_wall_kW': 18.50
            },
            'targets': {
                'dry_composition_pct': {
                    'H2': 48.27,
                    'CO': 23.79,
                    'CO2': 4.19,
                    'CH4': 3.76,
                    'N2': 0.65
                },
                'outlet_temperature_C': 1201.0
            }
        }
        
        # Case 4 configuration
        case4_config = {
            'name': 'Case 4 (50 bar(g), target 1400°C)',
            'pressure_gauge_bar': 50.0,
            'target_temperature_C': 1400.0,
            'feeds': {
                'primary_steam': {
                    'mass_flow_kg_h': 64.39,
                    'temperature_C': 277.5,
                    'composition': {'H2O': 1.0}
                },
                'secondary_steam': {
                    'mass_flow_kg_h': 22.77,
                    'temperature_C': 239.9,
                    'composition': {'H2O': 1.0}
                },
                'oxygen': {
                    'mass_flow_kg_h': 280.25,
                    'temperature_C': 239.9,
                    'composition': {'O2': 1.0}
                },
                'natural_gas': {
                    'mass_flow_kg_h': 195.90,
                    'temperature_C': 355.2,
                    'composition': {
                        'CH4': 96.20,
                        'CO2': 0.26,
                        'C2H6': 2.09,
                        'C3H8': 0.47,
                        'nC4H10': 0.05,
                        'iC4H10': 0.06,
                        'N2': 0.87
                    }
                }
            },
            'heat_loss': {
                'burner_cooling_kW': 36.36,
                'reactor_wall_kW': 22.92
            },
            'targets': {
                'dry_composition_pct': {
                    'H2': 48.06,
                    'CO': 25.61,
                    'CO2': 3.89,
                    'CH4': 0.06,
                    'N2': 0.67
                },
                'outlet_temperature_C': 1401.0
            }
        }
        
        # Write configuration files
        with open(self.config_dir / 'case1.yaml', 'w') as f:
            yaml.dump(case1_config, f, default_flow_style=False, indent=2)
        
        with open(self.config_dir / 'case4.yaml', 'w') as f:
            yaml.dump(case4_config, f, default_flow_style=False, indent=2)
    
    def _create_default_geometry_config(self) -> None:
        """Create default geometry configuration."""
        geometry_config = {
            'psr': {
                'volume_m3': 0.016136
            },
            'pfr': {
                'length_m': 2.32,
                'diameter_m': 0.49
            },
            'wall': {
                'temperature_K': 418.0
            }
        }
        
        with open(self.config_dir / 'geometry.yaml', 'w') as f:
            yaml.dump(geometry_config, f, default_flow_style=False, indent=2)
    
    def print_feed_header(self, case_config: Dict[str, Any]) -> None:
        """Print detailed feed header."""
        print("\n" + "="*80)
        print(f"HP-POX {case_config['name'].upper()}")
        print("="*80)
        print("Oxygen-blown configuration (no air, no N2 ballast)")
        
        pressure_abs = (case_config['pressure_gauge_bar'] + 1.0) * 1e5
        print(f"Pressure: {case_config['pressure_gauge_bar']} bar(g) = {pressure_abs/1e5:.1f} bar(abs)")
        print(f"Target temperature: {case_config['target_temperature_C']}°C")
        
        # Natural gas composition
        ng_comp = case_config['feeds']['natural_gas']['composition']
        print(f"\nNatural Gas Composition (vol%):")
        for species, pct in ng_comp.items():
            print(f"  {species}: {pct}%")
        
        # Calculate molar flows
        feeds = case_config['feeds']
        ng_mass_kg_s = feeds['natural_gas']['mass_flow_kg_h'] / 3600.0
        o2_mass_kg_s = feeds['oxygen']['mass_flow_kg_h'] / 3600.0
        h2o_mass_kg_s = (feeds['primary_steam']['mass_flow_kg_h'] + 
                        feeds['secondary_steam']['mass_flow_kg_h']) / 3600.0
        
        # Approximate molecular weights for molar flow calculation
        ng_mw = 16.0  # Approximate for CH4-dominated NG
        ng_mol_s = ng_mass_kg_s / ng_mw
        o2_mol_s = o2_mass_kg_s / 32.0
        h2o_mol_s = h2o_mass_kg_s / 18.0
        
        print(f"\nFeed Streams (kmol/s):")
        print(f"  NG: {ng_mol_s:.4f} kmol/s ({feeds['natural_gas']['mass_flow_kg_h']:.1f} kg/h)")
        print(f"  O2: {o2_mol_s:.4f} kmol/s ({feeds['oxygen']['mass_flow_kg_h']:.1f} kg/h)")
        print(f"  H2O: {h2o_mol_s:.4f} kmol/s ({(feeds['primary_steam']['mass_flow_kg_h'] + feeds['secondary_steam']['mass_flow_kg_h']):.1f} kg/h)")
        
        # Calculate ratios
        ch4_mol_s = ng_mol_s * ng_comp['CH4'] / 100.0
        o2_c_ratio = o2_mol_s / ch4_mol_s
        steam_c_ratio = h2o_mol_s / ch4_mol_s
        
        print(f"\nKey Ratios:")
        print(f"  O2:CH4 = {o2_c_ratio:.3f}")
        print(f"  Steam:C = {steam_c_ratio:.3f}")
        
        # Stream temperatures
        print(f"\nStream Temperatures:")
        for name, feed in feeds.items():
            print(f"  {name}: {feed['temperature_C']:.1f}°C")
        
        print(f"\nHeat Loss (from Richter benchmark):")
        print(f"  Burner cooling: {case_config['heat_loss']['burner_cooling_kW']:.1f} kW")
        print(f"  Reactor wall: {case_config['heat_loss']['reactor_wall_kW']:.1f} kW")
        print("="*80)
    
    def save_axial_profiles(self, pfr_result: Dict[str, Any], output_dir: Path) -> None:
        """Save axial profiles to CSV."""
        csv_file = output_dir / 'axial_profiles.csv'
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            header = ['z_m', 'T_K', 'P_Pa', 'u_m_s', 'residence_time_s']
            header.extend([f'Y_{species}' for species in pfr_result['species_names']])
            header.extend([f'X_{species}' for species in pfr_result['species_names']])
            writer.writerow(header)
            
            # Data
            n_cells = len(pfr_result['z_m'])
            for i in range(n_cells):
                row = [
                    pfr_result['z_m'][i],
                    pfr_result['temperature_K'][i],
                    pfr_result['pressure_Pa'][i],
                    pfr_result['velocity_m_s'][i],
                    pfr_result['residence_time_s'][i]
                ]
                row.extend(pfr_result['mass_fractions'][:, i])
                row.extend(pfr_result['mole_fractions'][:, i])
                writer.writerow(row)
    
    def save_outlet_composition(self, kpi_result: Dict[str, Any], output_dir: Path) -> None:
        """Save outlet composition table to CSV."""
        csv_file = output_dir / 'outlet_composition_table.csv'
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Wet composition
            writer.writerow(['Wet Composition (vol%)'])
            writer.writerow(['Species', 'Mole Fraction (%)'])
            for species, frac in kpi_result['wet_composition_pct'].items():
                writer.writerow([species, f"{frac:.2f}"])
            
            writer.writerow([])
            
            # Dry composition
            writer.writerow(['Dry Composition (vol%)'])
            writer.writerow(['Species', 'Mole Fraction (%)', 'Target (%)', 'Error (pp)'])
            for species, frac in kpi_result['dry_composition_pct'].items():
                target = kpi_result['targets']['dry_composition_pct'].get(species, 0.0)
                error = abs(frac - target)
                writer.writerow([species, f"{frac:.2f}", f"{target:.2f}", f"{error:.2f}"])
    
    def save_kpis(self, kpi_result: Dict[str, Any], output_dir: Path) -> None:
        """Save KPIs to JSON."""
        json_file = output_dir / 'kpis.json'
        
        with open(json_file, 'w') as f:
            json.dump(kpi_result, f, indent=2)
    
    def save_readme(self, case_config: Dict[str, Any], geometry_config: Dict[str, Any], 
                   output_dir: Path) -> None:
        """Save README with configuration details."""
        readme_file = output_dir / 'readme.txt'
        
        with open(readme_file, 'w') as f:
            f.write("HP-POX Model Results\n")
            f.write("===================\n\n")
            f.write(f"Case: {case_config['name']}\n")
            f.write(f"Pressure: {case_config['pressure_gauge_bar']} bar(g)\n")
            f.write(f"Target Temperature: {case_config['target_temperature_C']}°C\n\n")
            
            f.write("Geometry:\n")
            f.write(f"  PSR Volume: {geometry_config['psr']['volume_m3']:.6f} m³\n")
            f.write(f"  PFR Length: {geometry_config['pfr']['length_m']:.2f} m\n")
            f.write(f"  PFR Diameter: {geometry_config['pfr']['diameter_m']:.2f} m\n")
            f.write(f"  Wall Temperature: {geometry_config['wall']['temperature_K']:.1f} K\n\n")
            
            f.write("Heat Loss:\n")
            f.write(f"  Burner Cooling: {case_config['heat_loss']['burner_cooling_kW']:.1f} kW\n")
            f.write(f"  Reactor Wall: {case_config['heat_loss']['reactor_wall_kW']:.1f} kW\n\n")
            
            f.write("Files:\n")
            f.write("  axial_profiles.csv - Axial profiles (z, T, P, u, Y_i, X_i)\n")
            f.write("  axial_profiles_plot.png - 4-panel plot of axial profiles\n")
            f.write("  ignition_length_markers.png - Temperature with ignition markers\n")
            f.write("  outlet_composition_table.csv - Wet and dry compositions\n")
            f.write("  outlet_barplot.png - Outlet composition bar chart\n")
            f.write("  kpis.json - Complete KPI dictionary\n")
