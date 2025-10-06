"""
Main HP-POX model integrating all physics modules.
Implements the complete 0D-1D reactor model with proper physics.
"""

import numpy as np
from pathlib import Path
from typing import Dict, Any, Tuple
import yaml

from .thermo import ThermoManager
from .mixing import mix_streams
from .psr import PSR
from .pfr import PFR
from .validation import ValidationGates
from .config_cases import CaseConstants


class HPPOXModel:
    """Main HP-POX model with proper physics integration."""
    
    def __init__(self, mechanism_path: str = "gri30.yaml", transport_model: str = "mixture-averaged"):
        """Initialize HP-POX model.
        
        Args:
            mechanism_path: Path to mechanism file
            transport_model: Transport model ('mixture-averaged' for laptop, 'multicomponent' for HPC)
        """
        self.thermo = ThermoManager(mechanism_path, transport_model)
        self.validation = ValidationGates()
        self.case_constants = CaseConstants()
    
    def load_inlet_streams(self, config_path: str) -> Dict[str, Any]:
        """Load inlet stream configuration.
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            Inlet stream configuration
        """
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        return config['inlet']
    
    def create_inlet_state(self, inlet_config: Dict[str, Any], case_id: str) -> Tuple[Any, float]:
        """Create inlet state from configuration.
        
        Args:
            inlet_config: Inlet configuration
            case_id: Case identifier to select correct inlet data
            
        Returns:
            Tuple of (gas_state, total_mass_flow)
        """
        # Get case-specific inlet data
        if 'case1' in case_id.lower():
            case_data = inlet_config['case1']
        elif 'case4' in case_id.lower():
            case_data = inlet_config['case4']
        else:
            raise ValueError(f"Unknown case: {case_id}")
        
        # Convert inlet streams to list format for mixing
        streams = []
        
        # Process each stream
        stream_names = ['primary_steam', 'secondary_steam', 'oxygen', 'nitrogen_optisox', 'natural_gas']
        
        for stream_name in stream_names:
            if stream_name in case_data:
                stream_data = case_data[stream_name]
                
                # Get stream properties
                m_dot = stream_data['mass_kg_per_h'] / 3600.0  # Convert kg/h to kg/s
                T = stream_data['T_C'] + 273.15  # Convert °C to K
                P = self.case_constants.PRESSURE
                
                # Get composition
                if 'comp' in stream_data:
                    # Pure stream
                    Y = stream_data['comp']
                elif 'composition_volpct' in stream_data:
                    # Natural gas with volume percentages
                    Y = stream_data['composition_volpct']
                else:
                    raise ValueError(f"No composition found for stream {stream_name}")
                
                streams.append({
                    'm_dot': m_dot,
                    'T': T,
                    'P': P,
                    'Y': Y
                })
        
        if not streams:
            raise ValueError("No inlet streams found in configuration")
        
        # Mix streams with enthalpy-consistent temperature
        T_mix, Y_mix, m_dot_total = mix_streams(streams, self.thermo)
        
        # Create gas state
        gas_state = self.thermo.create_gas_state(T_mix, self.case_constants.PRESSURE, Y_mix)
        
        print(f"Inlet mixing:")
        print(f"  T_mix = {T_mix-273.15:.1f}°C")
        print(f"  m_dot_total = {m_dot_total:.3f} kg/s")
        print(f"  Y_mix = {Y_mix}")
        
        return gas_state, m_dot_total
    
    def run_case(self, case_id: str, inlet_config: Dict[str, Any], 
                 output_dir: Path, chemistry: bool = True) -> Dict[str, Any]:
        """Run complete HP-POX case.
        
        Args:
            case_id: Case identifier
            inlet_config: Inlet configuration
            output_dir: Output directory
            chemistry: Whether to include chemistry (True) or just wall cooling (False)
            
        Returns:
            Results dictionary
        """
        print(f"\n{'='*60}")
        print(f"RUNNING CASE: {case_id}")
        print(f"{'='*60}")
        
        # Get case configuration
        case_config = self.case_constants.get_case_config(case_id)
        
        # Create inlet state
        inlet_state, m_dot_total = self.create_inlet_state(inlet_config, case_id)
        
        # Create PSR with correct volume from Richter benchmark
        psr_volume = 16.136e-6  # m³ (16136 cm³ from Richter benchmark)
        psr = PSR(psr_volume, self.thermo)
        
        # Solve PSR with fixed burner heat loss using continuation
        q_burner_kw = case_config['burner_heat_loss']
        psr_outlet, psr_converged, psr_diagnostics = psr.solve_psr_continuation(
            inlet_state, m_dot_total, case_config['pressure'], q_burner_kw
        )
        
        # HARD PRECONDITION: Check PSR convergence
        if not psr_converged:
            print(f"ERROR: PSR did not converge—PFR not started")
            print(f"  PSR diagnostics: {psr_diagnostics}")
            raise RuntimeError(f"PSR convergence failed for case {case_id}")
        
        # Validate PSR
        T_target = case_config['target_temperature']
        psr_validation = self.validation.validate_psr(psr_outlet, T_target)
        print(f"PSR validation: {'PASS' if psr_validation['psr_temperature_gate'] else 'FAIL'}")
        
        # Extract Q_burner for compatibility
        Q_burner = psr_diagnostics['Q_burner']
        
        # Create PFR
        geometry = case_config['geometry']
        pfr = PFR(geometry['length_m'], geometry['diameter_m'], self.thermo)
        
        # Step 0: Test PFR chemistry-off
        print(f"\n{'='*40}")
        print("STEP 0: PFR CHEMISTRY-OFF TEST")
        print(f"{'='*40}")
        
        pfr_results_off = pfr.solve_pfr_operator_split(
            psr_outlet, m_dot_total, case_config['pressure'],
            case_config['wall_heat_loss'], chemistry_on=False
        )
        
        # Validate chemistry-off test
        cp_avg = np.mean([self.thermo.create_gas_state(T, case_config['pressure'], Y).cp_mass 
                          for T, Y in zip(pfr_results_off['temperature_K'], pfr_results_off['mass_fractions'])])
        
        chemistry_off_validation = self.validation.validate_pfr_chemistry_off(
            pfr_results_off['z_m'], pfr_results_off['temperature_K'],
            m_dot_total, cp_avg, case_config['wall_heat_loss']
        )
        
        print(f"Chemistry-off validation: {'PASS' if chemistry_off_validation['chemistry_off_gate'] else 'FAIL'}")
        
        # Step 1: PFR chemistry-on/off based on parameter
        print(f"\n{'='*40}")
        print(f"STEP 1: PFR CHEMISTRY-{'ON' if chemistry else 'OFF'}")
        print(f"{'='*40}")
        
        pfr_results_on = pfr.solve_pfr_operator_split(
            psr_outlet, m_dot_total, case_config['pressure'],
            case_config['wall_heat_loss'], chemistry_on=chemistry, Y_inlet=psr_outlet.mass_fractions
        )
        
        # Validate results
        outlet_state = self.thermo.create_gas_state(
            pfr_results_on['temperature_K'][-1],
            case_config['pressure'],
            pfr_results_on['mass_fractions'][-1, :]
        )
        
        # Dry basis composition
        dry_comp = outlet_state.get_dry_basis_composition()
        h2_co_ratio = dry_comp.get('H2', 0.0) / dry_comp.get('CO', 1e-10)
        
        print(f"\nOutlet composition (dry basis):")
        for species, fraction in dry_comp.items():
            print(f"  {species}: {fraction:.4f}")
        print(f"  H2/CO ratio: {h2_co_ratio:.2f}")
        
        # KPI validation
        kpi_validation = self.validation.validate_kpis(outlet_state, case_config['target_kpis'])
        
        # Compile results
        results = {
            'case_id': case_id,
            'psr_outlet': psr_outlet,
            'psr_validation': psr_validation,
            'psr_diagnostics': psr_diagnostics,
            'Q_burner': Q_burner,
            'pfr_results_off': pfr_results_off,
            'pfr_results_on': pfr_results_on,
            'chemistry_off_validation': chemistry_off_validation,
            'outlet_state': outlet_state,
            'dry_composition': dry_comp,
            'h2_co_ratio': h2_co_ratio,
            'kpi_validation': kpi_validation
        }
        
        # Print validation summary
        self.validation.print_validation_summary({
            'PSR Temperature': psr_validation['psr_temperature_gate'],
            'Chemistry-Off Test': chemistry_off_validation['chemistry_off_gate'],
            'H2/CO Ratio': h2_co_ratio,
            'Target H2/CO': case_config['target_kpis']['h2_co_ratio']
        })
        
        return results
    
    def save_results(self, results: Dict[str, Any], output_dir: Path):
        """Save results to output directory.
        
        Args:
            results: Results dictionary
            output_dir: Output directory
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save axial profiles
        pfr_results = results['pfr_results_on']
        
        # Create profiles CSV
        profiles_data = {
            'z_m': pfr_results['z_m'],
            'T_C': pfr_results['temperature_K'] - 273.15,
            'rho_kg_m3': pfr_results['density_kg_m3'],
            'u_m_s': pfr_results['velocity_m_s'],
            'residence_time_s': pfr_results['residence_time_s']
        }
        
        # Add species profiles
        for i, species in enumerate(pfr_results['species_names']):
            profiles_data[f'Y_{species}'] = pfr_results['mass_fractions'][:, i]
            profiles_data[f'X_{species}'] = pfr_results['mole_fractions'][:, i]
        
        # Save to CSV
        import pandas as pd
        df = pd.DataFrame(profiles_data)
        df.to_csv(output_dir / 'axial_profiles.csv', index=False)
        
        # Save KPIs (augment with outlet dry-basis computed from axial_profiles.csv)
        kpis = {
            'case_id': results['case_id'],
            'h2_co_ratio': results['h2_co_ratio'],
            'psr_temperature_C': results['psr_outlet'].temperature - 273.15,
            'pfr_outlet_temperature_C': results['outlet_state'].temperature - 273.15,
            'ignition_length_peak_m': pfr_results['ignition_length_peak_m'],
            'ignition_length_threshold_m': pfr_results['ignition_length_threshold_m'],
            'residence_time_s': pfr_results['residence_time_s'][-1],
            'Q_burner_kW': results['Q_burner'] / 1000.0
        }
        
        import json
        # Compute outlet dry-basis from X_* columns in the just-written axial_profiles.csv
        try:
            df_outlet = df.tail(1).reset_index(drop=True)
            # Collect X_* columns actually present
            x_cols = [c for c in df_outlet.columns if c.startswith('X_')]
            # Map species name without prefix
            species_mapping = {c[2:]: c for c in x_cols}
            # Build dry set excluding H2O
            dry_pairs = [(s, df_outlet[c].iloc[0]) for s, c in species_mapping.items() if s != 'H2O']
            total = float(sum(v for _, v in dry_pairs)) if dry_pairs else 0.0
            sim_dry = {}
            if total > 0.0 and all(np.isfinite(v) for _, v in dry_pairs):
                for s, v in dry_pairs:
                    sim_dry[s] = float(v) / total
                dry_basis_computed = True
            else:
                dry_basis_computed = False
            # Augment KPIs JSON
            kpis_aug = dict(kpis)
            kpis_aug.update({
                'species_columns': x_cols,
                'species_mapping': species_mapping,
                'dry_basis_computed': dry_basis_computed,
                'sim_dry': sim_dry if dry_basis_computed else {},
                'targets_available': False,
                'files_used': {
                    'axial_profiles': str((output_dir / 'axial_profiles.csv').as_posix()),
                    'kpis': str((output_dir / 'kpis.json').as_posix()),
                    'outlet_composition': None
                }
            })
            with open(output_dir / 'kpis.json', 'w') as f:
                json.dump(kpis_aug, f, indent=2)
        except Exception:
            # Fallback to original KPIs without dry-basis
            with open(output_dir / 'kpis.json', 'w') as f:
                json.dump(kpis, f, indent=2)
        
        # Save PSR diagnostics
        if 'psr_diagnostics' in results:
            with open(output_dir / 'diagnostics.json', 'w') as f:
                json.dump(results['psr_diagnostics'], f, indent=2)
        
        print(f"Results saved to: {output_dir}")


def run_case(case_id: str, chemistry: bool = True, outdir: str = "hp_pox_results", mechanism: str = "gri30.yaml", transport: str = "mixture-averaged", procs: int = 1) -> dict:
    """Run a single HP-POX case.
    
    Args:
        case_id: Case identifier (e.g., "A_case1_richter")
        chemistry: Whether to include chemistry (True) or just wall cooling (False)
        outdir: Output directory
        mechanism: Path to chemical mechanism file
        transport: Transport model ('mixture-averaged' for laptop, 'multicomponent' for HPC)
        procs: Number of processes for parallel execution (default: 1 for laptop)
        
    Returns:
        Dictionary with KPIs and results
    """
    model = HPPOXModel(mechanism, transport)
    inlet_config = model.load_inlet_streams("config.yaml")
    results = model.run_case(case_id, inlet_config, Path(outdir), chemistry)
    model.save_results(results, Path(outdir) / case_id)
    return results


def run_all(cases: list[str] | None = None, chemistry: bool = True, outdir: str = "hp_pox_results", mechanism: str = "gri30.yaml") -> dict:
    """Run all HP-POX cases.
    
    Args:
        cases: List of case IDs to run (None = run all available)
        chemistry: Whether to include chemistry (True) or just wall cooling (False)
        outdir: Output directory
        mechanism: Path to chemical mechanism file
        
    Returns:
        Dictionary with summary of all cases
    """
    from .config_cases import get_available_cases
    
    if cases is None:
        cases = get_available_cases()
    
    summary = {}
    for case_id in cases:
        print(f"\n{'='*60}")
        print(f"RUNNING CASE: {case_id}")
        print(f"{'='*60}")
        
        try:
            results = run_case(case_id, chemistry, outdir, mechanism)
            summary[case_id] = {
                'status': 'success',
                'psr_temperature_C': results['psr_outlet'].temperature - 273.15,
                'pfr_outlet_temperature_C': results['outlet_state'].temperature - 273.15,
                'h2_co_ratio': results['h2_co_ratio'],
                'ignition_length_m': results['pfr_results_on']['ignition_length_peak_m']
            }
        except Exception as e:
            summary[case_id] = {
                'status': 'failed',
                'error': str(e)
            }
            print(f"ERROR in {case_id}: {e}")
    
    return summary


def main():
    """Main function for testing."""
    # Test the model
    model = HPPOXModel("gri30.yaml")
    
    # Load configuration
    inlet_config = model.load_inlet_streams("config.yaml")
    
    # Run case
    case_id = "A_case1_richter"
    results = model.run_case(case_id, inlet_config, Path("test_output"))
    
    # Save results
    model.save_results(results, Path("test_output"))


if __name__ == "__main__":
    main()
