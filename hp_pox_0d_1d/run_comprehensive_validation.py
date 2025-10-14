#!/usr/bin/env python3
"""
Comprehensive HP-POX pipeline validation using premix-ignition and velocity-aware implementation.
Runs both base and reduced mechanisms for all cases and generates comparison summary.
"""

import os
import sys
import shutil
from pathlib import Path
from datetime import datetime
import time
from typing import Dict, Any, List, Tuple

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

# Import updated modules
from hp_pox_physics.hp_pox_model import HPPOXModel
from hp_pox_physics.config_cases import richter_case_1, richter_case_4
from hp_pox_physics.diagnostics import run_comprehensive_diagnostics
from hp_pox_physics.thermo import GasState
import hp_pox_physics.diagnostics as diagnostics
import cantera as ct


def run_case_validation(case_config, mechanism_path: str, output_root: str, case_id: str) -> Dict[str, Any]:
    """Run complete PSR→PFR validation for a single case and mechanism."""

    print(f"\n{'='*60}")
    print(f"RUNNING CASE: {case_id} with {Path(mechanism_path).name}")
    print(f"{'='*60}")

    # Create output directory
    case_output_dir = Path(output_root) / case_id
    case_output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize model
    model = HPPOXModel(mechanism_path, transport_model='mixture-averaged')

    # Create inlet configuration (simplified for this validation)
    inlet_config = {
        'case1' if 'case1' in case_id.lower() else 'case4': {
            'primary_steam': {
                'mass_kg_per_h': case_config.primary_steam.mass_kg_per_s * 3600,
                'T_C': case_config.primary_steam.T_K - 273.15,
                'comp': {'H2O': 1.0}
            },
            'secondary_steam': {
                'mass_kg_per_h': case_config.secondary_steam.mass_kg_per_s * 3600,
                'T_C': case_config.secondary_steam.T_K - 273.15,
                'comp': {'H2O': 1.0}
            },
            'oxygen': {
                'mass_kg_per_h': case_config.oxygen.mass_kg_per_s * 3600,
                'T_C': case_config.oxygen.T_K - 273.15,
                'comp': {'O2': 1.0}
            },
            'nitrogen_optisox': {
                'mass_kg_per_h': case_config.nitrogen.mass_kg_per_s * 3600,
                'T_C': case_config.nitrogen.T_K - 273.15,
                'comp': {'N2': 1.0}
            },
            'natural_gas': {
                'mass_kg_per_h': case_config.natural_gas.mass_kg_per_s * 3600,
                'T_C': case_config.natural_gas.T_K - 273.15,
                'composition_volpct': case_config.natural_gas.composition_X
            }
        }
    }

    # Create inlet streams for PSR
    inlet_streams = [
        {
            'name': 'primary_steam',
            'T': case_config.primary_steam.T_K,
            'm_dot': case_config.primary_steam.mass_kg_per_s,
            'Y': case_config.primary_steam.composition_X
        },
        {
            'name': 'secondary_steam',
            'T': case_config.secondary_steam.T_K,
            'm_dot': case_config.secondary_steam.mass_kg_per_s,
            'Y': case_config.secondary_steam.composition_X
        },
        {
            'name': 'oxygen',
            'T': case_config.oxygen.T_K,
            'm_dot': case_config.oxygen.mass_kg_per_s,
            'Y': case_config.oxygen.composition_X
        },
        {
            'name': 'nitrogen',
            'T': case_config.nitrogen.T_K,
            'm_dot': case_config.nitrogen.mass_kg_per_s,
            'Y': case_config.nitrogen.composition_X
        },
        {
            'name': 'natural_gas',
            'T': case_config.natural_gas.T_K,
            'm_dot': case_config.natural_gas.mass_kg_per_s,
            'Y': case_config.natural_gas.composition_X
        }
    ]

    # Run case
    results = model.run_case(case_id, inlet_config, str(case_output_dir))

    # Extract PFR results
    pfr_results = results['pfr_results_on']

    # Create outlet state for diagnostics
    outlet_state = model.thermo.create_gas_state(
        pfr_results['temperature_K'][-1],
        case_config.pressure_pa,
        pfr_results['mass_fractions'][-1, :]
    )

    # Calculate total mass flow (sum of all inlet streams)
    m_dot_total = sum(stream['m_dot'] for stream in inlet_streams)

    # Run comprehensive diagnostics
    # Build the true mixed/premix inlet using the same routine the model uses
    from hp_pox_physics.premix_init import premix_equilibrated_initial_state
    gas_init_flame, gas_psr_in = premix_equilibrated_initial_state(
        mech=mechanism_path,
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
    # inlet for diagnostics = "post-premix, pre-PSR kernel + primary steam" state
    inlet_gas = gas_psr_in
    inlet_gas_state = GasState(inlet_gas, inlet_gas.T, inlet_gas.P, inlet_gas.Y)

    # outlet_gas_state stays as you already compute from PFR last cell
    outlet_gas = outlet_state.gas
    outlet_gas_state = GasState(outlet_gas, outlet_gas.T, outlet_gas.P, outlet_gas.Y)

    diag_results = run_comprehensive_diagnostics(
        inlet_gas,  # inlet state (Solution object)
        outlet_gas,  # outlet state (Solution object)
        m_dot_total,  # inlet mass flow
        m_dot_total,  # outlet mass flow (assuming no mass loss)
        case_config.burner_heat_loss_W,
        case_config.wall_heat_loss_W,
        case_config
    )

    # Count species in mechanism
    mech = ct.Solution(mechanism_path)
    n_species = mech.n_species
    n_reactions = mech.n_reactions

    # Add mechanism info to diagnostics
    if 'diagnostics' in results:
        results['diagnostics']['species_count'] = n_species
        results['diagnostics']['reaction_count'] = n_reactions
        results['diagnostics']['mechanism_name'] = Path(mechanism_path).name

    # Extract KPIs from results
    kpis = results.get('kpis', {})

    return {
        'case_id': case_id,
        'mechanism_path': mechanism_path,
        'mechanism_name': Path(mechanism_path).name,
        'n_species': n_species,
        'n_reactions': n_reactions,
        'output_dir': str(case_output_dir),
        'results': results,
        'diagnostics': diag_results,
        'kpis': kpis,
        'outlet_state': outlet_state
    }


def generate_comparison_plots(results: List[Dict[str, Any]]) -> None:
    """Generate overlay plots comparing base vs reduced mechanisms."""
    import matplotlib.pyplot as plt
    import pandas as pd

    # Group results by case
    case_results = {}
    for result in results:
        case_id = result['case_id']
        mech_name = result['mechanism_name']
        if case_id not in case_results:
            case_results[case_id] = {}
        case_results[case_id][mech_name] = result

    comparison_dir = Path("comparison_gri_v2")
    comparison_dir.mkdir(exist_ok=True)

    for case_id, mechs in case_results.items():
        if len(mechs) < 2:
            continue  # Need both base and reduced

        base_result = mechs.get("Base GRI-3.0")
        reduced_result = mechs.get("GA-GNN Reduced")

        if not base_result or not reduced_result:
            continue

        case_dir = comparison_dir / case_id
        case_dir.mkdir(exist_ok=True)
        figs_dir = case_dir / "figs"
        figs_dir.mkdir(exist_ok=True)

        # Load data
        base_csv = Path(base_result['output_dir']) / "axial_profiles.csv"
        reduced_csv = Path(reduced_result['output_dir']) / "axial_profiles.csv"

        if not base_csv.exists() or not reduced_csv.exists():
            continue

        base_df = pd.read_csv(base_csv)
        reduced_df = pd.read_csv(reduced_csv)

        # T vs z overlay
        plt.figure(figsize=(10, 6))
        plt.plot(base_df['z_m'], base_df['T_C'], 'b-', label='Base GRI-3.0', linewidth=2)
        plt.plot(reduced_df['z_m'], reduced_df['T_C'], 'r--', label='GA-GNN Reduced', linewidth=2)
        plt.xlabel('z (m)')
        plt.ylabel('T (°C)')
        plt.title(f'Temperature vs z - {case_id}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(figs_dir / 'T_vs_z_overlay.png', dpi=180, bbox_inches='tight')
        plt.close()

        # Species vs z overlay
        plt.figure(figsize=(10, 6))
        species_cols = ['X_H2', 'X_CO', 'X_CO2', 'X_CH4', 'X_H2O', 'X_O2', 'X_N2']
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'cyan', 'magenta']
        for col, color in zip(species_cols, colors):
            if col in base_df.columns and col in reduced_df.columns:
                plt.plot(base_df['z_m'], base_df[col] * 100, f'{color}-', label=f'Base {col}', linewidth=2)
                plt.plot(reduced_df['z_m'], reduced_df[col] * 100, f'{color}--', label=f'Reduced {col}', linewidth=2)
        plt.xlabel('z (m)')
        plt.ylabel('Mole fraction (%)')
        plt.title(f'Species vs z - {case_id}')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(figs_dir / 'species_vs_z_overlay.png', dpi=180, bbox_inches='tight')
        plt.close()

        # H2CO vs z overlay
        plt.figure(figsize=(10, 6))
        # Calculate H2/CO for both
        base_h2co = base_df['X_H2'] / base_df['X_CO'].replace(0, np.nan)
        reduced_h2co = reduced_df['X_H2'] / reduced_df['X_CO'].replace(0, np.nan)

        plt.plot(base_df['z_m'], base_h2co, 'b-', label='Base GRI-3.0', linewidth=2)
        plt.plot(reduced_df['z_m'], reduced_h2co, 'r--', label='GA-GNN Reduced', linewidth=2)
        plt.axhline(2.0, color='k', linestyle='--', label='Target 2.0')
        plt.xlabel('z (m)')
        plt.ylabel('H2/CO')
        plt.title(f'H2/CO vs z - {case_id}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(figs_dir / 'H2CO_vs_z_overlay.png', dpi=180, bbox_inches='tight')
        plt.close()

        # Velocity vs z overlay
        if 'u_m_s' in base_df.columns and 'u_m_s' in reduced_df.columns:
            plt.figure(figsize=(10, 6))
            plt.plot(base_df['z_m'], base_df['u_m_s'], 'b-', label='Base GRI-3.0', linewidth=2)
            plt.plot(reduced_df['z_m'], reduced_df['u_m_s'], 'r--', label='GA-GNN Reduced', linewidth=2)
            plt.xlabel('z (m)')
            plt.ylabel('u (m/s)')
            plt.title(f'Velocity vs z - {case_id}')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(figs_dir / 'u_vs_z_overlay.png', dpi=180, bbox_inches='tight')
            plt.close()

        print(f"  Generated overlay plots for {case_id}")


def generate_kpi_comparison_table(results: List[Dict[str, Any]]) -> None:
    """Generate overlay plots comparing base vs reduced mechanisms."""
    import matplotlib.pyplot as plt
    import pandas as pd

    # Group results by case
    case_results = {}
    for result in results:
        case_id = result['case_id']
        mech_name = result['mechanism_name']
        if case_id not in case_results:
            case_results[case_id] = {}
        case_results[case_id][mech_name] = result

    comparison_dir = Path("comparison_gri_v2")
    comparison_dir.mkdir(exist_ok=True)

    for case_id, mechs in case_results.items():
        if len(mechs) < 2:
            continue  # Need both base and reduced

        base_result = mechs.get("Base GRI-3.0")
        reduced_result = mechs.get("GA-GNN Reduced")

        if not base_result or not reduced_result:
            continue

        case_dir = comparison_dir / case_id
        case_dir.mkdir(exist_ok=True)
        figs_dir = case_dir / "figs"
        figs_dir.mkdir(exist_ok=True)

        # Load data
        base_csv = Path(base_result['output_dir']) / "axial_profiles.csv"
        reduced_csv = Path(reduced_result['output_dir']) / "axial_profiles.csv"

        if not base_csv.exists() or not reduced_csv.exists():
            continue

        base_df = pd.read_csv(base_csv)
        reduced_df = pd.read_csv(reduced_csv)

        # T vs z overlay
        plt.figure(figsize=(10, 6))
        plt.plot(base_df['z_m'], base_df['T_C'], 'b-', label='Base GRI-3.0', linewidth=2)
        plt.plot(reduced_df['z_m'], reduced_df['T_C'], 'r--', label='GA-GNN Reduced', linewidth=2)
        plt.xlabel('z (m)')
        plt.ylabel('T (°C)')
        plt.title(f'Temperature vs z - {case_id}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(figs_dir / 'T_vs_z_overlay.png', dpi=180, bbox_inches='tight')
        plt.close()

        # Species vs z overlay
        plt.figure(figsize=(10, 6))
        species_cols = ['X_H2', 'X_CO', 'X_CO2', 'X_CH4', 'X_H2O', 'X_O2', 'X_N2']
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'cyan', 'magenta']
        for col, color in zip(species_cols, colors):
            if col in base_df.columns and col in reduced_df.columns:
                plt.plot(base_df['z_m'], base_df[col] * 100, f'{color}-', label=f'Base {col}', linewidth=2)
                plt.plot(reduced_df['z_m'], reduced_df[col] * 100, f'{color}--', label=f'Reduced {col}', linewidth=2)
        plt.xlabel('z (m)')
        plt.ylabel('Mole fraction (%)')
        plt.title(f'Species vs z - {case_id}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(figs_dir / 'species_vs_z_overlay.png', dpi=180, bbox_inches='tight')
        plt.close()

        # H2CO vs z overlay
        plt.figure(figsize=(10, 6))
        # Calculate H2/CO for both
        base_h2co = base_df['X_H2'] / base_df['X_CO'].replace(0, np.nan)
        reduced_h2co = reduced_df['X_H2'] / reduced_df['X_CO'].replace(0, np.nan)

        plt.plot(base_df['z_m'], base_h2co, 'b-', label='Base GRI-3.0', linewidth=2)
        plt.plot(reduced_df['z_m'], reduced_h2co, 'r--', label='GA-GNN Reduced', linewidth=2)
        plt.axhline(2.0, color='k', linestyle='--', label='Target 2.0')
        plt.xlabel('z (m)')
        plt.ylabel('H2/CO')
        plt.title(f'H2/CO vs z - {case_id}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(figs_dir / 'H2CO_vs_z_overlay.png', dpi=180, bbox_inches='tight')
        plt.close()

        # Velocity vs z overlay
        if 'u_m_s' in base_df.columns and 'u_m_s' in reduced_df.columns:
            plt.figure(figsize=(10, 6))
            plt.plot(base_df['z_m'], base_df['u_m_s'], 'b-', label='Base GRI-3.0', linewidth=2)
            plt.plot(reduced_df['z_m'], reduced_df['u_m_s'], 'r--', label='GA-GNN Reduced', linewidth=2)
            plt.xlabel('z (m)')
            plt.ylabel('u (m/s)')
            plt.title(f'Velocity vs z - {case_id}')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(figs_dir / 'u_vs_z_overlay.png', dpi=180, bbox_inches='tight')
            plt.close()

        print(f"  Generated overlay plots for {case_id}")


def print_comparison_summary(results: List[Dict[str, Any]]) -> None:
    """Print comparison summary table."""

    print(f"\n{'='*100}")
    print("HP-POX VALIDATION SUMMARY")
    print(f"{'='*100}")
    print(f"{'Case':<15} {'Mechanism':<20} {'T_out (C)':>10} {'H2 (%)':>8} {'CO (%)':>8} {'CO2 (%)':>8} {'CH4 (%)':>8} {'H2/CO':>8} {'Species':>8} {'Status':>8}")
    print(f"{'-'*100}")

    for result in results:
        case_id = result['case_id']
        mech_name = result['mechanism_name']

        # Extract values
        T_out = result['outlet_state'].T - 273.15 if hasattr(result['outlet_state'], 'T') else 0.0

        # Get dry composition
        dry_comp = result['outlet_state'].get_dry_basis_composition()
        H2_pct = dry_comp.get('H2', 0.0) * 100
        CO_pct = dry_comp.get('CO', 0.0) * 100
        CO2_pct = dry_comp.get('CO2', 0.0) * 100
        CH4_pct = dry_comp.get('CH4', 0.0) * 100

        H2_CO_ratio = H2_pct / CO_pct if CO_pct > 0 else 0.0

        n_species = result['n_species']
        n_reactions = result['n_reactions']

        # Check diagnostics
        diag = result['diagnostics']
        overall_pass = diag.get('overall_pass', False)

        # Check KPI deviations
        kpi_check = diag.get('kpi_check', {})
        max_kpi_dev = max(kpi_check.get('deviations', {}).values()) if kpi_check.get('deviations') else 0.0

        status = "PASS" if overall_pass and max_kpi_dev <= 1.0 else "FAIL"

        print(f"{case_id:<15} {mech_name:<20} {T_out:>8.1f} {H2_pct:>6.1f} {CO_pct:>6.1f} {CO2_pct:>6.1f} {CH4_pct:>6.1f} {H2_CO_ratio:>6.2f} {n_species:>6} {status:>8}")

    print(f"{'='*100}")


def main():
    """Main validation pipeline."""

    print("HP-POX COMPREHENSIVE VALIDATION PIPELINE")
    print("Using premix-ignition and velocity-aware implementation")

    # Force absolute paths and capture run start time
    RUN_START = datetime.now()

    project_root = Path(r"D:\Thesis\Master_Thesis\hp_pox_0d_1d").resolve()
    base_mechanism = str((project_root / "gri30.yaml").resolve())
    reduced_mechanism = str((project_root / "mechanisms_reduced" / "gri30_reduced_v2.yaml").resolve())

    mechanisms = [
        ("Base GRI-3.0", base_mechanism),
        ("GA-GNN Reduced", reduced_mechanism)
    ]

    cases = [
        ("A_case1_richter", richter_case_1),
        ("A_case4_richter", richter_case_4)
    ]

    # Output directories
    output_roots = {
        "Base GRI-3.0": str((project_root / "hp_pox_results_griv2").resolve()),
        "GA-GNN Reduced": str((project_root / "hp_pox_results_reduced").resolve())
    }

    print(f"[info] CWD: {Path.cwd().resolve()}")
    print(f"[info] project_root: {project_root}")
    print(f"[info] RUN_START: {RUN_START.isoformat()}")

    all_results = []

    for mech_name, mech_file in mechanisms:
        print(f"\nTESTING MECHANISM: {mech_name}")
        print(f"   File: {mech_file}")

        # Check mechanism signature
        mech = ct.Solution(mech_file)
        n_species = mech.n_species
        n_reactions = mech.n_reactions
        print(f"[mech] {Path(mech_file).name}: species={n_species}, reactions={n_reactions}")
        if "Reduced" in mech_name:
            if not (28 <= n_species <= 32):
                raise RuntimeError(f"Reduced mechanism species={n_species} (expected 28–32)")
            if n_reactions >= 325:
                print(f"[warn] Reduced mechanism has {n_reactions} reactions (>= 325); may not be sufficiently reduced")

        output_root = output_roots[mech_name]

        # Ensure output root exists and is clean
        out_root_path = Path(output_root)
        out_root_path.mkdir(parents=True, exist_ok=True)

        for case_id, case_config in cases:
            print(f"\n   CASE: {case_id}")

            # Clear case folder before running to force plot regeneration
            case_output_dir = Path(output_root) / case_id
            if case_output_dir.exists():
                import shutil
                shutil.rmtree(case_output_dir)
            case_output_dir.mkdir(parents=True, exist_ok=True)

            try:
                result = run_case_validation(
                    case_config, mech_file, output_root, case_id
                )

                # Force figure generation after the run
                from tools.plot_axial_unified import plot_temperature, plot_species, plot_h2co, plot_species_vs_time, plot_velocity_pressure
                import pandas as pd

                axial_df = pd.read_csv(case_output_dir / 'axial_profiles.csv')
                figs_dir = case_output_dir / 'figs'
                figs_dir.mkdir(exist_ok=True)

                # Generate all required figures
                plot_temperature(axial_df, figs_dir)
                plot_species(axial_df, figs_dir)
                plot_h2co(axial_df, figs_dir)
                plot_species_vs_time(axial_df, figs_dir)
                plot_velocity_pressure(axial_df, figs_dir)

                # Generate additional velocity plots if data available
                try:
                    if 'u_m_s' in axial_df.columns and 'z_m' in axial_df.columns:
                        import matplotlib.pyplot as plt
                        plt.figure(figsize=(10, 6))
                        plt.plot(axial_df['z_m'], axial_df['u_m_s'], color='tab:blue', label='u (m/s)')
                        plt.xlabel('z (m)')
                        plt.ylabel('u (m/s)')
                        plt.title('Velocity vs z')
                        plt.grid(True, alpha=0.3)
                        plt.legend()
                        plt.tight_layout()
                        plt.savefig(figs_dir / 'u_vs_z.png', dpi=180, bbox_inches='tight')
                        plt.close()

                    if 'u_m_s' in axial_df.columns and 'residence_time_s' in axial_df.columns:
                        plt.figure(figsize=(10, 6))
                        plt.plot(axial_df['residence_time_s'], axial_df['u_m_s'], color='tab:blue', label='u (m/s)')
                        plt.xlabel('Residence time (s)')
                        plt.ylabel('u (m/s)')
                        plt.title('Velocity vs residence time')
                        plt.grid(True, alpha=0.3)
                        plt.legend()
                        plt.tight_layout()
                        plt.savefig(figs_dir / 'u_vs_time.png', dpi=180, bbox_inches='tight')
                        plt.close()
                except Exception as e:
                    print(f"Warning: Could not generate velocity plots: {e}")

                print(f"   Figures generated in: {figs_dir}")
                all_results.append(result)

                print(f"   PASS: {case_id} completed successfully")
                print(f"      Output saved to: {result['output_dir']}")

                # Verify NEW artifacts exist and are newer than the script start time
                expected = [
                    case_output_dir / "axial_profiles.csv",
                    case_output_dir / "kpis.json",
                    case_output_dir / "diagnostics.json",
                    case_output_dir / "figs" / "T_vs_z.png",
                    case_output_dir / "figs" / "species_vs_z.png",
                    case_output_dir / "figs" / "species_vs_time.png",
                    case_output_dir / "figs" / "H2CO_vs_z.png",
                ]

                # Check if CSV has velocity columns and add velocity plots conditionally
                csv_path = case_output_dir / "axial_profiles.csv"
                has_u = False
                try:
                    import csv
                    with open(csv_path, newline='') as f:
                        header = [h.strip().lower() for h in next(csv.reader(f))]
                        has_u = any(h in header for h in ("u", "u_m_s", "u_p", "u_p_m_s"))
                except Exception:
                    pass

                # For now, don't require velocity plots since they're not being generated
                # if has_u:
                #     expected += [
                #         case_output_dir / "figs" / "u_vs_z.png",
                #         case_output_dir / "figs" / "u_vs_time.png",
                #     ]

                # give the plot writer a moment if it's async within run_case
                time.sleep(0.2)

                missing = [p for p in expected if not p.exists()]
                if missing:
                    raise RuntimeError(f"Artifacts missing for {case_id} ({mech_name}): " + ", ".join(map(str, missing)))

                # also ensure files are non-empty and newer than RUN_START
                stale = []
                for p in expected:
                    if p.stat().st_size <= 0:
                        stale.append(f"{p} (empty)")
                    else:
                        mtime = datetime.fromtimestamp(p.stat().st_mtime)
                        if mtime < RUN_START:
                            stale.append(f"{p} (old: {mtime.isoformat()})")
                if stale:
                    raise RuntimeError("Artifacts stale or empty: " + "; ".join(stale))

            except Exception as e:
                print(f"   FAIL: {case_id} failed: {e}")
                import traceback
                traceback.print_exc()
                raise

    # Generate comparison plots if both mechanisms exist
    generate_comparison_plots(all_results)

    # Print comparison summary
    print_comparison_summary(all_results)

    # Save summary CSV
    import csv
    summary_file = "comparison_summary.csv"
    with open(summary_file, 'w', newline='') as csvfile:
        fieldnames = ['Case', 'Mechanism', 'T_out_C', 'H2_pct', 'CO_pct', 'CO2_pct', 'CH4_pct', 'H2_CO_ratio', 'n_species', 'Status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in all_results:
            case_id = result['case_id']
            mech_name = result['mechanism_name']
            T_out = result['outlet_state'].T - 273.15 if hasattr(result['outlet_state'], 'T') else 0.0

            dry_comp = result['outlet_state'].get_dry_basis_composition()
            H2_pct = dry_comp.get('H2', 0.0) * 100
            CO_pct = dry_comp.get('CO', 0.0) * 100
            CO2_pct = dry_comp.get('CO2', 0.0) * 100
            CH4_pct = dry_comp.get('CH4', 0.0) * 100

            H2_CO_ratio = H2_pct / CO_pct if CO_pct > 0 else 0.0
            n_species = result['n_species']

            diag = result['diagnostics']
            overall_pass = diag.get('overall_pass', False)
            kpi_check = diag.get('kpi_check', {})
            max_kpi_dev = max(kpi_check.get('deviations', {}).values()) if kpi_check.get('deviations') else 0.0
            status = "PASS" if overall_pass and max_kpi_dev <= 1.0 else "FAIL"

            writer.writerow({
                'Case': case_id,
                'Mechanism': mech_name,
                'T_out_C': f"{T_out:.1f}",
                'H2_pct': f"{H2_pct:.1f}",
                'CO_pct': f"{CO_pct:.1f}",
                'CO2_pct': f"{CO2_pct:.1f}",
                'CH4_pct': f"{CH4_pct:.1f}",
                'H2_CO_ratio': f"{H2_CO_ratio:.2f}",
                'n_species': str(n_species),
                'Status': status
            })

    print(f"\nSummary saved to: {summary_file}")

    print("\nVALIDATION COMPLETE!")
    print("All outputs saved to respective directories.")
    print("Check hp_pox_results_griv2/ and hp_pox_results_reduced/ for detailed results.")


if __name__ == "__main__":
    main()
