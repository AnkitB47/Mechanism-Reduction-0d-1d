#!/usr/bin/env python3
"""
Debug single case execution.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from hp_pox_physics.hp_pox_model import HPPOXModel
from hp_pox_physics.config_cases import richter_case_1
import cantera as ct


def debug_case():
    print("=== DEBUGGING CASE EXECUTION ===")

    # Initialize model
    model = HPPOXModel(mechanism_path="gri30.yaml")
    print(f"Loaded mechanism: {len(model.thermo.species_names)} species, {len(model.thermo.gas.reaction_equations())} reactions")

    # Create inlet configuration
    inlet_config = {
        'case1': {
            'primary_steam': {
                'mass_kg_per_h': richter_case_1.primary_steam.mass_kg_per_s * 3600,
                'T_C': richter_case_1.primary_steam.T_K - 273.15,
                'comp': {'H2O': 1.0}
            },
            'secondary_steam': {
                'mass_kg_per_h': richter_case_1.secondary_steam.mass_kg_per_s * 3600,
                'T_C': richter_case_1.secondary_steam.T_K - 273.15,
                'comp': {'H2O': 1.0}
            },
            'oxygen': {
                'mass_kg_per_h': richter_case_1.oxygen.mass_kg_per_s * 3600,
                'T_C': richter_case_1.oxygen.T_K - 273.15,
                'comp': {'O2': 1.0}
            },
            'nitrogen_optisox': {
                'mass_kg_per_h': richter_case_1.nitrogen.mass_kg_per_s * 3600,
                'T_C': richter_case_1.nitrogen.T_K - 273.15,
                'comp': {'N2': 1.0}
            },
            'natural_gas': {
                'mass_kg_per_h': richter_case_1.natural_gas.mass_kg_per_s * 3600,
                'T_C': richter_case_1.natural_gas.T_K - 273.15,
                'composition_volpct': richter_case_1.natural_gas.composition_X
            }
        }
    }

    try:
        print(f"Running case: {richter_case_1.case_id}")

        # Create inlet streams for PSR
        inlet_streams = [
            {
                'name': 'primary_steam',
                'T': richter_case_1.primary_steam.T_K,
                'm_dot': richter_case_1.primary_steam.mass_kg_per_s,
                'Y': richter_case_1.primary_steam.composition_X
            },
            {
                'name': 'secondary_steam',
                'T': richter_case_1.secondary_steam.T_K,
                'm_dot': richter_case_1.secondary_steam.mass_kg_per_s,
                'Y': richter_case_1.secondary_steam.composition_X
            },
            {
                'name': 'oxygen',
                'T': richter_case_1.oxygen.T_K,
                'm_dot': richter_case_1.oxygen.mass_kg_per_s,
                'Y': richter_case_1.oxygen.composition_X
            },
            {
                'name': 'nitrogen',
                'T': richter_case_1.nitrogen.T_K,
                'm_dot': richter_case_1.nitrogen.mass_kg_per_s,
                'Y': richter_case_1.nitrogen.composition_X
            },
            {
                'name': 'natural_gas',
                'T': richter_case_1.natural_gas.T_K,
                'm_dot': richter_case_1.natural_gas.mass_kg_per_s,
                'Y': richter_case_1.natural_gas.composition_X
            }
        ]

        # Run the case
        case_out_dir = Path("debug_output") / richter_case_1.case_id
        case_out_dir.mkdir(parents=True, exist_ok=True)

        results = model.run_case(richter_case_1.case_id, inlet_config, str(case_out_dir))

        print(f"Case completed successfully")
        print(f"Results keys: {list(results.keys())}")

        if 'pfr_results_on' in results:
            pfr = results['pfr_results_on']
            print(f"PFR final T: {pfr['temperature_K'][-1] - 273.15:.1f}°C")

            # Check if axial_profiles.csv was created
            axial_file = case_out_dir / 'axial_profiles.csv'
            if axial_file.exists():
                print(f"Axial profiles saved: {axial_file}")
            else:
                print("ERROR: axial_profiles.csv not created")

        return True

    except Exception as e:
        print(f"Case failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    debug_case()
