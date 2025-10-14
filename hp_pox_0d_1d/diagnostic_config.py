#!/usr/bin/env python3
"""
Config-driven diagnostic for HP-POX baseline pipeline.
"""

import os
import sys
sys.path.insert(0, '.')

def main():
    # 0) Canonical sources verification
    print('=== CANONICAL SOURCES VERIFICATION ===')
    print('[whoami] config_cases.py:', os.path.abspath('hp_pox_0d_1d/hp_pox_physics/config_cases.py'))
    print('[whoami] gri30.yaml:', os.path.abspath('hp_pox_0d_1d/gri30.yaml'))
    print('[whoami] hp_pox_model.py:', os.path.abspath('hp_pox_0d_1d/hp_pox_physics/hp_pox_model.py'))
    print('[whoami] simulate_baseline.py:', os.path.abspath('hp_pox_0d_1d/tools/simulate_baseline.py'))

    # Import and check what's actually loaded
    import hp_pox_0d_1d.hp_pox_physics.config_cases as cc
    print('[loaded] config_cases module:', os.path.abspath(cc.__file__))

    # Check case constants
    print()
    print('=== CASE CONSTANTS DUMP ===')
    from hp_pox_0d_1d.hp_pox_physics.config_cases import CaseConstants

    cases = ['A_case1_richter', 'A_case4_richter', 'B_case1_134L', 'B_case4_134L']
    for case in cases:
        print()
        print(f'--- {case} ---')
        try:
            config = CaseConstants.get_case_config(case)
            print(f'  Pressure: {config["pressure"]} Pa')
            print(f'  Target T: {config["target_temperature"]} K')
            print(f'  Q_burner: {config["burner_heat_loss"]} W')
            print(f'  Wall heat loss: {config["wall_heat_loss"]} W/m')
            print(f'  Geometry: L={config["geometry"]["length"]}m, D={config["geometry"]["diameter"]}m')
        except Exception as e:
            print(f'  ERROR: {e}')

if __name__ == '__main__':
    main()
