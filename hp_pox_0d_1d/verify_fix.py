#!/usr/bin/env python3
from hp_pox_0d_1d.hp_pox_physics.config_cases import CaseConstants

config = CaseConstants.get_case_config('A_case1_richter')
print('=== FIXED CONFIG VERIFICATION ===')
print(f'Q_burner: {config["burner_heat_loss"]} W (should be 34290 W)')
print(f'Wall heat loss: {config["wall_heat_loss"]} W/m')
print(f'Geometry: {config["geometry"]}')
print()
print('Expected: 34.29 kW = 34290 W')
print('Actual:', config['burner_heat_loss'], 'W')
print('Match:', abs(config['burner_heat_loss'] - 34290) < 1)
