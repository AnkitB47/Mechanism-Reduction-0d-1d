#!/usr/bin/env python3
"""
Cross-check code/configs against REF and write Table A coverage CSV.
"""

from pathlib import Path
import yaml
from hp_pox_physics import config_reference
from hp_pox_physics import config_cases


def main():
    repo_root = Path(__file__).resolve().parents[1]
    cfg_path = repo_root / 'config.yaml'
    out_dir = repo_root / 'report'
    out_dir.mkdir(parents=True, exist_ok=True)
    table_a_csv = out_dir / 'table_A_coverage.csv'

    with open(cfg_path, 'r') as f:
        config_yaml = yaml.safe_load(f)

    items = config_reference.cross_check(config_yaml, config_cases)
    config_reference.write_table_a(items, str(table_a_csv))
    print(f"Table A written to: {table_a_csv.as_posix()}")


if __name__ == '__main__':
    main()


