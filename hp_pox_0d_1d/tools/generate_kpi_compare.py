#!/usr/bin/env python3
"""
Aggregate KPIs from hp_pox_results_* folders and compare vs REF outlet targets.
Writes report/kpi_compare_<timestamp>.csv.
"""

from pathlib import Path
import json
import csv
from datetime import datetime
from typing import Dict, Any, Union
from .common_io import load_kpis, dry_value, perr


def build_kpi_row(case: str, base_kpi_path: Union[str, Path], red_kpi_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Return dict with H2, CO, CO2, CH4, H2/CO, T_out_K and % errors.
    
    Args:
        case: Case name
        base_kpi_path: Path to baseline KPIs JSON
        red_kpi_path: Path to reduced KPIs JSON
        
    Returns:
        Dictionary with comparison data
    """
    base_kpis = load_kpis(base_kpi_path)
    red_kpis = load_kpis(red_kpi_path)
    
    if not base_kpis or not red_kpis:
        return {}
    
    # Extract dry basis values
    h2_base = dry_value(base_kpis, 'H2')
    h2_red = dry_value(red_kpis, 'H2')
    co_base = dry_value(base_kpis, 'CO')
    co_red = dry_value(red_kpis, 'CO')
    co2_base = dry_value(base_kpis, 'CO2')
    co2_red = dry_value(red_kpis, 'CO2')
    ch4_base = dry_value(base_kpis, 'CH4')
    ch4_red = dry_value(red_kpis, 'CH4')
    
    # Calculate H2/CO ratios
    h2co_base = h2_base / max(co_base, 1e-12) if co_base > 0 else 0.0
    h2co_red = h2_red / max(co_red, 1e-12) if co_red > 0 else 0.0
    
    # Temperature
    t_base = float(base_kpis.get('pfr_outlet_temperature_C', 0.0)) + 273.15
    t_red = float(red_kpis.get('pfr_outlet_temperature_C', 0.0)) + 273.15
    
    return {
        'case': case,
        'H2_red': h2_red,
        'H2_base': h2_base,
        'H2_err_%': perr(h2_base, h2_red),
        'CO_red': co_red,
        'CO_base': co_base,
        'CO_err_%': perr(co_base, co_red),
        'CO2_red': co2_red,
        'CO2_base': co2_base,
        'CO2_err_%': perr(co2_base, co2_red),
        'CH4_red': ch4_red,
        'CH4_base': ch4_base,
        'CH4_err_%': perr(ch4_base, ch4_red),
        'H2/CO_red': h2co_red,
        'H2/CO_base': h2co_base,
        'H2/CO_err_%': perr(h2co_base, h2co_red),
        'T_out_red_K': t_red,
        'T_out_base_K': t_base,
        'dT_out_K': t_red - t_base
    }


def main():
    """Original CLI functionality - kept for backward compatibility."""
    try:
        from hp_pox_physics.config_reference import REF
    except ImportError:
        print("Warning: REF targets not available, using simplified comparison")
        REF = None
    
    repo_root = Path(__file__).resolve().parents[1]
    out_dir = repo_root / 'report'
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_csv = out_dir / f'kpi_compare_{ts}.csv'

    # enumerate mechanisms/cases
    mech_dirs = [repo_root / 'hp_pox_results_gri', repo_root / 'hp_pox_results_sandiego']
    cases = ['A_case1_richter','A_case4_richter','B_case1_134L','B_case4_134L']

    rows = []
    for mdir in mech_dirs:
        mech = mdir.name.split('_')[-1]
        for case in cases:
            kpis_path = mdir / case / 'kpis.json'
            k = load_kpis(kpis_path)
            if not k:
                continue
            # outlet dry basis fractions if available
            sim_dry = k.get('sim_dry', {})
            h2 = sim_dry.get('H2', 0.0)
            co = sim_dry.get('CO', 0.0)
            co2 = sim_dry.get('CO2', 0.0)
            ch4 = sim_dry.get('CH4', 0.0)
            h2co = (h2 / max(co, 1e-12)) if co > 0 else 0.0
            t_out_k = (k.get('pfr_outlet_temperature_C', 0.0) + 273.15) if 'pfr_outlet_temperature_C' in k else 0.0
            
            # REF targets (if available)
            if REF:
                tref = None
                if 'case1' in case:
                    tref = REF['outlet_targets']['case1']
                elif 'case4' in case:
                    tref = REF['outlet_targets']['case4']
                h2_ref = (tref.get('H2_volpct', 0.0) / 100.0) if tref else 0.0
                co_ref = (tref.get('CO_volpct', 0.0) / 100.0) if tref else 0.0
                co2_ref = (tref.get('CO2_volpct', 0.0) / 100.0) if tref else 0.0
                ch4_ref = (tref.get('CH4_volpct', 0.0) / 100.0) if tref else 0.0
                t_ref_k = (tref.get('T_out_C', 0.0) + 273.15) if tref else 0.0
            else:
                h2_ref = co_ref = co2_ref = ch4_ref = t_ref_k = 0.0

            rows.append({
                'mechanism': mech,
                'case': case,
                'H2_out': h2, 'CO_out': co, 'CO2_out': co2, 'CH4_out': ch4,
                'H2_CO': h2co, 'T_out_K': t_out_k,
                'H2_ref': h2_ref, 'CO_ref': co_ref, 'CO2_ref': co2_ref, 'CH4_ref': ch4_ref, 'T_ref_K': t_ref_k,
                'dH2': h2 - h2_ref, 'dCO': co - co_ref, 'dCO2': co2 - co2_ref, 'dCH4': ch4 - ch4_ref, 'dT_K': t_out_k - t_ref_k
            })

    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"KPI comparison written to: {out_csv.as_posix()}")


if __name__ == '__main__':
    main()


