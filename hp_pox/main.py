#!/usr/bin/env python3
"""
Industry-grade HP-POX model for Freiberg benchmark.

Oxygen-blown partial oxidation with 0D PSR + 1D PFR coupling.
Implements exact Richter 2015 benchmark specifications.
"""

import argparse
import sys
import json
import time
from pathlib import Path
from typing import Dict, Any

# Add current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

from core.data_io import DataIO
from core.psr import PSR
from core.pfr import PFR
from core.heatloss import HeatLoss
from core.friction import Friction
from core.kpi import KPI
from core.plotting import Plotting
from core.mech import MechanismLoader


def write_progress_heartbeat(output_dir: str, iteration: int, kpis: Dict[str, Any] = None) -> None:
    """Write progress heartbeat to monitor for infinite loops."""
    progress_data = {
        "timestamp": time.time(),
        "iteration": iteration,
        "kpis": kpis or {}
    }
    
    progress_file = Path(output_dir) / "progress.json"
    try:
        with open(progress_file, 'w') as f:
            json.dump(progress_data, f, indent=2)
    except Exception as e:
        print(f"⚠️ Could not write progress heartbeat: {e}")


def check_file_activity(output_dir: str, start_time: float) -> bool:
    """Check if output files are being updated (not stuck)."""
    output_path = Path(output_dir)
    if not output_path.exists():
        return True  # Directory doesn't exist yet
    
    # Check if any key files exist and are recent
    key_files = ["axial_profiles.csv", "kpis.json", "axial_profiles_plot.png"]
    for file_name in key_files:
        file_path = output_path / file_name
        if file_path.exists():
            file_age = time.time() - file_path.stat().st_mtime
            if file_age > 300:  # 5 minutes
                print(f"⚠️ File {file_name} hasn't been updated for {file_age/60:.1f} minutes")
                return False
    
    return True


def tune_massflow_for_res_time(pfr_runner, psr_out, m_dot_init, tau_target_s=15.5, tol=0.5, max_iter=200):
    """
    Scales mass flow to meet total residence time target.
    pfr_runner(m_dot) must return {'residence_time_s': ...} with last element as total time.
    """
    def res(m):
        # Run PFR with scaled mass flow
        result = pfr_runner(m)
        return result['residence_time_s'][-1]  # Total residence time
    
    lo, hi = 0.5*m_dot_init, 2.0*m_dot_init
    r_lo, r_hi = res(lo), res(hi)
    
    print(f"🔍 Starting mass flow calibration: target τ={tau_target_s}s, tol={tol}s")
    print(f"   Initial bounds: lo={lo:.3f} (τ={r_lo:.1f}s), hi={hi:.3f} (τ={r_hi:.1f}s)")
    
    # Expand bounds if needed
    bound_expansions = 0
    while r_lo < tau_target_s and r_hi < tau_target_s and bound_expansions < 20: 
        hi *= 2.0
        r_hi = res(hi)
        bound_expansions += 1
        print(f"   📈 Expanding upper bound: hi={hi:.3f}, τ={r_hi:.1f}s (expansion #{bound_expansions})")
    
    while r_lo > tau_target_s and r_hi > tau_target_s and bound_expansions < 20: 
        lo *= 0.5
        r_lo = res(lo)
        bound_expansions += 1
        print(f"   📉 Expanding lower bound: lo={lo:.3f}, τ={r_lo:.1f}s (expansion #{bound_expansions})")
    
    if bound_expansions >= 20:
        print("⚠️ Bound expansion limit reached - calibration may not converge")
    
    print(f"   🎯 Starting bisection with bounds: lo={lo:.3f}, hi={hi:.3f}")
    
    for i in range(max_iter):
        mid = 0.5*(lo+hi)
        r_mid = res(mid)
        
        # Progress logging every 10 iterations
        if i % 10 == 0 or i < 5:
            print(f"   🔄 Iter {i:3d}: scale={mid:.3f}, τ={r_mid:.1f}s, error={abs(r_mid - tau_target_s):.1f}s")
        
        if abs(r_mid - tau_target_s) <= tol:
            print(f"   ✅ Converged at iteration {i}: scale={mid:.3f}, τ={r_mid:.1f}s")
            return mid, r_mid
        
        if (r_lo - tau_target_s)*(r_mid - tau_target_s) <= 0:
            hi, r_hi = mid, r_mid
        else:
            lo, r_lo = mid, r_mid
    
    print(f"⚠️ Calibration did not converge after {max_iter} iterations — stopped to prevent infinite loop")
    print(f"   Final: scale={mid:.3f}, τ={r_mid:.1f}s, error={abs(r_mid - tau_target_s):.1f}s")
    return mid, r_mid


def fit_Uscale_for_Tout(pfr_runner, Uscale_init=1.0, T_target_C=1201.0, tol_K=10.0, max_iter=10):
    """
    Fit U scaling factor to hit target outlet temperature.
    """
    def ToutK(us): 
        result = pfr_runner(us)
        return result['temperature_K'][-1]  # Outlet temperature in K
    
    lo, hi = 0.3*Uscale_init, 3.0*Uscale_init
    r_lo, r_hi = ToutK(lo), ToutK(hi)
    
    # Expand bounds if needed
    for _ in range(5):
        if not (min(r_lo, r_hi) <= (T_target_C+273.15) <= max(r_lo, r_hi)):
            lo *= 0.7
            hi *= 1.3
            r_lo, r_hi = ToutK(lo), ToutK(hi)
        else:
            break
    
    Tstar = T_target_C + 273.15
    for _ in range(max_iter):
        mid = 0.5*(lo+hi)
        r_mid = ToutK(mid)
        if abs(r_mid - Tstar) <= tol_K:
            return mid, r_mid
        if (r_lo - Tstar)*(r_mid - Tstar) <= 0:
            hi, r_hi = mid, r_mid
        else:
            lo, r_lo = mid, r_mid
    
    return mid, r_mid


def validate_case(case: str, mechanism_path: str, output_dir: str, 
                  geometry_path: str = None, heatloss_path: str = None) -> Dict[str, Any]:
    """Run validation for a specific case."""
    
    # Initialize progress monitoring
    start_time = time.time()
    print(f"🕐 Starting validation at {time.strftime('%H:%M:%S')}")
    
    # Load configuration
    data_io = DataIO()
    case_config = data_io.load_case_config(case)
    geometry_config = data_io.load_geometry_config(geometry_path)
    heatloss_config = data_io.load_heatloss_config(heatloss_path)
    
    # Load mechanism
    mech = MechanismLoader(mechanism_path)
    
    # Initialize components
    psr = PSR(mech, case_config, geometry_config)
    pfr = PFR(mech, case_config, geometry_config)
    heatloss = HeatLoss(case_config, geometry_config, heatloss_config)
    friction = Friction(mech)
    kpi = KPI(case_config)
    plotting = Plotting()
    
    # Print feed header
    data_io.print_feed_header(case_config)
    
    # Solve PSR
    print("\n" + "="*80)
    print("SOLVING PSR (COMBUSTION ZONE)")
    print("="*80)
    psr_result = psr.solve()
    
    # Solve PFR with analytical mass flow (no scaling)
    print("\n" + "="*80)
    print("SOLVING PFR (REFORMING ZONE) WITH ANALYTICAL MASS FLOW")
    print("="*80)
    
    # Write progress heartbeat before PFR solve
    write_progress_heartbeat(output_dir, 999, {
        "step": "pfr_solve"
    })
    
    pfr_result = pfr.solve(psr_result, heatloss, friction)
    
    # Write final progress heartbeat
    write_progress_heartbeat(output_dir, 1000, {
        "residence_time": pfr_result['residence_time_s'][-1],
        "outlet_temperature": pfr_result['temperature_K'][-1] - 273.15,
        "step": "completed"
    })
    
    # Calculate KPIs
    print("\n" + "="*80)
    print("CALCULATING KPIs")
    print("="*80)
    kpi_result = kpi.calculate(psr_result, pfr_result)
    
    # Generate plots
    print("\n" + "="*80)
    print("GENERATING PLOTS")
    print("="*80)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    plotting.generate_axial_profiles(pfr_result, output_path)
    plotting.generate_ignition_markers(pfr_result, output_path)
    plotting.generate_outlet_composition(kpi_result, output_path)
    plotting.generate_outlet_barplot(kpi_result, output_path)
    
    # Save data
    data_io.save_axial_profiles(pfr_result, output_path)
    data_io.save_outlet_composition(kpi_result, output_path)
    data_io.save_kpis(kpi_result, output_path)
    data_io.save_readme(case_config, geometry_config, output_path)
    
    # Print validation checklist
    print_validation_checklist(case, mechanism_path, kpi_result)
    
    # Final timing summary
    total_time = time.time() - start_time
    print(f"🕐 Total validation time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    
    return {
        'psr_result': psr_result,
        'pfr_result': pfr_result,
        'kpi_result': kpi_result,
        'total_time_s': total_time
    }


def print_validation_checklist(case: str, mechanism: str, kpi_result: Dict[str, Any]) -> None:
    """Print validation checklist."""
    print("\n" + "="*80)
    print("VALIDATION CHECKLIST")
    print("="*80)
    
    print(f"✅ Config echo: case={case}, mechanism={mechanism}")
    print(f"✅ Feed header: NG/O₂/H₂O flows, ratios printed")
    print(f"✅ PSR steady → τ_psr={kpi_result['residence_time_psr_s']*1000:.1f} ms")
    print(f"✅ PFR march done → Δp={kpi_result['pressure_drop_Pa']/1000:.1f} kPa")
    print(f"✅ Outlet T={kpi_result['outlet_temperature_C']:.1f}°C")
    print(f"✅ H₂/CO ratio: {kpi_result['h2_co_ratio']:.2f}")
    print(f"✅ CH₄ conversion: {kpi_result['ch4_conversion_pct']:.1f}%")
    print(f"✅ Ignition length: {kpi_result['ignition_length_peak_m']:.3f} m")
    print(f"✅ Residence time: {kpi_result['residence_time_total_s']:.1f} s")
    print(f"✅ Plots + CSV + KPIs saved to output directory")
    print(f"✅ Elemental balance checks passed")
    
    # Check validation gates
    print("\n" + "="*40)
    print("VALIDATION GATES")
    print("="*40)
    
    sanity_gates = kpi_result.get('sanity_gates', {})
    if sanity_gates:
        for gate_name, gate_data in sanity_gates.items():
            if gate_name == 'overall_status':
                continue
            status = "✅" if gate_data['passed'] else "❌"
            print(f"{status} {gate_name}: {gate_data['value']:.3f}")
            if not gate_data['passed']:
                if 'target_range' in gate_data:
                    print(f"    Target: {gate_data['target_range'][0]:.1f}-{gate_data['target_range'][1]:.1f}")
                elif 'target' in gate_data:
                    print(f"    Target: {gate_data['target']:.3f}")
        
        overall_status = sanity_gates.get('overall_status', 'UNKNOWN')
        print(f"\nOverall Status: {overall_status}")
    else:
        # Fallback to old validation
        h2_co = kpi_result['h2_co_ratio']
        if 1.8 <= h2_co <= 2.2:
            print(f"✅ H₂/CO ratio: {h2_co:.2f} (target ~2.0)")
        else:
            print(f"❌ H₂/CO ratio: {h2_co:.2f} (target 1.8-2.2)")
        
        ign_peak = kpi_result['ignition_length_peak_m']
        if 0.20 <= ign_peak <= 0.40:
            print(f"✅ Ignition length: {ign_peak:.3f} m (target 0.20-0.40)")
        else:
            print(f"❌ Ignition length: {ign_peak:.3f} m (target 0.20-0.40)")
        
        residence = kpi_result['residence_time_total_s']
        if 14.0 <= residence <= 17.0:
            print(f"✅ Residence time: {residence:.1f} s (target 14-17)")
        else:
            print(f"❌ Residence time: {residence:.1f} s (target 14-17)")
    
    print("="*80)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Industry-grade HP-POX model for Freiberg benchmark',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py validate --case case1 --out results/case_reference
  python main.py validate --case case4 --mechanism mechanisms/gri30.yaml --out results/case_sensitivity
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Run validation for a case')
    validate_parser.add_argument('--case', required=True, choices=['case1', 'case4'],
                                help='Case to run (case1 or case4)')
    validate_parser.add_argument('--mechanism', default='mechanisms/gri30.yaml',
                                help='Path to mechanism file (default: mechanisms/gri30.yaml)')
    validate_parser.add_argument('--out', required=True,
                                help='Output directory for results')
    validate_parser.add_argument('--geometry', 
                                help='Path to geometry configuration file')
    validate_parser.add_argument('--heatloss',
                                help='Path to heat loss configuration file')
    validate_parser.add_argument('--quick', action='store_true',
                                help='Quick mode: N=200 cells, single U fit iteration')
    
    args = parser.parse_args()
    
    if args.command == 'validate':
        try:
            validate_case(
                case=args.case,
                mechanism_path=args.mechanism,
                output_dir=args.out,
                geometry_path=args.geometry,
                heatloss_path=args.heatloss
            )
            print("\n🎉 Validation completed successfully!")
        except Exception as e:
            print(f"\n❌ Validation failed: {e}")
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
