#!/usr/bin/env python3
"""
HP-POX Model CLI

Clean command-line interface for running HP-POX 0D-1D reactor simulations.
Uses the hp_pox_physics module as the single source of truth.
"""

import argparse
from hp_pox_physics import run_case, run_all


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="HP-POX 0D-1D Reactor Model",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --case A_case1_richter --chemistry on
  python main.py --case A_case1_richter --chemistry off
  python main.py --chemistry on
        """
    )
    
    parser.add_argument(
        "--case", 
        type=str, 
        default="", 
        help="Case ID (e.g., A_case1_richter). Empty = run all cases."
    )
    parser.add_argument(
        "--chemistry", 
        choices=["on", "off"], 
        default="on",
        help="Enable chemistry (on) or wall cooling only (off)"
    )
    parser.add_argument(
        "--outdir", 
        type=str, 
        default="hp_pox_results",
        help="Output directory for results"
    )
    parser.add_argument(
        "--mech", 
        type=str, 
        default="SAN_DIEGO_MECH/c0_mech.yaml",
        help="Chemical mechanism file (default: San Diego for laptop, then Aramco3.0 for HPC)"
    )
    parser.add_argument(
        "--transport",
        type=str,
        default="mix",
        choices=["multi", "mix"],
        help="Transport model: 'multi' for multicomponent, 'mix' for mixture-averaged (default: mix for laptop)"
    )
    parser.add_argument(
        "--procs",
        type=int,
        default=1,
        help="Number of processes for parallel execution (default: 1 for laptop)"
    )
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()
    chemistry_on = args.chemistry == "on"
    transport_model = "multicomponent" if args.transport == "multi" else "mixture-averaged"
    
    print(f"HP-POX 0D-1D Reactor Model")
    print(f"Chemistry: {'ON' if chemistry_on else 'OFF'}")
    print(f"Mechanism: {args.mech}")
    print(f"Transport: {transport_model}")
    print(f"Output directory: {args.outdir}")
    print(f"{'='*60}")
    
    try:
        if args.case:
            # Run single case
            print(f"Running case: {args.case}")
            # Enforce required species for methane POX mechanisms
            from hp_pox_physics.thermo import ThermoManager
            tm = ThermoManager(args.mech, 'mixture-averaged' if transport_model == 'mixture-averaged' else 'multicomponent')
            REQUIRED = ["CO", "CO2", "CH4", "H2", "O2", "H2O", "N2"]
            missing = [s for s in REQUIRED if s not in tm.gas.species_names]
            if missing:
                raise RuntimeError(f"[{args.case}] Mechanism missing required species: {missing}")
            print(f"Mechanism OK for {args.case}: {args.mech} | species={len(tm.gas.species_names)}")

            kpis = run_case(args.case, chemistry=chemistry_on, outdir=args.outdir, mechanism=args.mech, transport=transport_model, procs=args.procs)
            
            print(f"\n{'='*60}")
            print("CASE SUMMARY")
            print(f"{'='*60}")
            print(f"Case: {args.case}")
            print(f"PSR Temperature: {kpis['psr_outlet'].temperature - 273.15:.1f}°C")
            print(f"PFR Outlet Temperature: {kpis['outlet_state'].temperature - 273.15:.1f}°C")
            print(f"H2/CO Ratio: {kpis['h2_co_ratio']:.2f}")
            print(f"Ignition Length: {kpis['pfr_results_on']['ignition_length_peak_m']:.3f}m")
            print(f"Results saved to: {args.outdir}/{args.case}/")
            
        else:
            # Run all cases
            print("Running all available cases")
            summary = run_all(cases=None, chemistry=chemistry_on, outdir=args.outdir, mechanism=args.mech)
            
            print(f"\n{'='*60}")
            print("SUMMARY OF ALL CASES")
            print(f"{'='*60}")
            for case_id, result in summary.items():
                if result['status'] == 'success':
                    print(f"{case_id:20s}: SUCCESS - T_out={result['pfr_outlet_temperature_C']:.1f}°C, H2/CO={result['h2_co_ratio']:.2f}")
                else:
                    print(f"{case_id:20s}: FAILED - {result['error']}")
            
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
