#!/usr/bin/env python3
"""
Sanity test for adaptive chemistry PFR with 10-20 cells.
Tests the new adaptive chemistry implementation with realistic conditions.
"""

import sys
import os
from pathlib import Path

# Add the hp_pox_physics module to path
sys.path.append(str(Path(__file__).parent / "hp_pox_physics"))

from hp_pox_physics.thermo import ThermoManager
from hp_pox_physics.psr import PSR
from hp_pox_physics.pfr import PFR
from hp_pox_physics.mixing import mix_streams

def run_sanity_test():
    """Run a short PFR test with 10-20 cells."""
    print("🧪 Running PFR Sanity Test (10-20 cells)")
    print("=" * 50)
    
    # Test parameters
    dt_min = 1e-7  # Start with larger dt_min for sanity test
    dz_init = 0.1  # 10 cm cells
    target_cells = 15  # 1.5m total length
    
    # Initialize components
    print("Initializing components...")
    thermo = ThermoManager("SAN_DIEGO_MECH/c0_mech.yaml", transport_model="mixture-averaged")
    psr = PSR(16.136e-6, thermo)  # Small volume for sanity test
    pfr = PFR(1.5, 0.49, thermo)  # Short reactor
    
    # Create inlet state (simplified)
    print("Creating inlet state...")
    streams = [
        {
            'm_dot': 0.1,  # kg/s
            'T': 500 + 273.15,  # K
            'Y': {'CH4': 0.96, 'C2H6': 0.02, 'N2': 0.02},
            'P': 51.013e5  # Pa
        },
        {
            'm_dot': 0.3,  # kg/s
            'T': 200 + 273.15,  # K
            'Y': {'O2': 1.0},
            'P': 51.013e5  # Pa
        },
        {
            'm_dot': 0.05,  # kg/s
            'T': 300 + 273.15,  # K
            'Y': {'H2O': 1.0},
            'P': 51.013e5  # Pa
        }
    ]
    
    T_mixed, Y_mixed, P_mixed = mix_streams(streams, thermo)
    
    # Create gas state
    inlet_state = thermo.create_gas_state(T_mixed, P_mixed, Y_mixed)
    print(f"  Inlet T: {inlet_state.temperature:.1f}K")
    print(f"  Inlet P: {inlet_state.pressure/1e5:.1f} bar")
    
    # Run PSR (simplified)
    print("Running PSR...")
    try:
        psr_outlet, converged, psr_diagnostics = psr.solve_psr_continuation(
            inlet_state, 
            m_dot=0.45,  # Total mass flow
            P=51.013e5,  # 51 bar
            q_burner_kw=10.0  # Small burner heat
        )
        print(f"  PSR outlet T: {psr_outlet.temperature:.1f}K")
    except Exception as e:
        print(f"  PSR failed: {e}, using inlet state")
        psr_outlet = inlet_state
    
    # Run PFR with sanity test parameters
    print("Running PFR...")
    try:
        pfr_results = pfr.solve_pfr(
            psr_outlet,
            m_dot=0.45,  # Total mass flow
            P=51.013e5,  # 51 bar
            q_wall_per_m=1000.0,  # 1 kW/m wall loss
            chemistry_on=True
        )
        
        print(f"  PFR completed: {len(pfr_results['z_m'])} cells")
        print(f"  Outlet T: {pfr_results['temperature_K'][-1]:.1f}K")
        
        # Check unified data
        if hasattr(pfr, 'unified_data') and pfr.unified_data:
            print(f"  Unified data: {len(pfr.unified_data)} rows")
            
            # Print summary
            df = pd.DataFrame(pfr.unified_data)
            total_cells = len(df)
            fallback_cells = df['used_fallback'].sum()
            
            print(f"  Chemistry Summary:")
            print(f"    Total cells: {total_cells}")
            print(f"    Fallback cells: {fallback_cells} ({100*fallback_cells/total_cells:.1f}%)")
            print(f"    Mean substeps: {df['substeps'].mean():.1f}")
            print(f"    Mean halvings: {df['halvings'].mean():.1f}")
            
            if 'H2_CO' in df.columns and df['H2_CO'].max() > 0:
                print(f"    H2/CO range: {df['H2_CO'].min():.3f} - {df['H2_CO'].max():.3f}")
        
        # Save unified CSV
        os.makedirs("outputs/pfr_run", exist_ok=True)
        pfr.save_unified_csv("outputs/pfr_run/test_profile.csv")
        
        # Generate plots
        print("Generating plots...")
        os.system("python tools/plot_pfr_unified.py outputs/pfr_run/test_profile.csv --out figures")
        
        print("✅ Sanity test completed successfully!")
        
        # Acceptance checks
        print("\n📋 Acceptance Checks:")
        if hasattr(pfr, 'unified_data') and pfr.unified_data:
            df = pd.DataFrame(pfr.unified_data)
            
            # Check fallback map is not all ones
            fallback_ratio = df['used_fallback'].mean()
            if 0 < fallback_ratio < 1:
                print(f"  ✅ Fallback map mixed (not all ones): {100*fallback_ratio:.1f}% fallback")
            else:
                print(f"  ⚠️  Fallback map uniform: {100*fallback_ratio:.1f}% fallback")
            
            # Check H2/CO variation
            if 'H2_CO' in df.columns and df['H2_CO'].std() > 0.01:
                print(f"  ✅ H2/CO shows variation: std={df['H2_CO'].std():.3f}")
            else:
                print(f"  ⚠️  H2/CO variation limited: std={df['H2_CO'].std():.3f}")
            
            # Check temperature evolution
            if len(df) > 1 and abs(df['T'].iloc[-1] - df['T'].iloc[0]) > 10:
                print(f"  ✅ Temperature evolves: ΔT={abs(df['T'].iloc[-1] - df['T'].iloc[0]):.1f}K")
            else:
                print(f"  ⚠️  Limited temperature evolution: ΔT={abs(df['T'].iloc[-1] - df['T'].iloc[0]):.1f}K")
            
            # Check step sizes
            if df['dt_effective'].std() > 0:
                print(f"  ✅ Step sizes vary: std={df['dt_effective'].std():.2e}")
            else:
                print(f"  ⚠️  Step sizes constant: std={df['dt_effective'].std():.2e}")
            
            # Check convergence metrics
            if df['substeps'].std() > 0 or df['halvings'].std() > 0:
                print(f"  ✅ Convergence metrics vary: substeps_std={df['substeps'].std():.1f}, halvings_std={df['halvings'].std():.1f}")
            else:
                print(f"  ⚠️  Convergence metrics constant")
        
    except Exception as e:
        print(f"❌ PFR failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    import pandas as pd
    success = run_sanity_test()
    sys.exit(0 if success else 1)
