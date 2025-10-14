#!/usr/bin/env python3
"""
Test the enhanced PSR robustness implementation.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from hp_pox_physics.hp_pox_model import HPPOXModel

def test_psr_robustness():
    """Test the enhanced PSR implementation."""
    print("Testing enhanced PSR robustness implementation...")
    
    # Load model with GRI-3.0
    model = HPPOXModel("gri30.yaml", transport_model='mixture-averaged')
    
    # Load inlet streams
    inlet_streams = model.load_inlet_streams('config.yaml')
    
    # Test one case to verify PSR improvements
    case_id = 'A_case1_richter'
    print(f"\nTesting case: {case_id}")
    
    try:
        # Run the case
        results = model.run_case(case_id, inlet_streams, Path('outputs/psr_test'))
        
        # Check results
        print(f"\nPSR Results:")
        print(f"  Converged: {results['psr_diagnostics']['converged']}")
        print(f"  T_out: {results['psr_diagnostics']['T_out']:.1f}K")
        print(f"  Q_chem: {results['psr_diagnostics']['Q_chem']/1000:.2f} kW")
        print(f"  Branch: {results['psr_diagnostics'].get('branch_status', 'UNKNOWN')}")
        print(f"  Residual: {results['psr_diagnostics']['residual']:.2e}")
        
        # Validate PSR
        psr_validation = results['psr_validation']
        print(f"\nPSR Validation:")
        print(f"  Temperature gate: {psr_validation['psr_temperature_gate']}")
        
        if results['psr_diagnostics']['converged'] and psr_validation['psr_temperature_gate']:
            print("✅ PSR robustness test PASSED")
        else:
            print("❌ PSR robustness test FAILED")
            
    except Exception as e:
        print(f"❌ PSR test failed with error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_psr_robustness()
