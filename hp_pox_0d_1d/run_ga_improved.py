#!/usr/bin/env python3
"""
Run GA reduction with improved CH4 pathway protection and fitness function.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from reduction.cli import main

if __name__ == "__main__":
    # Run GA reduction with improved constraints and fitness
    main([
        'reduce',
        '--base-mech', 'gri30.yaml',
        '--constraints', 'reduction/configs/constraints.yaml',
        '--out-mech', 'mechanisms_reduced/gri30_reduced_v3.yaml',
        '--cases-root', 'hp_pox_results_gri',
        '--cases', 'A_case1_richter', 'A_case4_richter', 'B_case1_134L', 'B_case4_134L',
        '--workers', '4',
        '--pop', '40',
        '--gens', '20',
        '--elite', '4',
        '--mut', '0.12',
        '--cx', '0.8'
    ])
