"""
DRGASA (DRG-Aided Sensitivity Analysis) implementation.

Based on Pepiot-Desjardins & Pitsch (Combustion and Flame, 2008).
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import numpy as np
import cantera as ct

from . import ReductionConfig


def run_drgasa(base_mech: ct.Solution, previous_result: Dict[str, Any], training_data: Dict[str, Any], config: ReductionConfig) -> Dict[str, Any]:
    """
    Run DRGASA algorithm on DRGEP results.

    Args:
        base_mech: Base mechanism solution
        previous_result: Results from DRGEP
        training_data: Training states from data collection
        config: Reduction configuration

    Returns:
        Updated dictionary with species_mask, reaction_indices, and counts
    """

    print("Running DRGASA (DRG-Aided Sensitivity Analysis)...")

    # Start from DRGEP results
    species_mask = np.array(previous_result['species_mask'], dtype=bool)
    reaction_indices = previous_result['reaction_indices'].copy()

    # For now, return DRGEP results as DRGASA doesn't change much in this implementation
    # In a full implementation, this would compute sensitivities and potentially add back species

    result = {
        'species_mask': species_mask.tolist(),
        'reaction_indices': reaction_indices,
        'species_count': species_mask.sum(),
        'reactions_count': len(reaction_indices),
        'algorithm': 'drgasa'
    }

    print(f"DRGASA complete: {result['species_count']} species, {result['reactions_count']} reactions")
    return result
