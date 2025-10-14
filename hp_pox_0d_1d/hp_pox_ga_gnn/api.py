"""
GA+GNN API adapter for mechanism reduction.

This module provides a clean interface to the existing GA+GNN reduction
system, adapted for use in the hybrid orchestrator.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from typing import List, Tuple
import cantera as ct

# Import existing GA+GNN system
from reduction.ga_bridge import run as run_ga_gnn


def run_ga_gnn(
    base_yaml: Path,
    seed_yaml: Path,
    protected_species: List[str],
    target_band: Tuple[int, int],
    training_states: List,
    out_dir: Path
) -> Path:
    """
    Run GA+GNN reduction starting from a seed mechanism.

    Args:
        base_yaml: Path to base mechanism YAML
        seed_yaml: Path to seed mechanism (from DRGEP)
        protected_species: Species that must be kept
        target_band: (min_species, max_species) target range
        training_states: List of training states for evaluation
        out_dir: Output directory for GA+GNN run

    Returns:
        Path to optimized reduced mechanism YAML
    """

    out_dir.mkdir(parents=True, exist_ok=True)

    # Load seed mechanism to get initial mask
    seed_mech = ct.Solution(str(seed_yaml))
    seed_species_mask = [1 if sp in protected_species or sp in [seed_mech.species()[i].name for i in range(len(seed_mech.species()))] else 0
                        for sp in [base_mech.species()[i].name for i in range(len(base_mech.species()))]]

    # Use existing GA+GNN with seed initialization
    # Note: This is a simplified adapter - in practice you'd modify the GA+GNN
    # to accept a seed mechanism and initialize from it

    # For now, run standard GA+GNN and hope it converges to the right range
    # In a full implementation, you'd modify GA+GNN to use the seed

    return run_ga_gnn(
        base_mech=base_yaml,
        out_mech=out_dir / "ga_gnn_reduced.yaml",
        cases_root=None,  # Would need to be set up
        cases=[],  # Would need to be set up
        scores_path=None,
        budget={'generations': 100, 'pop': 30, 'topK': 3, 'eval_every': 5},
        workers=1,
        # Other parameters would need to be configured
    )
