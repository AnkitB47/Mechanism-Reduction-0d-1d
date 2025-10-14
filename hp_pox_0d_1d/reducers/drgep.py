"""
DRGEP (Directed Relation Graph with Error Propagation) implementation.

Based on Lu & Law (Combustion and Flame, 2005/2006).
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import numpy as np
import cantera as ct

from . import ReductionConfig


def run_drgep(
    base_yaml: Path,
    training_states: List,
    tau_candidates: List[float],
    protected_species: List[str],
    stop_when_species_leq: int,
    out_dir: Path
) -> Tuple[List[int], Path]:
    """
    Run DRGEP algorithm on the training dataset.

    Args:
        base_yaml: Path to base mechanism YAML
        training_states: List of State objects for training
        tau_candidates: List of tau_s values to try
        protected_species: Species that must be kept
        stop_when_species_leq: Stop when species count <= this value
        out_dir: Output directory

    Returns:
        Tuple of (species_mask, reduced_yaml_path)
    """

    print("Running DRGEP (Directed Relation Graph with Error Propagation)...")

    # Load base mechanism
    base_mech = ct.Solution(str(base_yaml))
    n_species = len(base_mech.species())
    n_reactions = len(base_mech.reactions())

    # Start with all species
    species_mask = np.ones(n_species, dtype=bool)

    # Protect required species
    protected_indices = set()
    for species_name in protected_species:
        idx = base_mech.species_index(species_name)
        if idx >= 0:
            protected_indices.add(idx)

    print(f"Protected {len(protected_indices)} species: {[base_mech.species()[i].name for i in protected_indices]}")

    # Build species graph from training states
    species_graph = _build_species_graph(base_mech, training_states)

    # Apply DRGEP pruning
    for tau_s in tau_candidates:
        print(f"  Applying DRGEP with tau_s = {tau_s}")

        # Compute species importance scores
        importance_scores = _compute_drgep_importance(base_mech, species_graph, protected_indices)

        # Prune species below threshold
        kept_species = set(protected_indices)
        for i in range(n_species):
            if i in protected_indices:
                continue
            if importance_scores[i] >= tau_s:
                kept_species.add(i)

        # Update mask
        new_mask = np.zeros(n_species, dtype=bool)
        for i in kept_species:
            new_mask[i] = True

        print(f"    Species count: {new_mask.sum()} (threshold {tau_s})")

        # Check if we're close to target
        if new_mask.sum() <= stop_when_species_leq:
            print(f"    Reached target with tau_s = {tau_s}")
            species_mask = new_mask
            break

        species_mask = new_mask

    # Filter reactions based on kept species
    kept_species_names = {base_mech.species()[i].name for i in range(n_species) if species_mask[i]}
    reaction_indices = []

    for i, reaction in enumerate(base_mech.reactions()):
        # Check if all species in this reaction are kept
        rxn_species = set(reaction.reactants.keys()) | set(reaction.products.keys())
        if rxn_species.issubset(kept_species_names):
            reaction_indices.append(i)

    # Write reduced mechanism
    out_dir.mkdir(parents=True, exist_ok=True)
    reduced_yaml = out_dir / "seed.yaml"

    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    from reduction.mechanism_io import write_reduced_mechanism_ct
    write_reduced_mechanism_ct(
        base_yaml, reaction_indices, reduced_yaml, species_mask.tolist(), enforce_band=False
    )

    print(f"DRGEP complete: {species_mask.sum()} species, {len(reaction_indices)} reactions")
    print(f"Seed mechanism: {reduced_yaml}")

    return species_mask.tolist(), reduced_yaml


def _build_species_graph(base_mech: ct.Solution, training_states: List) -> Dict[int, Dict[int, float]]:
    """Build species interaction graph from training data."""

    n_species = len(base_mech.species())
    species_graph = {i: {} for i in range(n_species)}

    # Process training states
    for state in training_states:
        _add_state_to_graph(base_mech, state, species_graph)

    return species_graph


def _add_state_to_graph(base_mech: ct.Solution, state: Dict[str, Any], species_graph: Dict[int, Dict[int, float]]) -> None:
    """Add a single state to the species graph."""

    n_species = len(base_mech.species())

    # Get ROP data (placeholder for now)
    rop = state.get('rop', {})

    # For each reaction, update species couplings
    for rxn_idx, rate in rop.items():
        if rxn_idx < len(base_mech.reactions()):
            reaction = base_mech.reactions()[rxn_idx]

            # Get species involved in this reaction
            rxn_species = set(reaction.reactants.keys()) | set(reaction.products.keys())

            # Update graph edges between these species
            species_indices = [base_mech.species_index(sp) for sp in rxn_species if base_mech.species_index(sp) >= 0]

            # Add bidirectional edges with reaction rate magnitude
            for i in species_indices:
                for j in species_indices:
                    if i != j:
                        edge_weight = abs(rate)
                        species_graph[i][j] = max(species_graph[i].get(j, 0), edge_weight)


def _compute_drgep_importance(base_mech: ct.Solution, species_graph: Dict[int, Dict[int, float]], protected_indices: Set[int]) -> np.ndarray:
    """Compute DRGEP importance scores for all species."""

    n_species = len(base_mech.species())
    importance = np.zeros(n_species)

    # Start from protected species (importance = 1.0)
    for i in protected_indices:
        importance[i] = 1.0

    # Propagate importance through graph
    changed = True
    max_iterations = 10

    for iteration in range(max_iterations):
        if not changed:
            break

        changed = False
        new_importance = importance.copy()

        for i in range(n_species):
            if i in protected_indices:
                continue

            # Compute propagated importance from neighbors
            max_neighbor_importance = 0
            total_weight = 0

            for j, weight in species_graph[i].items():
                max_neighbor_importance = max(max_neighbor_importance, importance[j])
                total_weight += weight

            if total_weight > 0:
                # Error propagation formula from Lu & Law
                propagated = max_neighbor_importance * (1 - np.exp(-total_weight / 1000))
                new_importance[i] = propagated

                if abs(new_importance[i] - importance[i]) > 1e-6:
                    changed = True

        importance = new_importance

    return importance
