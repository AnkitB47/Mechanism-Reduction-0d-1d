"""
Fixed-point reaction↔species closure implementation.

Ensures consistency between species and reactions by iterative pruning.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import numpy as np
import cantera as ct

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from .__init__ import ReductionConfig
from reduction.mechanism_io import write_reduced_mechanism_ct


def fixed_point_closure(
    yaml_in: Path,
    protected_species: List[str],
    out_yaml: Path,
    max_iter: int = 50,
    verbose: bool = True
) -> Path:
    """
    Run fixed-point reaction↔species closure.

    Args:
        yaml_in: Input mechanism YAML path
        protected_species: Species that must be kept
        out_yaml: Output mechanism YAML path
        max_iter: Maximum iterations
        verbose: Enable verbose logging

    Returns:
        Path to closed mechanism YAML
    """

    if verbose:
        print("Running fixed-point reaction<->species closure...")

    # Load input mechanism
    mech = ct.Solution(str(yaml_in))
    n_species = len(mech.species())
    n_reactions = len(mech.reactions())
    species_names = [mech.species()[i].name for i in range(n_species)]

    # Get protected species indices
    protected_indices = set()
    for species_name in protected_species:
        idx = mech.species_index(species_name)
        if idx >= 0:
            protected_indices.add(idx)

    # Start with all species from input mechanism
    S = set(species_names[i] for i in range(n_species)) | set(protected_species)

    changed = True
    iteration = 0

    while changed and iteration < max_iter:
        changed = False
        iteration += 1

        # Keep only reactions whose (reactants ∪ products) ⊆ S
        kept_reaction_indices = []
        for i in range(n_reactions):
            reaction = mech.reactions()[i]
            rxn_species = set(reaction.reactants.keys()) | set(reaction.products.keys())
            if rxn_species.issubset(S):
                kept_reaction_indices.append(i)

        # Compute live species from kept reactions + protected
        live_species = set()
        for i in kept_reaction_indices:
            reaction = mech.reactions()[i]
            live_species.update(reaction.reactants.keys())
            live_species.update(reaction.products.keys())
        live_species |= set(protected_species)

        # Check if live_species differs from S
        if live_species != S:
            S = live_species
            changed = True

        if verbose:
            print(f"  iter {iteration}: species={len(S)}, reactions={len(kept_reaction_indices)}")

    # Create final species mask
    final_species_mask = np.zeros(n_species, dtype=bool)
    for i in range(n_species):
        if species_names[i] in S:
            final_species_mask[i] = True

    # Write closed mechanism
    from reduction.mechanism_io import write_reduced_mechanism_ct
    write_reduced_mechanism_ct(
        yaml_in, kept_reaction_indices, out_yaml, final_species_mask.tolist(), enforce_band=False
    )

    if verbose:
        print(f"Fixed-point closure complete: {final_species_mask.sum()} species, {len(kept_reaction_indices)} reactions")
        print(f"Closed mechanism: {out_yaml}")

    return out_yaml
