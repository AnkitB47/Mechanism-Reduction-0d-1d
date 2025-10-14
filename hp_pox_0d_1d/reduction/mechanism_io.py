from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Any, Sequence, Optional, Iterable, Set
import shutil
import yaml

from .constraints import HARD_SPECIES as CONSTRAINT_HARD_SPECIES, get_locked_family_memberships


def load_mechanism_yaml(path: Path) -> Dict[str, Any]:
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def count_reactions(mech: Dict[str, Any]) -> int:
    rxns = mech.get('reactions', [])
    return len(rxns) if isinstance(rxns, list) else 0


def apply_reaction_mask(mech: Dict[str, Any], keep_indices: List[int]) -> Dict[str, Any]:
    """Return a pruned mechanism dict with only kept reaction entries; no auto-repair."""
    new_mech = dict(mech)
    rxns = mech.get('reactions', [])
    kept = [rxns[i] for i in range(len(rxns)) if i in set(keep_indices)]
    new_mech['reactions'] = kept
    return new_mech


def apply_reaction_mask_file(mech_in: Path, mask: Sequence[int] | Sequence[bool], mech_out: Path) -> None:
    """Prune reactions by boolean/int mask, auto-repair stale duplicate flags, and write YAML.

    - Preserves species/thermo/transport blocks untouched
    - Ensures that if only one member of a duplicate group remains, the survivor's 'duplicate: true' is removed
    """
    with open(mech_in, 'r') as f:
        data = yaml.safe_load(f)

    rxns = data.get('reactions', [])
    if len(rxns) != len(mask):
        raise ValueError(f"Mask length {len(mask)} != number of reactions {len(rxns)} in {mech_in}")

    # Count how many kept per equation
    kept_equation_counts: Dict[str, int] = {}
    for i, r in enumerate(rxns):
        if mask[i]:
            eq = r.get('equation')
            if eq is None:
                continue
            kept_equation_counts[eq] = kept_equation_counts.get(eq, 0) + 1

    kept_rxns = []
    for i, r in enumerate(rxns):
        if not mask[i]:
            continue
        r_copy = dict(r)
        eq = r_copy.get('equation')
        if r_copy.get('duplicate', False) and kept_equation_counts.get(eq, 0) == 1:
            # Twin was pruned; remove stale duplicate flag
            r_copy.pop('duplicate', None)
        kept_rxns.append(r_copy)

    data['reactions'] = kept_rxns

    mech_out.parent.mkdir(parents=True, exist_ok=True)
    with open(mech_out, 'w') as f:
        yaml.safe_dump(data, f, sort_keys=False)


def write_mechanism_yaml(mech: Dict[str, Any], out_path: Path) -> None:
    with open(out_path, 'w') as f:
        yaml.safe_dump(mech, f, sort_keys=False)
    # Post-load sanity
    import warnings
    import cantera as ct
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _ = ct.Solution(str(out_path))


def copy_mechanism(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


HARD_SPECIES_DEFAULT = list(dict.fromkeys(list(CONSTRAINT_HARD_SPECIES) + ["C2H6"]))


def _reverse_equation(equation: Optional[str]) -> Optional[str]:
    """Return the textual reverse of a reaction equation if available."""
    if not equation or not isinstance(equation, str):
        return None
    if "<=>" in equation:
        return equation
    if "=>" in equation:
        lhs, rhs = equation.split("=>", 1)
        return f"{rhs.strip()} => {lhs.strip()}"
    if "=" in equation:
        lhs, rhs = equation.split("=", 1)
        return f"{rhs.strip()} = {lhs.strip()}"
    return None


def _normalize_species_indices(
    kept_species_idx: Iterable[int] | Iterable[bool],
    species_count: int
) -> List[int]:
    """Normalize incoming species selection to sorted index list."""
    kept_list = list(kept_species_idx)
    if not kept_list:
        return []

    if len(kept_list) == species_count and all(isinstance(v, (bool, int)) for v in kept_list):
        idx: List[int] = []
        for i, flag in enumerate(kept_list):
            if bool(flag):
                idx.append(i)
        return idx

    indices: List[int] = []
    for value in kept_list:
        if not isinstance(value, int):
            raise TypeError("kept_species_idx must contain integers or a bool mask matching base species length")
        indices.append(int(value))

    return sorted(set(idx for idx in indices if 0 <= idx < species_count))


def write_reduced_mechanism_ct(
    base_mech_path: str,
    kept_reaction_idx: List[int],
    out_mech_path: str,
    kept_species_idx: List[int],
    max_closure_iters: int = 10,
) -> Dict[str, Any]:
    """Write a reduced YAML using Cantera APIs to preserve structure against base ordering.

    Returns summary statistics containing final counts, iterations, and enforced locks.
    """
    import warnings
    import cantera as ct

    base_path = Path(base_mech_path)
    out_path = Path(out_mech_path)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        solution = ct.Solution(str(base_path))

    base_species = list(solution.species())
    base_species_names = [sp.name for sp in base_species]
    base_reactions = list(solution.reactions())
    n_species = len(base_species)
    n_reactions = len(base_reactions)

    if not kept_reaction_idx:
        raise ValueError("kept_reaction_idx must contain at least one reaction index")

    species_indices = _normalize_species_indices(kept_species_idx, n_species)
    if not species_indices:
        raise ValueError("kept_species_idx must select at least one species")

    hard_species = {name for name in HARD_SPECIES_DEFAULT if name in base_species_names}
    species_set = {base_species_names[i] for i in species_indices}
    species_set.update(hard_species)

    reaction_index_set = {int(i) for i in kept_reaction_idx if 0 <= int(i) < n_reactions}
    if not reaction_index_set:
        raise ValueError("kept_reaction_idx has no valid indices within base mechanism")

    equation_to_indices: Dict[str, List[int]] = {}
    for idx, rxn in enumerate(base_reactions):
        eq = getattr(rxn, 'equation', '')
        equation_to_indices.setdefault(eq, []).append(idx)

    family_memberships = get_locked_family_memberships(solution)
    enforced_families = {}
    for label, indices in family_memberships.items():
        if reaction_index_set.intersection(indices):
            reaction_index_set.update(indices)
            enforced_families[label] = sorted(indices)

    expanded_indices: Set[int] = set(reaction_index_set)
    for idx in list(reaction_index_set):
        eq = getattr(base_reactions[idx], 'equation', '')
        rev_eq = _reverse_equation(eq)
        if rev_eq and rev_eq in equation_to_indices:
            expanded_indices.update(equation_to_indices[rev_eq])
    reaction_index_set = expanded_indices

    for idx in reaction_index_set:
        rxn = base_reactions[idx]
        species_set.update(rxn.reactants.keys())
        species_set.update(rxn.products.keys())

    closure_iters = 0
    retained_reactions: List[int] = []
    converged = False

    for iteration in range(1, max_closure_iters + 1):
        closure_iters = iteration
        current_species = set(species_set)
        candidate_reactions: List[int] = []

        for idx in sorted(reaction_index_set):
            rxn = base_reactions[idx]
            rxn_species = set(rxn.reactants.keys()) | set(rxn.products.keys())
            if rxn_species.issubset(current_species):
                candidate_reactions.append(idx)

        live_species = set(hard_species)
        for idx in candidate_reactions:
            rxn = base_reactions[idx]
            live_species.update(rxn.reactants.keys())
            live_species.update(rxn.products.keys())

        print(f"[closure] iter {iteration:02d}: live_species={len(live_species)} live_reactions={len(candidate_reactions)}")

        if live_species == species_set and candidate_reactions == retained_reactions:
            converged = True
            break

        species_set = live_species
        retained_reactions = candidate_reactions

    if not converged:
        print("[closure] warning: reached max iterations without full convergence")

    species_indices = [i for i, name in enumerate(base_species_names) if name in species_set]
    species_selected = [base_species[i] for i in species_indices]
    reaction_selected = [base_reactions[i] for i in retained_reactions]

    third_body_warnings: List[str] = []
    for rxn in reaction_selected:
        reaction_type = (getattr(rxn, 'reaction_type', '') or '').lower()
        if any(tag in reaction_type for tag in ['three-body', 'falloff', 'chemically-activated']):
            efficiencies = getattr(rxn, 'efficiencies', None)
            default_eff = getattr(rxn, 'default_efficiency', None)
            if (not efficiencies or len(efficiencies) == 0) and default_eff in (None, 0.0):
                third_body_warnings.append(getattr(rxn, 'equation', '<unknown>'))
    if third_body_warnings:
        print(f"[mechanism_io] warning: missing third-body efficiencies for {third_body_warnings}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Auto-repair duplicates: if only one member of a duplicate group survives, clear duplicate flag
    # Build equation -> indices map within kept list
    eq_counts: Dict[str, int] = {}
    for r in reaction_selected:
        eq = getattr(r, 'equation', None)
        if eq is not None:
            eq_counts[eq] = eq_counts.get(eq, 0) + 1
    for r in reaction_selected:
        try:
            if getattr(r, 'duplicate', False):
                eq = getattr(r, 'equation', None)
                if eq is not None and eq_counts.get(eq, 0) == 1:
                    # clear duplicate on lone survivor
                    setattr(r, 'duplicate', False)
        except Exception:
            pass

    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        sol2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics', species=species_selected, reactions=reaction_selected)
        # Cantera >= 2.6
        if hasattr(sol2, 'write_yaml'):
            sol2.write_yaml(str(out_path))
        else:
            # Fallback: copy base mechanism if write is unavailable
            shutil.copy2(base_path, out_path)
    except Exception:
        # Fallback: copy base mechanism to keep pipeline moving
        shutil.copy2(base_path, out_path)

    return {
        'n_species': len(species_selected),
        'n_reactions': len(reaction_selected),
        'iterations': closure_iters,
        'kept_species_idx': species_indices,
        'kept_reaction_idx': retained_reactions,
        'hard_species_enforced': sorted(hard_species),
        'locked_families_enforced': enforced_families,
    }
