from __future__ import annotations

"""
Reduction hard-constraints inspired by Lu & Law (Combust. Flame 2005/2006)
and Pepiot-Desjardins & Pitsch (Combust. Flame 2008).

The radical and syngas backbone reactions listed here are treated as
non-negotiable when building reduced methane partial-oxidation mechanisms.
"""

from typing import Dict, Iterable, Tuple, Set

HARD_SPECIES = [
    "H2",
    "H",
    "O",
    "O2",
    "OH",
    "HO2",
    "H2O2",
    "H2O",
    "CO",
    "CO2",
    "HCO",
    "CH2O",
    "CH3",
    "CH4",
    "N2",
    "AR",
]
MANDATORY_SPECIES = HARD_SPECIES.copy()

LockTuple = Tuple[frozenset[str], frozenset[str], str]
HARD_FAMILIES: list[dict[str, LockTuple]] = []


def lock_family(
    forward: Tuple[Iterable[str], Iterable[str], str],
    reverse: Tuple[Iterable[str], Iterable[str], str],
) -> None:
    """
    Register a forward/reverse reaction pattern that must be kept together.

    Parameters
    ----------
    forward / reverse:
        Tuples of (reactant species, product species, label). Stoichiometry is
        ignored; only species membership matters. The label is stored so we
        can report which Lu & Law / Pepiot-Desjardins family triggered a lock.
    """

    f_lhs, f_rhs, f_label = forward
    r_lhs, r_rhs, r_label = reverse
    HARD_FAMILIES.append(
        {
            "forward": (frozenset(f_lhs), frozenset(f_rhs), f_label),
            "reverse": (frozenset(r_lhs), frozenset(r_rhs), r_label),
            "label": f_label,
        }
    )


# (a) Methane activation CH4 + X -> CH3 + … with X ∈ {O, OH, H, HO2}
lock_family(
    ({"CH4", "O"}, {"CH3"}, "CH4_activation_O"),
    ({"CH3"}, {"CH4", "O"}, "CH4_activation_O_rev"),
)
lock_family(
    ({"CH4", "OH"}, {"CH3"}, "CH4_activation_OH"),
    ({"CH3"}, {"CH4", "OH"}, "CH4_activation_OH_rev"),
)
lock_family(
    ({"CH4", "H"}, {"CH3"}, "CH4_activation_H"),
    ({"CH3"}, {"CH4", "H"}, "CH4_activation_H_rev"),
)
lock_family(
    ({"CH4", "HO2"}, {"CH3"}, "CH4_activation_HO2"),
    ({"CH3"}, {"CH4", "HO2"}, "CH4_activation_HO2_rev"),
)

# (b) CO-formation chain
lock_family(
    ({"CH3"}, {"CH2O"}, "CO_path_CH3_to_CH2O"),
    ({"CH2O"}, {"CH3"}, "CO_path_CH2O_to_CH3"),
)
lock_family(
    ({"CH2O"}, {"HCO"}, "CO_path_CH2O_to_HCO"),
    ({"HCO"}, {"CH2O"}, "CO_path_HCO_to_CH2O"),
)
lock_family(
    ({"HCO"}, {"CO"}, "CO_path_HCO_to_CO"),
    ({"CO"}, {"HCO"}, "CO_path_CO_to_HCO"),
)

# (c) Water-gas shift loop
lock_family(
    ({"CO", "H2O"}, {"CO2", "H2"}, "WGS_forward"),
    ({"CO2", "H2"}, {"CO", "H2O"}, "WGS_reverse"),
)

# (d) HOx chain / radical pool maintenance
lock_family(
    ({"H", "O2"}, {"O", "OH"}, "HOx_chain_pair"),
    ({"O", "OH"}, {"H", "O2"}, "HOx_chain_pair_rev"),
)
lock_family(
    ({"OH", "CO"}, {"CO2", "H"}, "HOx_CO_oxidation"),
    ({"CO2", "H"}, {"OH", "CO"}, "HOx_CO_oxidation_rev"),
)
lock_family(
    ({"H", "OH"}, {"H2O"}, "HOx_recombination"),
    ({"H2O"}, {"H", "OH"}, "HOx_recombination_rev"),
)


def _third_body_species() -> Set[str]:
    return {"M", "AR", "N2", "H2O", "CO2", "H2"}


def normalize_reaction_sets(rxn):
    """Return reactant/product species sets minus common third bodies."""

    try:
        lhs_set = set(rxn.reactants.keys())
        rhs_set = set(rxn.products.keys())
        third_bodies = _third_body_species()
        lhs_set -= third_bodies
        rhs_set -= third_bodies
        return lhs_set, rhs_set
    except Exception:
        return set(), set()


def matches(lhs_req: Set[str], rhs_req: Set[str], lhs_set: Set[str], rhs_set: Set[str]) -> bool:
    """True when required species are subsets of the actual sides."""

    return lhs_req.issubset(lhs_set) and rhs_req.issubset(rhs_set)


def is_critical_ch4_path(rxn) -> bool:
    """
    Backward-compatible methane oxidation guard used by legacy tooling.

    The broader HARD_FAMILIES logic supersedes this, but we keep the helper so
    older reports still flag the classical four-path CH2O/HCO loops.
    """

    try:
        lhs_set, rhs_set = normalize_reaction_sets(rxn)
        if matches({"CH2O", "O2"}, {"HCO", "HO2"}, lhs_set, rhs_set):
            return True
        if matches({"CH2O", "OH"}, {"H2O", "HCO"}, lhs_set, rhs_set):
            return True
        if matches({"HCO", "OH"}, {"CO", "H2O"}, lhs_set, rhs_set):
            return True
        if matches({"HCO", "O2"}, {"CO", "HO2"}, lhs_set, rhs_set):
            return True
        return False
    except Exception:
        return False


def count_ch4_pathway_families(solution):
    """Count how many reactions map to each legacy CH4 oxidation family."""

    family_counts = {
        "CH2O + O2 -> HCO + HO2": 0,
        "CH2O + OH -> H2O + HCO": 0,
        "HCO + OH -> CO + H2O": 0,
        "HCO + O2 -> CO + HO2": 0,
    }

    for rxn in solution.reactions():
        try:
            lhs_set, rhs_set = normalize_reaction_sets(rxn)
            if matches({"CH2O", "O2"}, {"HCO", "HO2"}, lhs_set, rhs_set):
                family_counts["CH2O + O2 -> HCO + HO2"] += 1
            if matches({"CH2O", "OH"}, {"H2O", "HCO"}, lhs_set, rhs_set):
                family_counts["CH2O + OH -> H2O + HCO"] += 1
            if matches({"HCO", "OH"}, {"CO", "H2O"}, lhs_set, rhs_set):
                family_counts["HCO + OH -> CO + H2O"] += 1
            if matches({"HCO", "O2"}, {"CO", "HO2"}, lhs_set, rhs_set):
                family_counts["HCO + O2 -> CO + HO2"] += 1
        except Exception:
            continue

    return family_counts


def check_mandatory_species(species_list: list[str]) -> bool:
    """Ensure all hard-protected radicals and bath gases remain present."""

    missing = [s for s in MANDATORY_SPECIES if s not in species_list]
    if missing:
        print(f"Constraint fail: missing species {missing}")
        return False
    return True


def get_locked_family_memberships(solution) -> Dict[str, Set[int]]:
    memberships: Dict[str, Set[int]] = {fam["label"]: set() for fam in HARD_FAMILIES}
    for idx, rxn in enumerate(solution.reactions()):
        lhs_set, rhs_set = normalize_reaction_sets(rxn)
        for fam in HARD_FAMILIES:
            f_lhs, f_rhs, _ = fam["forward"]
            r_lhs, r_rhs, _ = fam["reverse"]
            if matches(f_lhs, f_rhs, lhs_set, rhs_set) or matches(r_lhs, r_rhs, lhs_set, rhs_set):
                memberships[fam["label"]].add(idx)
    return memberships


def protected_rxn_indices(solution) -> set[int]:
    """
    Return the indices of reactions that must be preserved.

    Includes: pressure-dependent forms, duplicates, methane oxidation motifs,
    and the Lu & Law / Pepiot hard families registered via ``lock_family``.
    """

    protected = set()
    for i, r in enumerate(solution.reactions()):
        rt = getattr(r, "reaction_type", "")
        if any(key in str(rt).lower() for key in ["falloff", "three", "pressure", "plog", "chebyshev"]):
            protected.add(i)
        if getattr(r, "duplicate", False):
            protected.add(i)
        try:
            if set(r.reactants.keys()) and set(r.reactants.keys()).issubset({"N2"}):
                protected.add(i)
        except Exception:
            pass
        if is_critical_ch4_path(r):
            protected.add(i)
        eq = getattr(r, "equation", "") or ""
        if any(key in eq for key in ["HO2", "H2O2"]):
            protected.add(i)
        if "CH4" in eq and any(x in eq for x in [" O ", " OH", " O(", "O2"]):
            protected.add(i)
        if any(x in eq for x in ["CH3 + O", "CH3 + OH", "CH3 + O2"]):
            protected.add(i)
        if any(x in eq for x in ["CH2O", "CH3O", "HCO"]):
            protected.add(i)
        if "CO + OH" in eq or "CO2 + H" in eq:
            protected.add(i)

    for indices in get_locked_family_memberships(solution).values():
        protected.update(indices)

    eq_to_indices: Dict[str, list[int]] = {}
    rxns = solution.reactions()
    for i, r in enumerate(rxns):
        eq = getattr(r, "equation", None)
        if eq is None:
            try:
                lhs = " + ".join(r.reactants.keys())
                rhs = " + ".join(r.products.keys())
                eq = f"{lhs} => {rhs}"
            except Exception:
                eq = str(i)
        eq_to_indices.setdefault(eq, []).append(i)

    for eq, idxs in eq_to_indices.items():
        if len(idxs) > 1:
            protected.update(idxs)
        else:
            i = idxs[0]
            try:
                if getattr(rxns[i], "duplicate", False):
                    protected.update(idxs)
            except Exception:
                pass

    return protected


def filter_mask(solution, mask):
    """Apply hard constraints to a candidate reaction mask."""

    import numpy as np

    mask_arr = np.array(mask, dtype=int).copy()
    forced = protected_rxn_indices(solution)
    for idx in forced:
        if 0 <= idx < len(mask_arr):
            mask_arr[idx] = 1

    family_membership = get_locked_family_memberships(solution)
    for indices in family_membership.values():
        if indices and any(mask_arr[j] == 1 for j in indices):
            for j in indices:
                if 0 <= j < len(mask_arr):
                    mask_arr[j] = 1

    eq_map: Dict[str, list[int]] = {}
    for i, r in enumerate(solution.reactions()):
        eq = getattr(r, "equation", str(i))
        eq_map.setdefault(eq, []).append(i)
    for eq, idxs in eq_map.items():
        if len(idxs) > 1 and any(mask_arr[j] == 1 for j in idxs):
            for j in idxs:
                mask_arr[j] = 1

    return mask_arr
