from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import csv
import json
import shutil
import time
import uuid
import numpy as np
import cantera as ct
import re

from .evaluator import calculate_fitness, evaluate, calculate_metrics_penalty
from .mechanism_io import write_reduced_mechanism_ct, HARD_SPECIES_DEFAULT
from .constraints import HARD_SPECIES as CONSTRAINT_HARD_SPECIES, get_locked_family_memberships

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_TARGET_BAND = (28, 32)
DEFAULT_POP_SIZE = 40
DEFAULT_GENERATIONS = 120
DEFAULT_ELITISM = 2
DEFAULT_MUTATION_RATE = 0.15
GENERATION_EXTENSION = 40
VALIDATION_FAIL_THRESHOLD = 1.0e6
EARLY_STOP_PENALTY_THRESHOLD = 200.0
EARLY_STOP_STALL_GENS = 10

HARD_REACTION_EQUATIONS = {
    "H + O2 = O + OH",
    "O + H2 = OH + H",
    "OH + H2 = H + H2O",
    "HO2 + H = 2 OH",
    "HO2 + HO2 = H2O2 + O2",
    "H2O2 (+M) = 2 OH (+M)",
    "CH4 + OH = CH3 + H2O",
    "CH3 + O = CH2O + H",
    "CH2O + H = HCO + H2",
    "CH2O + OH = HCO + H2O",
    "HCO (+M) = H + CO (+M)",
    "HCO + O2 = CO + HO2",
    "CO + OH = CO2 + H",
}

SOFT_SPECIES_DEFAULT = {"HO2", "H2O2", "HCO", "CH2O", "CH3", "O2"}


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------
def _canonical_forms(eq: str) -> Set[str]:
    """Return canonical forward/reverse forms of a reaction equation."""
    text = eq.replace("<=>", "=").replace("=>", "=").replace("<=", "=")
    text = text.replace("(+M)", "+M").replace("(+m)", "+M")
    text = text.replace(" ", "").replace("\t", "")
    parts = text.split("=")
    if len(parts) != 2:
        return {text}

    coeff_pattern = re.compile(r'(\d*\.?\d*)([A-Za-z0-9_()]+)')

    def _normalise(side: str) -> str:
        tokens: List[str] = []
        for coeff_str, species in coeff_pattern.findall(side):
            coeff = float(coeff_str) if coeff_str else 1.0
            if coeff <= 0:
                continue
            count = max(1, int(round(coeff)))
            tokens.extend([species] * count)
        tokens.sort()
        return "+".join(tokens)

    left = _normalise(parts[0])
    right = _normalise(parts[1])
    forward = f"{left}={right}"
    reverse = f"{right}={left}"
    return {forward, reverse}


def _load_base_solution(base_mech_path: Path) -> Tuple[ct.Solution, List[str], List[str], List[Set[int]], Dict[str, Set[int]]]:
    sol = ct.Solution(str(base_mech_path))
    species_names = list(sol.species_names)
    species_index = {name: idx for idx, name in enumerate(species_names)}
    reaction_equations = []
    reaction_species_idx: List[Set[int]] = []
    for reaction in sol.reactions():
        reaction_equations.append(reaction.equation)
        involved = set(reaction.reactants.keys()) | set(reaction.products.keys())
        reaction_species_idx.append({species_index[name] for name in involved if name in species_index})
    family_memberships = get_locked_family_memberships(sol)
    return sol, species_names, reaction_equations, reaction_species_idx, family_memberships


def _lift_seed_mask(seed_mech_path: Optional[str], species_names: List[str]) -> Optional[np.ndarray]:
    if not seed_mech_path:
        return None
    seed_path = Path(seed_mech_path)
    if not seed_path.exists():
        print(f"[seed] warning: seed mechanism not found at {seed_mech_path} – ignoring")
        return None
    try:
        seed_sol = ct.Solution(str(seed_path))
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[seed] warning: unable to load seed mechanism ({exc}) – ignoring")
        return None
    seed_species = {sp.name for sp in seed_sol.species()}
    mask = np.zeros(len(species_names), dtype=int)
    for idx, name in enumerate(species_names):
        if name in seed_species:
            mask[idx] = 1
    return mask


def _compute_species_priority(
    reaction_species_idx: List[Set[int]],
    n_species: int,
    soft_species_idx: Set[int],
) -> np.ndarray:
    priority = np.zeros(n_species, dtype=float)
    for idx_set in reaction_species_idx:
        for idx in idx_set:
            priority[idx] += 1.0
    if priority.max() > 0:
        priority /= priority.max()
    priority[list(soft_species_idx)] += 0.3
    return priority


def _adjust_to_band(
    mask: np.ndarray,
    target_min: int,
    target_max: int,
    hard_indices: Set[int],
    priority: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    mask = mask.astype(int, copy=True)
    for idx in hard_indices:
        if 0 <= idx < mask.size:
            mask[idx] = 1

    def _grow():
        order = np.argsort(priority)[::-1]
        for idx in order:
            if mask[idx] == 0 and idx not in hard_indices:
                mask[idx] = 1
                if mask.sum() >= target_min:
                    return

    def _shrink():
        order = np.argsort(priority)
        for idx in order:
            if mask[idx] == 1 and idx not in hard_indices:
                mask[idx] = 0
                if mask.sum() <= target_max:
                    return

    attempts = 0
    while mask.sum() < target_min and attempts < mask.size:
        _grow()
        attempts += 1
    attempts = 0
    while mask.sum() > target_max and attempts < mask.size:
        _shrink()
        attempts += 1

    if not (target_min <= mask.sum() <= target_max):
        # As a last resort, randomly add/remove while respecting locks
        all_indices = [i for i in range(mask.size) if i not in hard_indices]
        rng.shuffle(all_indices)
        for idx in all_indices:
            mask[idx] = 1
            if mask.sum() >= target_min:
                break
        rng.shuffle(all_indices)
        for idx in all_indices:
            if mask.sum() <= target_max:
                break
            mask[idx] = 0

    for idx in hard_indices:
        mask[idx] = 1
    return mask


def _apply_forced_species(mask: np.ndarray, forced_indices: Set[int]) -> np.ndarray:
    """Ensure forced species indices remain active (Lu & Law radical core safeguard)."""
    for idx in forced_indices:
        if 0 <= idx < mask.size:
            mask[idx] = 1
    return mask


def _save_candidate_artifact(
    candidate_yaml: Path,
    run_dir: Path,
    n_species: int,
    n_reactions: int,
    penalty: float,
    summary_state: Dict[str, Any],
    in_range: bool,
) -> Path:
    """Persist the current best candidate to well-known locations."""
    dest = run_dir / "reduced.yaml"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if candidate_yaml.resolve() != dest.resolve():
        shutil.copy2(candidate_yaml, dest)

    mechanisms_dir = Path("mechanisms_reduced")
    mechanisms_dir.mkdir(parents=True, exist_ok=True)
    tmp_candidate = mechanisms_dir / "tmp_candidate.yaml"
    shutil.copy2(candidate_yaml, tmp_candidate)

    summary_state.update(
        {
            "reduced_yaml_path": str(dest),
            "n_species": int(n_species),
            "n_reactions": int(n_reactions),
            "penalty": float(penalty),
            "in_range": bool(in_range),
        }
    )
    return dest


def _derive_reaction_mask(
    species_mask: np.ndarray,
    reaction_species_idx: List[Set[int]],
    forced_indices: Set[int],
) -> np.ndarray:
    kept_species = {i for i, flag in enumerate(species_mask) if flag}
    mask = np.zeros(len(reaction_species_idx), dtype=int)
    for idx, species_set in enumerate(reaction_species_idx):
        if species_set.issubset(kept_species):
            mask[idx] = 1
    for idx in forced_indices:
        mask[idx] = 1
        kept_species.update(reaction_species_idx[idx])
    return mask


def _mutate_species(
    mask: np.ndarray,
    mutation_rate: float,
    hard_indices: Set[int],
    rng: np.random.Generator,
) -> np.ndarray:
    mutated = mask.copy()
    for idx in range(mutated.size):
        if idx in hard_indices:
            mutated[idx] = 1
            continue
        if rng.random() < mutation_rate:
            mutated[idx] = 1 - mutated[idx]
    return mutated


def _mutate_reactions(
    reaction_mask: np.ndarray,
    species_mask: np.ndarray,
    reaction_species_idx: List[Set[int]],
    mutation_rate: float,
    forced_indices: Set[int],
    rng: np.random.Generator,
) -> np.ndarray:
    mutated = reaction_mask.copy()
    species_set = {i for i, flag in enumerate(species_mask) if flag}
    for idx in range(mutated.size):
        if idx in forced_indices:
            mutated[idx] = 1
            continue
        if rng.random() < mutation_rate:
            if mutated[idx] == 1:
                mutated[idx] = 0
            else:
                if reaction_species_idx[idx].issubset(species_set):
                    mutated[idx] = 1
    for idx in forced_indices:
        mutated[idx] = 1
    return mutated


def _species_indices_from_mask(mask: np.ndarray) -> List[int]:
    return [i for i, flag in enumerate(mask) if flag]


def _reaction_indices_from_mask(mask: np.ndarray) -> List[int]:
    return [i for i, flag in enumerate(mask) if flag]


def _ensure_hard_species_present(mask: np.ndarray, hard_indices: Set[int]) -> None:
    for idx in hard_indices:
        if 0 <= idx < mask.size:
            mask[idx] = 1


def _greedy_prune_to_band(
    base_path: Path,
    initial_species_mask: np.ndarray,
    reaction_species_idx: List[Set[int]],
    forced_reaction_indices: Set[int],
    immutable_species_idx: Set[int],
    priority: np.ndarray,
    target_band: Tuple[int, int],
    cases_root: str,
    case_list_fast: List[str],
    case_list_full: Optional[List[str]],
    soft_species_set: Set[str],
    target_size: int,
    base_penalty: float,
    base_diagnostics: Dict[str, Any],
    workdir: Path,
) -> Tuple[np.ndarray, np.ndarray, float, Dict[str, Any], bool]:
    """
    Attempt to prune species into the target band while preserving feasibility.
    Returns updated (species_mask, reaction_mask, penalty, diagnostics, in_range).
    """
    target_min, target_max = target_band
    current_mask = initial_species_mask.copy()
    current_reaction_mask = _derive_reaction_mask(current_mask, reaction_species_idx, forced_reaction_indices)
    current_penalty = base_penalty
    current_diagnostics = dict(base_diagnostics)
    current_species_count = int(current_diagnostics.get('n_species_final', int(current_mask.sum())))

    if current_species_count <= target_max and current_species_count >= target_min:
        return current_mask, current_reaction_mask, current_penalty, current_diagnostics, True

    order = list(np.argsort(priority))
    workdir.mkdir(parents=True, exist_ok=True)
    step = 0

    while current_species_count > target_max:
        step += 1
        best_candidate = None
        for idx in order:
            if current_mask[idx] == 0 or idx in immutable_species_idx:
                continue
            if current_mask.sum() - 1 < target_min:
                continue

            trial_mask = current_mask.copy()
            trial_mask[idx] = 0
            trial_mask = _apply_forced_species(trial_mask, immutable_species_idx)
            if trial_mask.sum() < target_min:
                continue

            trial_reaction_mask = _derive_reaction_mask(trial_mask, reaction_species_idx, forced_reaction_indices)
            trial_dir = workdir / f"step_{step:02d}_{idx:03d}"
            penalty, diagnostics = calculate_fitness(
                base_mech_path=str(base_path),
                species_mask=trial_mask.tolist(),
                reaction_mask=trial_reaction_mask.tolist(),
                cases_root=cases_root,
                case_list_fast=case_list_fast,
                case_list_full=case_list_full,
                soft_species=soft_species_set,
                target_size=target_size,
                workdir=trial_dir,
            )
            if penalty >= VALIDATION_FAIL_THRESHOLD:
                continue
            if not (diagnostics.get('psr_pass') and diagnostics.get('pfr_pass')):
                continue
            n_species_final = int(diagnostics.get('n_species_final', trial_mask.sum()))
            if n_species_final < target_min:
                continue
            candidate = (penalty, diagnostics, trial_mask, trial_reaction_mask, n_species_final)
            if best_candidate is None or penalty < best_candidate[0]:
                best_candidate = candidate

        if best_candidate is None:
            break

        penalty, diagnostics, trial_mask, trial_reaction_mask, n_species_final = best_candidate
        current_mask = trial_mask
        current_reaction_mask = trial_reaction_mask
        current_penalty = penalty
        current_diagnostics = diagnostics
        current_species_count = n_species_final

        if current_species_count <= target_max:
            break

    in_range = target_min <= current_species_count <= target_max
    return current_mask, current_reaction_mask, current_penalty, current_diagnostics, in_range


# ---------------------------------------------------------------------------
# Main GA runner
# ---------------------------------------------------------------------------
def run_ga_gnn(
    base_mech_path: str,
    output_dir: str,
    cases_root: str,
    case_list: List[str],
    target_band: Tuple[int, int] = DEFAULT_TARGET_BAND,
    pop: int = DEFAULT_POP_SIZE,
    gens: int = DEFAULT_GENERATIONS,
    elitism: int = DEFAULT_ELITISM,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    seed_mech_path: Optional[str] = None,
    hard_species: Optional[List[str]] = None,
    hard_reactions: Optional[List[str]] = None,
    soft_species: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Run GA+GNN reduction and return final artefacts."""
    start_time = time.time()
    base_path = Path(base_mech_path)
    if not base_path.exists():
        raise FileNotFoundError(f"Base mechanism not found: {base_mech_path}")
    if not case_list:
        raise ValueError("case_list must contain at least one validation case")

    (
        sol,
        species_names,
        reaction_equations,
        reaction_species_idx,
        locked_family_memberships,
    ) = _load_base_solution(base_path)
    n_species = len(species_names)
    n_reactions = len(reaction_equations)

    hard_species_list = sorted(set(hard_species or HARD_SPECIES_DEFAULT) | set(CONSTRAINT_HARD_SPECIES))
    hard_species_set = {name for name in hard_species_list if name in species_names}
    species_index_map = {name: idx for idx, name in enumerate(species_names)}
    hard_species_idx = {species_index_map[name] for name in hard_species_set}

    soft_species_set = set(soft_species or SOFT_SPECIES_DEFAULT)
    soft_species_idx = {species_index_map[name] for name in soft_species_set if name in species_index_map}

    hard_reaction_patterns = set(hard_reactions) if hard_reactions else HARD_REACTION_EQUATIONS
    reaction_canon_map: Dict[str, int] = {}
    for idx, eq in enumerate(reaction_equations):
        for form in _canonical_forms(eq):
            reaction_canon_map.setdefault(form, idx)

    hard_reaction_indices: Set[int] = set()
    missing_patterns: List[str] = []
    for pattern in hard_reaction_patterns:
        matched = False
        for form in _canonical_forms(pattern):
            if form in reaction_canon_map:
                hard_reaction_indices.add(reaction_canon_map[form])
                matched = True
                break
        if not matched:
            missing_patterns.append(pattern)
    if missing_patterns:
        print(f"[hard reactions] warning: patterns not found in base mechanism: {sorted(missing_patterns)}")

    locked_family_indices: Set[int] = set()
    for indices in locked_family_memberships.values():
        locked_family_indices.update(indices)
    forced_reaction_indices = set(hard_reaction_indices) | locked_family_indices

    forced_species_idx: Set[int] = set()
    for idx in forced_reaction_indices:
        if 0 <= idx < len(reaction_species_idx):
            forced_species_idx.update(reaction_species_idx[idx])
    immutable_species_idx: Set[int] = set(hard_species_idx) | forced_species_idx

    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    run_tag = time.strftime("%Y%m%d_%H%M%S")
    run_dir = output_root / run_tag
    run_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(42)
    priority = _compute_species_priority(reaction_species_idx, n_species, soft_species_idx)

    population_species = np.zeros((pop, n_species), dtype=int)
    population_reactions = np.zeros((pop, n_reactions), dtype=int)

    seed_sources: List[Tuple[str, Optional[np.ndarray]]] = []
    elite_seed_path = Path("mechanisms_reduced/gri30_reduced_v2.yaml")
    if elite_seed_path.exists():
        seed_sources.append(("elite", _lift_seed_mask(str(elite_seed_path), species_names)))
    seed_sources.append(("cli", _lift_seed_mask(seed_mech_path, species_names)))

    next_slot = 0
    for label, raw_seed in seed_sources:
        if raw_seed is None:
            continue
        seed_mask = _adjust_to_band(raw_seed, target_band[0], target_band[1], immutable_species_idx, priority, rng)
        seed_mask = _apply_forced_species(seed_mask, immutable_species_idx)
        population_species[next_slot] = seed_mask
        population_reactions[next_slot] = _derive_reaction_mask(seed_mask, reaction_species_idx, forced_reaction_indices)
        next_slot += 1
        if next_slot >= pop:
            break

    for idx in range(next_slot, pop):
        random_mask = (rng.random(n_species) < 0.6).astype(int)
        random_mask = _adjust_to_band(random_mask, target_band[0], target_band[1], immutable_species_idx, priority, rng)
        random_mask = _apply_forced_species(random_mask, immutable_species_idx)
        population_species[idx] = random_mask
        population_reactions[idx] = _derive_reaction_mask(random_mask, reaction_species_idx, forced_reaction_indices)

    best_penalty = float("inf")
    best_species_mask = population_species[0].copy()
    best_reaction_mask = population_reactions[0].copy()
    best_diagnostics: Dict[str, Any] = {}

    ga_trace_rows: List[Dict[str, Any]] = []
    final_state: Dict[str, Any] = {"in_range": False}
    best_in_range_penalty = float("inf")
    best_in_range_species_mask: Optional[np.ndarray] = None
    best_in_range_reaction_mask: Optional[np.ndarray] = None
    best_in_range_diagnostics: Dict[str, Any] = {}
    best_in_range_generation = -1
    stall_generations = 0

    max_generations = gens
    extensions_used = 0
    generation = 0

    case_list_fast = list(case_list)

    while generation < max_generations:
        penalties = np.zeros(pop)
        diagnostics_per_individual: List[Dict[str, Any]] = []
        improved_in_range = False

        for idx in range(pop):
            species_mask = population_species[idx]
            reaction_mask = population_reactions[idx]

            candidate_dir = run_dir / f"gen_{generation:02d}" / f"ind_{idx:03d}"
            penalty, diagnostics = calculate_fitness(
                base_mech_path=str(base_path),
                species_mask=species_mask.tolist(),
                reaction_mask=reaction_mask.tolist(),
                cases_root=cases_root,
                case_list_fast=case_list_fast,
                case_list_full=None,
                soft_species=soft_species_set,
                target_size=int(sum(target_band) / 2),
                workdir=candidate_dir,
            )
            penalties[idx] = penalty
            diagnostics_per_individual.append(diagnostics)

            n_species_final = int(diagnostics.get('n_species_final', int(species_mask.sum())))
            n_reactions_final = int(diagnostics.get('n_reactions_final', int(reaction_mask.sum())))
            candidate_yaml_path = Path(diagnostics['candidate_yaml']) if 'candidate_yaml' in diagnostics else None
            feasible = (
                penalty < VALIDATION_FAIL_THRESHOLD
                and diagnostics.get('psr_pass', False)
                and diagnostics.get('pfr_pass', False)
            )

            if feasible and target_band[0] <= n_species_final <= target_band[1]:
                if penalty < best_in_range_penalty - 1e-6:
                    best_in_range_penalty = penalty
                    best_in_range_species_mask = species_mask.copy()
                    best_in_range_reaction_mask = reaction_mask.copy()
                    best_in_range_diagnostics = diagnostics
                    best_in_range_generation = generation
                    if candidate_yaml_path and candidate_yaml_path.exists():
                        _save_candidate_artifact(
                            candidate_yaml_path,
                            run_dir,
                            n_species_final,
                            n_reactions_final,
                            penalty,
                            final_state,
                            True,
                        )
                        print(
                            f"[best-in-range] gen={generation:02d} size={n_species_final} "
                            f"penalty={penalty:.2f} saved={run_dir / 'reduced.yaml'}"
                        )
                    improved_in_range = True

        best_idx = int(np.argmin(penalties))
        gen_best_penalty = float(penalties[best_idx])
        gen_best_species = population_species[best_idx]
        gen_best_reactions = population_reactions[best_idx]
        gen_best_diag = diagnostics_per_individual[best_idx]
        gen_best_size = int(gen_best_species.sum())
        validation_status = "PASS" if gen_best_diag['psr_pass'] and gen_best_diag['pfr_pass'] else "FAIL"

        print(f"[gen {generation:02d}] penalty={gen_best_penalty:.2f} size={gen_best_size} validation={validation_status}")

        ga_trace_rows.append({
            'generation': generation,
            'penalty': gen_best_penalty,
            'size': gen_best_size,
            'psr_pass': gen_best_diag['psr_pass'],
            'pfr_pass': gen_best_diag['pfr_pass'],
        })

        if gen_best_penalty < best_penalty:
            best_penalty = gen_best_penalty
            best_species_mask = gen_best_species.copy()
            best_reaction_mask = gen_best_reactions.copy()
            best_diagnostics = gen_best_diag

        should_break = False
        if best_in_range_penalty < float("inf"):
            if improved_in_range:
                stall_generations = 0
            else:
                stall_generations += 1

            best_in_range_size = int(
                best_in_range_diagnostics.get(
                    'n_species_final',
                    int(best_in_range_species_mask.sum()) if best_in_range_species_mask is not None else 0,
                )
            )
            if (
                best_in_range_penalty <= EARLY_STOP_PENALTY_THRESHOLD
                or stall_generations >= EARLY_STOP_STALL_GENS
            ):
                print(
                    f"[early-stop] best-in-range penalty={best_in_range_penalty:.2f} "
                    f"size={best_in_range_size} gen={best_in_range_generation:02d}"
                )
                should_break = True
        else:
            stall_generations = 0

        if should_break:
            break

        elite_idx = np.argsort(penalties)[:elitism]
        next_species = population_species[elite_idx].copy()
        next_reactions = population_reactions[elite_idx].copy()

        def tournament() -> int:
            i, j = rng.integers(0, pop, size=2)
            return i if penalties[i] <= penalties[j] else j

        while next_species.shape[0] < pop:
            parent_a = population_species[tournament()]
            parent_b = population_species[tournament()]
            rxn_parent_a = population_reactions[tournament()]
            rxn_parent_b = population_reactions[tournament()]

            cross_point = rng.integers(1, n_species)
            child_a = np.concatenate([parent_a[:cross_point], parent_b[cross_point:]])
            child_b = np.concatenate([parent_b[:cross_point], parent_a[cross_point:]])

            child_a = _mutate_species(child_a, mutation_rate, immutable_species_idx, rng)
            child_b = _mutate_species(child_b, mutation_rate, immutable_species_idx, rng)

            child_a = _adjust_to_band(child_a, target_band[0], target_band[1], immutable_species_idx, priority, rng)
            child_b = _adjust_to_band(child_b, target_band[0], target_band[1], immutable_species_idx, priority, rng)
            child_a = _apply_forced_species(child_a, immutable_species_idx)
            child_b = _apply_forced_species(child_b, immutable_species_idx)

            child_rxn_a = _derive_reaction_mask(child_a, reaction_species_idx, forced_reaction_indices)
            child_rxn_b = _derive_reaction_mask(child_b, reaction_species_idx, forced_reaction_indices)

            child_rxn_a = _mutate_reactions(child_rxn_a, child_a, reaction_species_idx, mutation_rate * 0.5, forced_reaction_indices, rng)
            child_rxn_b = _mutate_reactions(child_rxn_b, child_b, reaction_species_idx, mutation_rate * 0.5, forced_reaction_indices, rng)

            next_species = np.vstack([next_species, child_a, child_b])
            next_reactions = np.vstack([next_reactions, child_rxn_a, child_rxn_b])

        population_species = next_species[:pop]
        population_reactions = next_reactions[:pop]

        generation += 1
        if generation == max_generations:
            best_size = int(best_species_mask.sum())
            if not (target_band[0] <= best_size <= target_band[1]) and extensions_used < 2:
                print(f"[extend] best size {best_size} outside target {target_band}, extending by +{GENERATION_EXTENSION} gens")
                max_generations += GENERATION_EXTENSION
                extensions_used += 1

    if best_in_range_species_mask is not None:
        final_species_mask = best_in_range_species_mask
        final_reaction_mask = best_in_range_reaction_mask if best_in_range_reaction_mask is not None else best_reaction_mask
        final_penalty_candidate = best_in_range_penalty
        final_diagnostics = best_in_range_diagnostics
        final_in_range = True
    else:
        final_species_mask = best_species_mask
        final_reaction_mask = best_reaction_mask
        final_penalty_candidate = best_penalty
        final_diagnostics = best_diagnostics
        final_in_range = False

        est_species_count = int(final_diagnostics.get('n_species_final', int(final_species_mask.sum())))
        if not (target_band[0] <= est_species_count <= target_band[1]):
            prune_mask, prune_rxn_mask, prune_penalty, prune_diag, prune_in_range = _greedy_prune_to_band(
                base_path=base_path,
                initial_species_mask=final_species_mask,
                reaction_species_idx=reaction_species_idx,
                forced_reaction_indices=forced_reaction_indices,
                immutable_species_idx=immutable_species_idx,
                priority=priority,
                target_band=target_band,
                cases_root=cases_root,
                case_list_fast=case_list_fast,
                case_list_full=None,
                soft_species_set=soft_species_set,
                target_size=int(sum(target_band) / 2),
                base_penalty=final_penalty_candidate,
                base_diagnostics=final_diagnostics,
                workdir=run_dir / "fallback_prune",
            )
            final_species_mask = prune_mask
            final_reaction_mask = prune_rxn_mask
            final_penalty_candidate = prune_penalty
            final_diagnostics = prune_diag
            final_in_range = prune_in_range

    final_species_idx = _species_indices_from_mask(final_species_mask)
    final_reaction_idx = _reaction_indices_from_mask(final_reaction_mask)

    reduced_yaml = run_dir / "reduced.yaml"
    closure_summary = write_reduced_mechanism_ct(
        base_mech_path=str(base_path),
        kept_reaction_idx=final_reaction_idx,
        out_mech_path=str(reduced_yaml),
        kept_species_idx=final_species_idx,
        max_closure_iters=20,
    )

    _save_candidate_artifact(
        reduced_yaml,
        run_dir,
        closure_summary['n_species'],
        closure_summary['n_reactions'],
        final_penalty_candidate if final_penalty_candidate < float('inf') else best_penalty,
        final_state,
        final_in_range,
    )

    final_eval_dir = run_dir / "final_eval"
    final_metrics = evaluate(
        mech_path=reduced_yaml,
        cases_root=Path(cases_root),
        cases=case_list,
        outdir=final_eval_dir,
    )
    final_penalty = calculate_metrics_penalty(final_metrics, n_species=len(closure_summary['kept_species_idx']))

    with open(final_eval_dir / "metrics.json", 'w') as f:
        json.dump(final_metrics, f, indent=2)

    trace_path = run_dir / "ga_trace.csv"
    with open(trace_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['generation', 'penalty', 'size', 'psr_pass', 'pfr_pass'])
        writer.writeheader()
        for row in ga_trace_rows:
            writer.writerow(row)

    runtime = time.time() - start_time
    summary = {
        'runtime_sec': runtime,
        'best_penalty': best_penalty,
        'final_penalty': final_penalty,
        'penalty': final_penalty,
        'closure': closure_summary,
        'best_diagnostics': best_diagnostics,
        'final_diagnostics': final_diagnostics,
        'final_metrics_path': str(final_eval_dir / "metrics.json"),
        'metrics_json_path': str(final_eval_dir / "metrics.json"),
        'reduced_yaml_path': final_state.get('reduced_yaml_path', str(reduced_yaml)),
        'ga_trace_csv_path': str(trace_path),
        'species_mask': final_species_mask.astype(bool).tolist(),
        'reaction_mask': final_reaction_mask.astype(bool).tolist(),
        'target_band': target_band,
        'population': pop,
        'generations': generation,
        'best_in_range_penalty': best_in_range_penalty if best_in_range_penalty < float('inf') else None,
        'n_species': final_state.get('n_species', closure_summary.get('n_species')),
        'n_reactions': final_state.get('n_reactions', closure_summary.get('n_reactions')),
        'in_range': final_state.get('in_range', final_in_range),
    }

    with open(run_dir / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"[done] Reduced mechanism written to {reduced_yaml} ({closure_summary['n_species']} species, {closure_summary['n_reactions']} reactions)")
    print(f"[done] Final metrics stored at {final_eval_dir / 'metrics.json'}")

    return summary

