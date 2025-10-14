from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Any, Iterable, Optional, Tuple, Set, NamedTuple, Union
import json
import uuid
import numpy as np
import pandas as pd

from hp_pox_physics.config_cases import CaseConstants
from hp_pox_physics.hp_pox_model import HPPOXModel
from .mechanism_io import write_reduced_mechanism_ct

PSR_TEMP_TOL_K = 50.0
PSR_SPECIES_TOL = 0.15
PFR_SPECIES_TOL = 0.10
CH4_SLIP_TOL = 0.10  # fraction
SIZE_BIAS_WEIGHT = 35.0
SOFT_SPECIES_WEIGHT = 5.0
VALIDATION_FAIL_PENALTY = 1.0e6


class Penalty(NamedTuple):
    """Structured feasibility penalty packet for GA consumption."""
    value: float
    reason: str
    summary: str
    meta: Optional[Dict[str, Any]] = None


def _format_eval_summary(
    case: str,
    psr_converged: bool,
    psr_temp_k: float,
    x_o2: float,
    x_ch4: float,
    x_co: float,
    x_h2: float,
    ratio: float,
    penalty_value: float,
    reason: str,
) -> str:
    psr_flag = "PASS" if psr_converged else "FAIL"
    temp_txt = "nan"
    if np.isfinite(psr_temp_k):
        temp_txt = f"{psr_temp_k:.0f}"
    return (
        f"[eval] case={case} PSR={psr_flag} T={temp_txt}K  "
        f"PFR_out: X_O2={x_o2:.2f} X_CH4={x_ch4:.2f} X_CO={x_co:.2f} X_H2={x_h2:.2f}  "
        f"H2/CO={ratio:.2f} penalty={penalty_value:.3g} [{reason}]"
    )


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a); b = np.asarray(b)
    if len(a) == 0 or len(a) != len(b):
        return float('nan')
    return float(np.sqrt(np.mean((a - b) ** 2)))


def evaluate(mech_path: Path, cases_root: Path, cases: List[str], outdir: Path) -> Union[Dict[str, Any], Penalty]:
    """
    Run HP-POX validation cases and return metrics or a hard feasibility penalty.

    Enforces Lu & Law / Pepiot-desjardins feasibility rules:
      * PSR must converge to T >= 1200 K.
      * Outlet must contain syngas (H2+CO) and avoid simultaneous CH4/O2 slip.
      * H2/CO stays near 2.0 with a soft penalty otherwise.
    """
    mech_path = Path(mech_path)
    outdir.mkdir(parents=True, exist_ok=True)
    model = HPPOXModel(mech_path, transport_model='mixture-averaged')

    per_case: Dict[str, Any] = {}
    h2co_vals: List[float] = []
    ch4_vals: List[float] = []
    t_vals: List[float] = []
    feasibility_penalty = 0.0

    inlet_cache = model.load_inlet_streams()

    for case in cases:
        base_csv = cases_root / case / 'axial_profiles.csv'
        if not base_csv.exists():
            continue
        df_base = pd.read_csv(base_csv)

        inlet_config = inlet_cache
        results = model.run_case(case, inlet_config, outdir, chemistry=True)
        model.save_results(results, outdir / case)

        df_red = pd.read_csv(outdir / case / 'axial_profiles.csv')

        dry_comp = results.get('dry_composition', {}) or {}
        outlet_state = results.get('outlet_state')
        t_out_k = float(getattr(outlet_state, "temperature", float("nan")))

        psr_diag = results.get('psr_diagnostics', {}) or {}
        psr_converged = bool(psr_diag.get('converged', True))
        psr_temp = float(psr_diag.get('T_out', getattr(results.get('psr_outlet'), "temperature", float("nan"))))

        x_h2 = float(dry_comp.get('H2', 0.0))
        x_co = float(dry_comp.get('CO', 0.0))
        x_o2 = float(dry_comp.get('O2', 0.0))
        x_ch4 = float(dry_comp.get('CH4', 0.0))

        ratio = x_h2 / max(x_co, 1e-12)

        if (not psr_converged) or (np.isfinite(psr_temp) and psr_temp < 1200.0):
            summary = _format_eval_summary(case, psr_converged, psr_temp, x_o2, x_ch4, x_co, x_h2, ratio, VALIDATION_FAIL_PENALTY, "COLD_PSR")
            print(summary)
            return Penalty(VALIDATION_FAIL_PENALTY, "COLD_PSR", summary, {
                "case": case,
                "psr_converged": psr_converged,
                "psr_temperature_K": psr_temp,
                "dry_composition": dry_comp,
            })

        if x_co + x_h2 < 1e-3:
            summary = _format_eval_summary(case, psr_converged, psr_temp, x_o2, x_ch4, x_co, x_h2, ratio, VALIDATION_FAIL_PENALTY, "NO_SYNGAS")
            print(summary)
            return Penalty(VALIDATION_FAIL_PENALTY, "NO_SYNGAS", summary, {
                "case": case,
                "dry_composition": dry_comp,
            })

        if x_o2 > 0.05 and x_ch4 > 0.05:
            summary = _format_eval_summary(case, psr_converged, psr_temp, x_o2, x_ch4, x_co, x_h2, ratio, VALIDATION_FAIL_PENALTY, "UNBURNT_MIX")
            print(summary)
            return Penalty(VALIDATION_FAIL_PENALTY, "UNBURNT_MIX", summary, {
                "case": case,
                "dry_composition": dry_comp,
            })

        if x_h2 <= 1e-6 or x_co <= 1e-6:
            summary = _format_eval_summary(case, psr_converged, psr_temp, x_o2, x_ch4, x_co, x_h2, ratio, VALIDATION_FAIL_PENALTY, "INVALID_RATIO")
            print(summary)
            return Penalty(VALIDATION_FAIL_PENALTY, "INVALID_RATIO", summary, {
                "case": case,
                "dry_composition": dry_comp,
            })

        ratio = x_h2 / max(x_co, 1e-12)
        case_penalty = 1.0e4 * abs(ratio - 2.0)
        feasibility_penalty += case_penalty

        h2co_vals.append(ratio)
        ch4_vals.append(x_ch4)
        t_vals.append(t_out_k)

        def get_col(df: pd.DataFrame, name: str) -> np.ndarray:
            return df[name].values if name in df.columns else np.array([])

        rmse_T = _rmse(get_col(df_red, 'T_C'), get_col(df_base, 'T_C'))
        rmse_H2 = _rmse(get_col(df_red, 'X_H2'), get_col(df_base, 'X_H2'))
        rmse_CO = _rmse(get_col(df_red, 'X_CO'), get_col(df_base, 'X_CO'))

        per_case[case] = {
            'H2_CO_out': ratio,
            'CH4_out': x_ch4,
            'T_out_K': t_out_k,
            'psr_converged': psr_converged,
            'psr_temperature_K': psr_temp,
            'feasibility_penalty': case_penalty,
            'dry_X': {'H2': x_h2, 'CO': x_co, 'O2': x_o2, 'CH4': x_ch4},
            'rmse': {'T': rmse_T, 'X_H2': rmse_H2, 'X_CO': rmse_CO}
        }

    metrics = {
        'H2_CO_out_avg': float(np.mean(h2co_vals)) if h2co_vals else float('nan'),
        'CH4_slip_pct_avg': float(np.mean(ch4_vals) * 100.0) if ch4_vals else float('nan'),
        'T_out_avg': float(np.mean(t_vals)) if t_vals else float('nan'),
        'profile_rmse': {},
        'per_case': per_case,
        'feasibility_penalty': feasibility_penalty,
    }
    metrics['mech_yaml'] = str(mech_path.resolve())
    try:
        metrics['n_species'] = int(model.thermo.gas.n_species)  # type: ignore[attr-defined]
        metrics['n_reactions'] = int(model.thermo.gas.n_reactions)  # type: ignore[attr-defined]
    except AttributeError:
        pass

    with open(outdir / 'metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    return metrics


def calculate_metrics_penalty(metrics: Dict[str, Any], alpha: float = 8.0, n_species: int = None, n_reactions: int = None) -> float:
    """
    Calculate KPI error (percent) combining outlet composition, temperature, CH4 slip and PSR/PFR penalties.

    Args:
        metrics: Validation metrics dictionary
        alpha: Penalty weight for CH4 slip above 10% (fractional)
        n_species: Number of species in candidate mechanism (optional)
        n_reactions: Number of reactions in candidate mechanism (optional)

    Returns:
        Positive KPI error expressed as percentage (lower is better).
    """
    h2co_avg = float(metrics.get('H2_CO_out_avg', 0.0) or 0.0)
    ch4_slip_pct = float(metrics.get('CH4_slip_pct_avg', 0.0) or 0.0)
    t_out_avg = float(metrics.get('T_out_avg', 0.0) or 0.0)

    h2co_target = 2.0

    per_case_metrics = metrics.get('per_case', {}) or {}
    temp_error_terms: List[float] = []
    temp_targets: List[float] = []
    for case_name, case_data in per_case_metrics.items():
        raw_temp = case_data.get('T_out_K')
        try:
            case_temp = float(raw_temp)
        except (TypeError, ValueError):
            continue
        try:
            case_target = CaseConstants.get_target_temperature(case_name)
        except ValueError:
            continue
        if not np.isfinite(case_temp) or not np.isfinite(case_target):
            continue
        temp_targets.append(case_target)
        temp_error_terms.append(abs(case_temp - case_target) / max(case_target, 1e-6))

    if temp_error_terms:
        t_out_target = float(np.mean(temp_targets))
        t_error = float(np.mean(temp_error_terms))
    else:
        fallback_target = CaseConstants.get_target_temperature('A_case1_richter')
        t_out_target = fallback_target
        t_error = abs(t_out_avg - fallback_target) / max(fallback_target, 1e-6)

    h2co_error = abs(h2co_avg - h2co_target) / max(h2co_target, 1e-6)

    ch4_slip_fraction = max(0.0, ch4_slip_pct / 100.0)
    ch4_penalty = max(0.0, ch4_slip_fraction - 0.10)

    psr_penalty = 0.0
    per_case = metrics.get('per_case', {})
    for case_name, case_data in per_case.items():
        case_t_out = float(case_data.get('T_out_K', 0.0) or 0.0)
        if case_t_out < 1200.0:
            psr_penalty += 0.5
        rmse = case_data.get('rmse', {})
        if isinstance(rmse, dict):
            for key in ('X_H2', 'X_CO'):
                val = rmse.get(key)
                if val is not None and not np.isnan(val):
                    psr_penalty += max(0.0, float(val) - 0.10)

    size_bias = 0.0
    if n_species is not None:
        size_bias += 0.02 * abs(n_species - 30)
    if n_reactions is not None:
        size_bias += 0.005 * abs(n_reactions - 200)

    feas_penalty = float(metrics.get('feasibility_penalty', 0.0) or 0.0)
    kpi_error = 100.0 * (h2co_error + 0.5 * t_error + alpha * ch4_penalty + psr_penalty + size_bias) + feas_penalty
    return max(0.0, kpi_error)


def validate_mechanism(mech_path: Path, cases_root: Path, cases: List[str], outdir: Path) -> int:
    metrics = evaluate(mech_path, cases_root, cases, outdir)
    if isinstance(metrics, Penalty):
        payload = {
            'penalty': metrics.value,
            'reason': metrics.reason,
            'summary': metrics.summary,
        }
        with open(outdir / 'metrics.json', 'w') as f:
            json.dump(payload, f, indent=2)
        print(metrics.summary)
        print(f"Fitness score: {metrics.value:.4f} (feasibility fail: {metrics.reason})")
        return 1
    fitness = calculate_metrics_penalty(metrics)

    metrics['fitness_score'] = fitness
    with open(outdir / 'metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    
    print(f"Validation complete: {outdir.resolve().as_posix()}")
    print(f"Fitness score: {fitness:.4f}")
    return 0


def _mask_to_indices(mask: Iterable[bool], expected_len: int) -> List[int]:
    mask_list = list(mask)
    if len(mask_list) != expected_len:
        raise ValueError(f"Mask length {len(mask_list)} does not match expected {expected_len}")
    return [i for i, flag in enumerate(mask_list) if bool(flag)]


def _case_tolerances(case_name: str) -> Tuple[float, float]:
    try:
        target_k = CaseConstants.get_target_temperature(case_name)
        target_c = target_k - 273.15
        temp_tol = max(PSR_TEMP_TOL_K, 0.02 * abs(target_c))
    except ValueError:
        temp_tol = PSR_TEMP_TOL_K

    upper = case_name.upper()
    species_tol = PSR_SPECIES_TOL if "PSR" in upper else PFR_SPECIES_TOL
    return temp_tol, species_tol


def calculate_fitness(
    base_mech_path: str,
    species_mask: List[bool],
    reaction_mask: List[bool],
    cases_root: str,
    case_list_fast: List[str],
    case_list_full: Optional[List[str]] = None,
    soft_species: Optional[Set[str]] = None,
    target_size: int = 30,
    workdir: Optional[Path] = None,
) -> Tuple[float, Dict[str, Any]]:
    """Evaluate a candidate genome and return (penalty, diagnostics)."""
    import cantera as ct

    base_path = Path(base_mech_path)
    if not base_path.exists():
        raise FileNotFoundError(f"Base mechanism not found: {base_mech_path}")

    sol = ct.Solution(str(base_path))
    base_species_names = list(sol.species_names)
    base_reactions = list(sol.reactions())

    n_species_base = len(base_species_names)
    n_reactions_base = len(base_reactions)

    if len(species_mask) != n_species_base:
        raise ValueError(f"species_mask length {len(species_mask)} != base species count {n_species_base}")
    if len(reaction_mask) != n_reactions_base:
        raise ValueError(f"reaction_mask length {len(reaction_mask)} != base reaction count {n_reactions_base}")

    species_idx = _mask_to_indices(species_mask, n_species_base)
    reaction_idx = _mask_to_indices(reaction_mask, n_reactions_base)

    tmp_root = workdir or Path("runs") / "tmp_ga"
    tmp_root.mkdir(parents=True, exist_ok=True)
    run_tag = uuid.uuid4().hex[:12]
    cand_dir = tmp_root / f"cand_{run_tag}"
    cand_dir.mkdir(parents=True, exist_ok=True)
    cand_yaml = cand_dir / "candidate.yaml"

    closure_stats = write_reduced_mechanism_ct(
        base_mech_path=str(base_path),
        kept_reaction_idx=reaction_idx,
        out_mech_path=str(cand_yaml),
        kept_species_idx=species_idx,
        max_closure_iters=12,
    )

    eval_root = Path(cases_root)
    if not eval_root.exists():
        raise FileNotFoundError(f"cases_root not found: {cases_root}")

    fast_eval_dir = cand_dir / "fast_eval"
    diagnostics: Dict[str, Any] = {
        'candidate_yaml': str(cand_yaml),
        'closure': closure_stats,
        'fast_metrics_path': str(fast_eval_dir / 'metrics.json'),
        'cases': {},
        'psr_pass': True,
        'pfr_pass': True,
    }

    try:
        metrics_fast = evaluate(
            mech_path=cand_yaml,
            cases_root=eval_root,
            cases=case_list_fast,
            outdir=fast_eval_dir,
        )
    except Exception as exc:
        diagnostics['psr_pass'] = False
        diagnostics['pfr_pass'] = False
        diagnostics['error'] = str(exc)
        return VALIDATION_FAIL_PENALTY, diagnostics

    if isinstance(metrics_fast, Penalty):
        diagnostics['psr_pass'] = False
        diagnostics['pfr_pass'] = False
        diagnostics['failure_reason'] = metrics_fast.reason
        diagnostics['failure_summary'] = metrics_fast.summary
        diagnostics['feasibility_penalty'] = metrics_fast.value
        return float(metrics_fast.value), diagnostics

    total_feasibility_penalty = float(metrics_fast.get('feasibility_penalty', 0.0) or 0.0)


    def analyse_metrics(metric_blob: Dict[str, Any], case_list: List[str], pass_holder: Dict[str, bool]) -> float:
        total_pen = 0.0
        per_case = metric_blob.get('per_case', {})
        for case in case_list:
            case_data = per_case.get(case, {})
            rmse = case_data.get('rmse', {})
            t_rmse = float(rmse.get('T', float('inf')))
            h2_rmse = float(rmse.get('X_H2', float('inf')))
            co_rmse = float(rmse.get('X_CO', float('inf')))

            temp_tol, species_tol = _case_tolerances(case)
            is_psr = "PSR" in case.upper()
            is_pfr = "PFR" in case.upper()

            temp_fail = not np.isfinite(t_rmse) or t_rmse > temp_tol
            species_fail = not np.isfinite(h2_rmse) or h2_rmse > species_tol or not np.isfinite(co_rmse) or co_rmse > species_tol

            if is_psr:
                pass_holder['psr'] &= not (temp_fail or species_fail)
            if is_pfr or not is_psr:
                pass_holder['pfr'] &= not species_fail

            diagnostics['cases'][case] = {
                'rmse_T': t_rmse,
                'rmse_X_H2': h2_rmse,
                'rmse_X_CO': co_rmse,
                'pass': not (temp_fail or species_fail),
            }

            temp_excess = max(0.0, t_rmse - temp_tol) / max(temp_tol, 1e-6)
            h2_excess = max(0.0, h2_rmse - species_tol) / max(species_tol, 1e-6)
            co_excess = max(0.0, co_rmse - species_tol) / max(species_tol, 1e-6)
            total_pen += temp_excess + h2_excess + co_excess
        return total_pen

    pass_flags = {'psr': True, 'pfr': True}
    fast_penalty = analyse_metrics(metrics_fast, case_list_fast, pass_flags)

    metrics_full = None
    if case_list_full:
        full_eval_dir = cand_dir / "full_eval"
        metrics_full = evaluate(
            mech_path=cand_yaml,
            cases_root=eval_root,
            cases=case_list_full,
            outdir=full_eval_dir,
        )
        if isinstance(metrics_full, Penalty):
            diagnostics['psr_pass'] = False
            diagnostics['pfr_pass'] = False
            diagnostics['failure_reason'] = metrics_full.reason
            diagnostics['failure_summary'] = metrics_full.summary
            diagnostics['feasibility_penalty'] = metrics_full.value
            return float(metrics_full.value), diagnostics
        diagnostics['full_metrics_path'] = str(full_eval_dir / 'metrics.json')
        fast_penalty += analyse_metrics(metrics_full, case_list_full, pass_flags)
        total_feasibility_penalty += float(metrics_full.get('feasibility_penalty', 0.0) or 0.0)

    diagnostics['psr_pass'] = pass_flags['psr']
    diagnostics['pfr_pass'] = pass_flags['pfr']

    n_species_final = closure_stats['n_species']
    n_reactions_final = closure_stats['n_reactions']

    ch4_slip_pct = float(metrics_fast.get('CH4_slip_pct_avg', 0.0) or 0.0)
    ch4_penalty = max(0.0, ch4_slip_pct / 100.0 - CH4_SLIP_TOL)

    size_offset = (n_species_final - target_size) / max(target_size, 1)
    size_penalty = SIZE_BIAS_WEIGHT * (size_offset ** 2)

    if soft_species:
        kept_species_names = {base_species_names[i] for i in closure_stats['kept_species_idx']}
        dropped_soft = [s for s in soft_species if s in base_species_names and s not in kept_species_names]
        soft_penalty = SOFT_SPECIES_WEIGHT * len(dropped_soft)
    else:
        soft_penalty = 0.0

    aggregate_penalty = (
        100.0 * fast_penalty
        + 100.0 * ch4_penalty
        + size_penalty
        + soft_penalty
        + total_feasibility_penalty
    )

    if not diagnostics['psr_pass'] or not diagnostics['pfr_pass']:
        aggregate_penalty = VALIDATION_FAIL_PENALTY

    diagnostics.update({
        'n_species_final': n_species_final,
        'n_reactions_final': n_reactions_final,
        'ch4_slip_pct': ch4_slip_pct,
        'size_penalty': size_penalty,
        'soft_penalty': soft_penalty,
        'feasibility_penalty': total_feasibility_penalty,
        'raw_penalty': aggregate_penalty,
    })

    return float(aggregate_penalty), diagnostics
