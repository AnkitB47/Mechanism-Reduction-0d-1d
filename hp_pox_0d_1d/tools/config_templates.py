"""
Configuration templating helpers for HP-POX orchestration.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TEMPLATE_PATH = REPO_ROOT / "hp_pox_0d_1d" / "config.yaml"

BUILTIN_TEMPLATE = {
    "mechanism": "",
    "cases_root": "hp_pox_results_full",
    "output_root": "hp_pox_results_full",
    "cases": ["A_case1_richter"],
}


def write_full_config(
    template_src: Optional[Path],
    out_path: Path,
    mechanism_path: Path,
    cases_root: Optional[Path] = None,
    output_root: Optional[Path] = None,
) -> Path:
    """
    Write a configuration YAML for the HP-POX runner.

    Parameters
    ----------
    template_src:
        Optional source configuration to clone. If ``None`` or missing, falls back to
        ``hp_pox_0d_1d/config.yaml`` or a built-in minimal template.
    out_path:
        Destination path for the generated configuration.
    mechanism_path:
        Mechanism YAML to embed in the configuration.
    cases_root / output_root:
        Optional overrides for result directories.
    """
    template_path: Optional[Path] = None
    if template_src and Path(template_src).exists():
        template_path = Path(template_src)
    elif DEFAULT_TEMPLATE_PATH.exists():
        template_path = DEFAULT_TEMPLATE_PATH

    data = None
    if template_path:
        with template_path.open("r") as f:
            data = yaml.safe_load(f)

    if not isinstance(data, dict):
        data = BUILTIN_TEMPLATE.copy()

    data["mechanism"] = str(mechanism_path.resolve())

    if cases_root is not None:
        data["cases_root"] = str(cases_root.resolve())
    else:
        data.setdefault("cases_root", "hp_pox_results_full")

    if output_root is not None:
        data["output_root"] = str(output_root.resolve())
    else:
        data.setdefault("output_root", "hp_pox_results_full")

    validation = data.setdefault("validation", {})
    validation.setdefault("h2co_target", 2.0)
    validation.setdefault("h2co_rel_tol", 0.05)
    validation.setdefault("max_abs_dT_K", 50)
    validation.setdefault("require_hot_branch", True)
    validation.setdefault("require_psr_converged", True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        yaml.safe_dump(data, f, sort_keys=False)
    return out_path
