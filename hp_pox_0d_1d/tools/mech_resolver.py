"""
Mechanism resolver utilities for the HP-POX tooling.

Provides alias-to-path resolution so users can reference mechanisms by short names
such as ``gri30`` or ``aramco`` as well as by explicit file paths.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]

SEARCH_ROOTS: Tuple[Path, ...] = (
    REPO_ROOT / "mechanisms",
    REPO_ROOT / "mechanism",
    REPO_ROOT / "mechanisms_reduced",
    REPO_ROOT / "hp_pox_0d_1d" / "mechanisms",
    REPO_ROOT / "hp_pox_0d_1d" / "mechanisms_reduced",
    REPO_ROOT / "hp_pox",
    REPO_ROOT / "hp_pox_0d_1d",
)

ALIAS_MAP: Dict[str, Dict[str, Iterable[str]]] = {
    "gri30": {
        "aliases": {"gri", "gri30", "gri_30", "gri-30", "gri30mech"},
        "patterns": (
            "gri30.yaml",
            "gri30*.yaml",
            "gri/gri30*.yaml",
            "**/gri30*.yaml",
        ),
    },
    "aramco": {
        "aliases": {"aramco", "aramco30", "aramco3", "aramco2", "aramco2.0", "aramco_mech", "aramcomech"},
        "patterns": (
            "aramco*.yaml",
            "aramco/*aramco*.yaml",
            "**/aramco*.yaml",
        ),
    },
    "san_diego": {
        "aliases": {"sandiego", "san_diego", "sd", "sdmech", "sdt", "sdmechanism"},
        "patterns": (
            "san_diego*.yaml",
            "sandiego*.yaml",
            "san_diego/*sd*.yaml",
            "**/san_diego*.yaml",
            "**/sd*.yaml",
        ),
    },
}


def _normalize_alias(alias: str) -> str:
    return alias.lower().replace("-", "").replace("_", "")


def _candidate_paths(root: Path, pattern: str) -> List[Path]:
    if not root.exists():
        return []
    return [path for path in root.rglob(pattern) if path.is_file() and path.suffix.lower() == ".yaml"]


def _search_alias(alias_key: str) -> List[Path]:
    entry = ALIAS_MAP[alias_key]
    patterns = entry["patterns"]
    matches: List[Path] = []
    for root in SEARCH_ROOTS:
        for pattern in patterns:
            matches.extend(_candidate_paths(root, pattern))
    # Deduplicate while preserving order
    seen = set()
    unique_matches = []
    for path in matches:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_matches.append(resolved)
    return unique_matches


def list_available_mechanisms() -> Dict[str, List[Path]]:
    """Return a mapping of alias families to discovered mechanism files."""
    available: Dict[str, List[Path]] = {}
    for family in ALIAS_MAP:
        results = _search_alias(family)
        if results:
            available[family] = results
    return available


def resolve_mechanism(mech_arg: str) -> Path:
    """
    Resolve a mechanism argument to an existing YAML file.

    Parameters
    ----------
    mech_arg:
        Either a filesystem path (absolute or relative to the repo root) or a known alias
        such as ``gri30``, ``aramco``, or ``san_diego``.

    Returns
    -------
    Path
        The resolved mechanism path.

    Raises
    ------
    FileNotFoundError
        If the mechanism cannot be resolved.
    """
    candidate = Path(mech_arg)
    if candidate.suffix.lower() == ".yaml" or candidate.exists() or "/" in mech_arg or "\\" in mech_arg:
        resolved = candidate if candidate.is_absolute() else (REPO_ROOT / candidate)
        if resolved.exists():
            print(f"[mech_resolver] discovered mechanism -> {resolved.resolve()}")
            return resolved.resolve()
        raise FileNotFoundError(
            f"Mechanism file '{mech_arg}' not found. Checked: {resolved.resolve()}"
        )

    normalized = _normalize_alias(mech_arg)
    for family, entry in ALIAS_MAP.items():
        if normalized in {_normalize_alias(alias) for alias in entry["aliases"]}:
            matches = _search_alias(family)
            if matches:
                chosen = None
                for candidate in matches:
                    name = candidate.name.lower()
                    if "reduced" not in name and "tmp" not in name:
                        chosen = candidate
                        break
                if chosen is None:
                    chosen = matches[0]
                print(f"[mech_resolver] alias '{mech_arg}' resolved to {chosen}")
                return chosen
            probed = [str(root) for root in SEARCH_ROOTS]
            raise FileNotFoundError(
                f"Mechanism alias '{mech_arg}' not resolved. "
                f"Searched patterns {list(entry['patterns'])} under roots {probed}."
            )

    raise FileNotFoundError(
        f"Unknown mechanism alias '{mech_arg}'. "
        f"Known aliases: {sorted({alias for entry in ALIAS_MAP.values() for alias in entry['aliases']})}. "
        f"Provide a full path to a YAML file or one of the supported aliases."
    )
