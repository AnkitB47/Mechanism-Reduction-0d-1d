#!/usr/bin/env python3
"""
Common I/O utilities for plotting and data handling.
"""

import matplotlib
from pathlib import Path
import json
import pandas as pd
from typing import Dict, Tuple, Any, Union


def ensure_noninteractive_backend() -> None:
    """Set matplotlib to non-interactive backend safely."""
    try:
        matplotlib.use("Agg")
    except Exception:
        pass  # Already set or failed to set


def read_axial_csv(path: Union[str, Path]) -> Tuple[pd.DataFrame, str]:
    """
    Read axial profiles CSV and return (dataframe, z_column_name).
    
    Tolerates both 'z_m' and 'z' column names.
    
    Args:
        path: Path to CSV file
        
    Returns:
        Tuple of (dataframe, z_column_name)
        
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    df = pd.read_csv(path)
    
    # Determine z column name
    if 'z_m' in df.columns:
        z_col = 'z_m'
    elif 'z' in df.columns:
        z_col = 'z'
    else:
        raise ValueError(f"No z column found in {path}. Expected 'z_m' or 'z'.")
    
    return df, z_col


def safe_mkdir(path: Union[str, Path]) -> None:
    """Create directory with parents if needed."""
    Path(path).mkdir(parents=True, exist_ok=True)


def load_kpis(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load KPIs JSON file.
    
    Args:
        path: Path to JSON file
        
    Returns:
        Dictionary with KPIs data, empty dict if file missing/invalid
    """
    try:
        with open(path, 'r') as f:
            return json.load(f)
    except Exception:
        return {}


def dry_value(kpis: Dict[str, Any], key: str) -> float:
    """
    Extract dry basis value from KPIs.
    
    Args:
        kpis: KPIs dictionary
        key: Species key (e.g., 'H2', 'CO')
        
    Returns:
        Dry basis mole fraction, 0.0 if not available
    """
    if not kpis.get('dry_basis_computed', False):
        return 0.0
    
    sim_dry = kpis.get('sim_dry', {})
    return float(sim_dry.get(key, 0.0))


def perr(a: float, b: float) -> float:
    """
    Calculate percentage error: 100 * (b - a) / a
    
    Args:
        a: Reference value
        b: Comparison value
        
    Returns:
        Percentage error, float('nan') if a == 0
    """
    if a == 0:
        return float('nan')
    return 100.0 * (b - a) / a
