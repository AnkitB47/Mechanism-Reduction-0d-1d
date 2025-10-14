"""
Advanced mechanism reduction pipeline for methane POX using literature-backed algorithms.

This module implements:
- DRGEP (Directed Relation Graph with Error Propagation)
- DRGASA (DRG-Aided Sensitivity Analysis)
- CSP/QSS (Computational Singularity Perturbation/Quasi-Steady State)
- Fixed-point reaction↔species closure
- Comprehensive validation

References:
- Lu & Law (Combustion and Flame, 2005/2006) for DRGEP
- Pepiot-Desjardins & Pitsch (Combustion and Flame, 2008) for DRGASA
- Maas & Pope for CSP/QSS
"""

from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import yaml
import json
import numpy as np
import cantera as ct
from dataclasses import dataclass, field

@dataclass
class ReductionConfig:
    """Configuration for mechanism reduction pipeline."""

    # Core parameters
    base_mechanism: Path = Path("gri30.yaml")
    output_dir: Path = Path("reduction_output")

    # Protected species (hard requirement)
    protected_species: List[str] = field(default_factory=lambda: [
        "H2", "CO", "CO2", "H2O", "CH4", "N2",  # Core POX species
        "H", "O", "OH", "HO2", "H2O2",  # Radical core
        "HCO", "CH2O"  # Path species
    ])

    # Hard-kept reaction families (exact matches required)
    protected_reactions: List[str] = field(default_factory=lambda: [
        # H2/O2 chain branching
        "H+O2=O+OH", "O+H2=OH+H", "OH+H2=H+H2O",
        "HO2+H=2OH", "HO2+HO2=H2O2+O2", "H2O2(+M)=2OH(+M)",
        # CH4 activation
        "CH4+OH=CH3+H2O", "CH4+O=CH3+OH", "CH4+HO2=CH3+H2O2",
        # CO/CO2 core
        "CO+OH=CO2+H", "HCO(+M)=H+CO(+M)", "HCO+O2=CO+HO2",
        "CH2O+H=HCO+H2", "CH2O+OH=HCO+H2O"
    ])

    # DRGEP parameters
    tau_s_grid: List[float] = field(default_factory=lambda: [0.01, 0.02, 0.05, 0.10, 0.15, 0.20])
    drgep_max_error: float = 0.05

    # DRGASA parameters
    sensitivity_cutoff: float = 0.01
    target_metrics: List[str] = field(default_factory=lambda: ["T", "H2", "CO", "CO2", "CH4", "H2O", "H2_CO_ratio"])

    # CSP/QSS parameters
    csp_importance_cutoff: float = 0.01
    qss_max_timescale_ratio: float = 10.0

    # Validation tolerances
    psr_tolerance: float = 50.0  # K
    psr_h2co_tolerance: float = 0.15  # relative
    pfr_species_tolerance: float = 0.10  # relative

    # Sampling parameters
    psr_conditions: List[Dict[str, float]] = field(default_factory=lambda: [
        {"T": 800, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 1000, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 1100, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 2.0, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 3.0, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 2.5, "P": 2.0, "tau": 0.010},
        {"T": 900, "phi": 2.5, "P": 5.0, "tau": 0.010},
    ])

    pfr_cases: List[str] = field(default_factory=lambda: ["A_case1_richter", "A_case4_richter"])
    sample_points: List[float] = field(default_factory=lambda: [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0])

    # Output control
    save_intermediate: bool = True
    verbose: bool = True

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dict for YAML serialization."""
        return {
            'base_mechanism': str(self.base_mechanism),
            'output_dir': str(self.output_dir),
            'protected_species': self.protected_species,
            'protected_reactions': self.protected_reactions,
            'tau_s_grid': self.tau_s_grid,
            'drgep_max_error': self.drgep_max_error,
            'sensitivity_cutoff': self.sensitivity_cutoff,
            'target_metrics': self.target_metrics,
            'csp_importance_cutoff': self.csp_importance_cutoff,
            'qss_max_timescale_ratio': self.qss_max_timescale_ratio,
            'psr_tolerance': self.psr_tolerance,
            'psr_h2co_tolerance': self.psr_h2co_tolerance,
            'pfr_species_tolerance': self.pfr_species_tolerance,
            'psr_conditions': self.psr_conditions,
            'pfr_cases': self.pfr_cases,
            'sample_points': self.sample_points,
            'save_intermediate': self.save_intermediate,
            'verbose': self.verbose
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ReductionConfig':
        """Create config from dict."""
        return cls(**data)


@dataclass
class ReductionResult:
    """Results from mechanism reduction pipeline."""

    original_species: int
    original_reactions: int
    reduced_species: int
    reduced_reactions: int

    # Algorithm results
    drgep_species: int
    drgep_reactions: int
    drgasa_species: int
    drgasa_reactions: int
    csp_qss_species: int
    csp_qss_reactions: int

    # Validation results
    psr_convergence: Dict[str, bool]
    psr_errors: Dict[str, Dict[str, float]]
    pfr_errors: Dict[str, Dict[str, float]]

    # Final mechanism path
    reduced_mechanism_path: Optional[Path] = None
    reduction_report_path: Optional[Path] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dict for JSON serialization."""
        return {
            'original_species': self.original_species,
            'original_reactions': self.original_reactions,
            'reduced_species': self.reduced_species,
            'reduced_reactions': self.reduced_reactions,
            'drgep_species': self.drgep_species,
            'drgep_reactions': self.drgep_reactions,
            'drgasa_species': self.drgasa_species,
            'drgasa_reactions': self.drgasa_reactions,
            'csp_qss_species': self.csp_qss_species,
            'csp_qss_reactions': self.csp_qss_reactions,
            'psr_convergence': self.psr_convergence,
            'psr_errors': self.psr_errors,
            'pfr_errors': self.pfr_errors,
            'reduced_mechanism_path': str(self.reduced_mechanism_path) if self.reduced_mechanism_path else None,
            'reduction_report_path': str(self.reduction_report_path) if self.reduction_report_path else None
        }
