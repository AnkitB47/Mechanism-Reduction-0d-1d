"""
Physics data API for mechanism reduction training.

Exposes PSR/PFR samplers for DRGEP training data collection.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from typing import List, Dict, Any
import cantera as ct

# Import existing physics systems
from .hp_pox_model import HPPOXModel
from .config_cases import CaseConfig


class State:
    """Training state for mechanism reduction."""

    def __init__(self, T: float, P: float, phi: float, tau: float, composition: Dict[str, float], source: str = "PSR"):
        self.T = T
        self.P = P
        self.phi = phi
        self.tau = tau
        self.composition = composition
        self.source = source
        self.rop = {}  # Placeholder for ROP data

    def get_solution(self) -> ct.Solution:
        """Create Cantera solution for this state."""
        # Create a dummy solution for now
        # In practice, this would set up the proper mixture
        gas = ct.Solution("gri30.yaml")
        gas.TPX = self.T, self.P, self.composition
        return gas

    def get(self, key, default=None):
        """Dict-like interface for DRGEP."""
        return getattr(self, key, default)


def collect_training_states() -> List[State]:
    """
    Collect PSR and PFR training states for mechanism reduction.

    Returns:
        List of State objects for DRGEP training
    """

    states = []

    # PSR conditions for training
    psr_conditions = [
        {"T": 800, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 1000, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 1100, "phi": 2.5, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 2.0, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 3.0, "P": 1.0, "tau": 0.010},
        {"T": 900, "phi": 2.5, "P": 2.0, "tau": 0.010},
        {"T": 900, "phi": 2.5, "P": 5.0, "tau": 0.010},
    ]

    # Add PSR states (placeholder - would run actual PSR simulations)
    for cond in psr_conditions:
        state = State(
            T=cond["T"],
            P=cond["P"] * ct.one_atm,
            phi=cond["phi"],
            tau=cond["tau"],
            composition={'H2': 0.1, 'CO': 0.1, 'CO2': 0.05, 'CH4': 0.01, 'H2O': 0.1, 'O2': 0.01, 'N2': 0.63},
            source="PSR"
        )
        states.append(state)

    # Add PFR states (placeholder - would sample from actual PFR runs)
    pfr_cases = ["A_case1_richter", "A_case4_richter"]
    sample_points = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]

    for case in pfr_cases:
        for z in sample_points:
            state = State(
                T=1000 + 200 * z,
                P=ct.one_atm,
                phi=2.5,
                tau=0.010,
                composition={'H2': 0.1, 'CO': 0.1, 'CO2': 0.05, 'CH4': 0.01, 'H2O': 0.1, 'O2': 0.01, 'N2': 0.63},
                source="PFR"
            )
            states.append(state)

    return states
