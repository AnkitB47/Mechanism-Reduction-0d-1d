"""
Data collection module for mechanism reduction training dataset.

Collects PSR and PFR states across various conditions for training
the reduction algorithms (DRGEP, DRGASA, CSP/QSS).
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import numpy as np
import cantera as ct
import json

from . import ReductionConfig


def collect_training_data(base_mech: ct.Solution, config: ReductionConfig) -> Dict[str, Any]:
    """
    Collect comprehensive training dataset for reduction algorithms.

    Args:
        base_mech: Base mechanism solution
        config: Reduction configuration

    Returns:
        Dictionary containing PSR and PFR training states
    """

    training_data = {
        'psr_states': [],
        'pfr_states': [],
        'metadata': {
            'base_mechanism': str(config.base_mechanism),
            'protected_species': config.protected_species,
            'sample_conditions': len(config.psr_conditions) + len(config.pfr_cases) * len(config.sample_points)
        }
    }

    if config.verbose:
        print(f"Collecting training data from {len(config.psr_conditions)} PSR conditions and {len(config.pfr_cases)} PFR cases...")

    # Collect PSR states
    for i, psr_cond in enumerate(config.psr_conditions):
        try:
            psr_state = _collect_psr_state(base_mech, psr_cond, i)
            training_data['psr_states'].append(psr_state)
            if config.verbose and (i + 1) % 5 == 0:
                print(f"  Collected {i+1}/{len(config.psr_conditions)} PSR states")
        except Exception as e:
            if config.verbose:
                print(f"  Warning: Failed to collect PSR state {i}: {e}")

    # Collect PFR states
    for case_name in config.pfr_cases:
        try:
            pfr_states = _collect_pfr_states(base_mech, case_name, config, len(training_data['pfr_states']))
            training_data['pfr_states'].extend(pfr_states)
            if config.verbose:
                print(f"  Collected {len(pfr_states)} states from PFR case: {case_name}")
        except Exception as e:
            if config.verbose:
                print(f"  Warning: Failed to collect PFR states for {case_name}: {e}")

    if config.verbose:
        print(f"Training data collection complete: {len(training_data['psr_states'])} PSR states, {len(training_data['pfr_states'])} PFR states")

    return training_data


def _collect_psr_state(base_mech: ct.Solution, condition: Dict[str, float], index: int) -> Dict[str, Any]:
    """Collect a single PSR state."""

    # Create PSR reactor
    gas = base_mech
    reactor = ct.IdealGasConstPressureReactor(gas)

    # Set inlet conditions
    T_in = condition['T']
    P = condition['P'] * ct.one_atm
    phi = condition['phi']
    tau = condition['tau']

    # Create inlet mixture (CH4/Air for now - should be generalized)
    # For POX, need to set proper fuel/oxidizer composition
    # This is a placeholder - need proper POX inlet creation
    gas.TPX = T_in, P, 'CH4:1, O2:0.5, N2:1.887'  # phi ≈ 2.5

    # Set reactor volume based on residence time
    # V = tau * Q / (P / (R * T)) but simplified for now
    reactor.volume = tau * 0.001  # Placeholder

    # Solve to steady state
    net = ct.ReactorNet([reactor])

    # Advance to steady state
    try:
        net.advance_to_steady_state(rtol=1e-6, atol=1e-12)
    except Exception:
        # If steady state fails, try time advance
        net.advance(tau * 10)

    # Extract state
    state = {
        'index': index,
        'condition': condition,
        'T': reactor.T,
        'P': reactor.P,
        'composition': dict(reactor.thermo.Y),
        'rop': dict(reactor.kinetics.net_rates_of_progress),
        'source': 'PSR'
    }

    return state


def _collect_pfr_states(base_mech: ct.Solution, case_name: str, config: ReductionConfig, start_index: int) -> List[Dict[str, Any]]:
    """Collect PFR states for a specific case."""

    states = []

    # This would integrate with the existing PFR simulation code
    # For now, return placeholder states

    for i, z in enumerate(config.sample_points):
        state = {
            'index': start_index + i,
            'case': case_name,
            'z': z,
            'T': 1000 + 200 * z,  # Placeholder temperature profile
            'P': ct.one_atm,
            'composition': {'H2': 0.1, 'CO': 0.1, 'CO2': 0.05, 'CH4': 0.01, 'H2O': 0.1, 'O2': 0.01, 'N2': 0.63},  # Placeholder
            'rop': {},  # Would need ROP computation
            'source': 'PFR'
        }
        states.append(state)

    return states
