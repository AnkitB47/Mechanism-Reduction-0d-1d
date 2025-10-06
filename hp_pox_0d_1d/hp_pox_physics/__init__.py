"""
HP-POX Physics Module

A clean implementation of the HP-POX 0D-1D reactor model with proper physics
and validation gates. Implements Richter-consistent behavior without artificial scaling.
"""

from .hp_pox_model import run_case, run_all

__all__ = ["run_case", "run_all"]
