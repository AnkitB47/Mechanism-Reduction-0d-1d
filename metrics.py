"""Distance metrics for reactor comparison."""

from __future__ import annotations

import numpy as np


def l2_distance(Y_full: np.ndarray, Y_red: np.ndarray, time: np.ndarray) -> float:
    """Integrate squared difference over time."""
    diff = Y_full - Y_red
    return float(np.trapz(np.sum(diff**2, axis=1), time))


def pv_error(full: np.ndarray, red: np.ndarray) -> float:
    """Relative L2 error between two progress variables."""
    return float(np.linalg.norm(full - red) / (np.linalg.norm(full) + 1e-12))


def ignition_delay(time: np.ndarray, T: np.ndarray) -> tuple[float, float]:
    """Return ignition delay and max dT/dt."""

    dTdt = np.gradient(T, time, edge_order=2)
    kmax = int(np.argmax(dTdt))
    return float(time[kmax]), float(dTdt[kmax])
