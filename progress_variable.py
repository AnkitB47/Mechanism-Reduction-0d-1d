import numpy as np
from typing import Sequence


def progress_variable(Y: np.ndarray, weights: Sequence[float]) -> np.ndarray:
    """Compute progress variable PV(t) = sum_i w_i * Y_i(t)."""
    W = np.array(weights)
    return (Y * W[None, :]).sum(axis=1)
