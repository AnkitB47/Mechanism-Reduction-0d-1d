import cantera as ct
import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple

@dataclass
class BatchResult:
    time: np.ndarray
    temperature: np.ndarray
    mass_fractions: np.ndarray


def run_constant_pressure(gas: ct.Solution, T0: float, P0: float, Y0: dict, tf: float, nsteps: int = 1000) -> BatchResult:
    gas.TPY = T0, P0, Y0
    reactor = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([reactor])
    time = []
    T = []
    Y = []
    dt = tf/nsteps
    for _ in range(nsteps):
        sim.advance(sim.time + dt)
        time.append(sim.time)
        T.append(reactor.T)
        Y.append(reactor.Y)
    return BatchResult(np.array(time), np.array(T), np.array(Y))


def run_constant_volume(gas: ct.Solution, T0: float, P0: float, Y0: dict, tf: float, nsteps: int = 1000) -> BatchResult:
    gas.TPY = T0, P0, Y0
    reactor = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([reactor])
    time = []
    T = []
    Y = []
    dt = tf/nsteps
    for _ in range(nsteps):
        sim.advance(sim.time + dt)
        time.append(sim.time)
        T.append(reactor.T)
        Y.append(reactor.Y)
    return BatchResult(np.array(time), np.array(T), np.array(Y))
