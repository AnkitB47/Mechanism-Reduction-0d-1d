import cantera as ct
import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple, Optional

@dataclass
class BatchResult:
    """Result of a batch reactor simulation."""

    time: np.ndarray
    temperature: np.ndarray
    mass_fractions: np.ndarray
    ignition_delay: Optional[float] = None


def run_constant_pressure(
    gas: ct.Solution,
    T0: float,
    P0: float,
    X0: dict,
    tf: float,
    nsteps: int = 1000,
    detect_ignition: bool = True,
) -> BatchResult:
    """Integrate an adiabatic constant-pressure reactor.

    Parameters
    ----------
    gas:
        Cantera ``Solution`` object representing the mechanism.
    T0, P0:
        Initial temperature and pressure.
    X0:
        Initial composition specified as mole fractions. Using mole fractions
        avoids accidental non-reactive mixtures when large values are passed
        for stoichiometric ratios.
    tf, nsteps:
        Final time and number of integration steps.
    detect_ignition:
        If ``True`` the ignition delay is estimated from the temperature
        derivative.
    """

    # Use TPX to ensure the provided composition is interpreted as mole
    # fractions.  Previous versions used TPY which silently normalised the
    # numbers as mass fractions and resulted in nearly non-reactive mixtures.
    gas.TPX = T0, P0, X0
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on")
    sim = ct.ReactorNet([reactor])

    time: list[float] = []
    T: list[float] = []
    Y: list[np.ndarray] = []
    delay = None
    last_dTdt = 0.0
    dt = tf / nsteps
    for _ in range(nsteps):
        sim.advance(sim.time + dt)
        time.append(sim.time)
        T.append(reactor.T)
        Y.append(reactor.Y)
        if detect_ignition:
            dTdt = (reactor.T - T[-2]) / dt if len(T) > 1 else 0.0
            if last_dTdt > 0 and dTdt <= 0 and delay is None:
                delay = time[-2]
            last_dTdt = dTdt

    return BatchResult(np.array(time), np.array(T), np.array(Y), ignition_delay=delay)


def run_constant_volume(
    gas: ct.Solution,
    T0: float,
    P0: float,
    Y0: dict,
    tf: float,
    nsteps: int = 1000,
    detect_ignition: bool = True,
) -> BatchResult:
    """Integrate a constant-volume adiabatic reactor."""

    # As with the constant-pressure case, interpret composition as mole
    # fractions to avoid near-inert mixtures when using stoichiometric ratios.
    gas.TPX = T0, P0, Y0
    reactor = ct.IdealGasReactor(gas, energy="on")
    sim = ct.ReactorNet([reactor])
    time = []
    T = []
    Y = []
    delay = None
    last_dTdt = 0.0
    dt = tf / nsteps
    for _ in range(nsteps):
        sim.advance(sim.time + dt)
        time.append(sim.time)
        T.append(reactor.T)
        Y.append(reactor.Y)
        if detect_ignition:
            dTdt = (reactor.T - T[-2]) / dt if len(T) > 1 else 0.0
            if last_dTdt > 0 and dTdt <= 0 and delay is None:
                delay = time[-2]
            last_dTdt = dTdt

    return BatchResult(np.array(time), np.array(T), np.array(Y), ignition_delay=delay)
