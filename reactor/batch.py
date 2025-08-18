"""Batch reactor integrations used throughout the project."""

import cantera as ct
import numpy as np
from dataclasses import dataclass
from typing import Optional

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
    X_or_Y0: dict,
    tf: float,
    nsteps: int = 1000,
    detect_ignition: bool = True,
    use_mole: bool = True,
    log_times: bool = False,
) -> BatchResult:
    """Integrate an adiabatic constant-pressure reactor.

    The integration samples the state at ``nsteps`` equally spaced times up to
    ``tf`` (or geometrically spaced when ``log_times`` is ``True``).  The
    initial state at ``t=0`` is always recorded.
    """

    # Set consistent fractions
    if use_mole:
        gas.TPX = T0, P0, X_or_Y0  # mole fractions
    else:
        gas.TPY = T0, P0, X_or_Y0  # mass fractions

    r = ct.IdealGasConstPressureReactor(gas, energy="on")
    net = ct.ReactorNet([r])
    net.atol, net.rtol = 1e-18, 1e-8

    # record t=0 **before** any advance
    times = [0.0]
    temps = [r.T]
    Ys = [r.thermo.Y.copy()]

    # advance on fixed times (optionally log-spaced)
    if log_times:
        ts = np.geomspace(1e-12, tf, nsteps)
        ts = np.insert(ts, 0, 0.0)
    else:
        dt = tf / nsteps
        ts = np.linspace(0.0, tf, nsteps + 1)

    for t in ts[1:]:
        net.advance(t)
        times.append(t)
        temps.append(r.T)
        Ys.append(r.thermo.Y.copy())

    times = np.array(times)
    temps = np.array(temps)
    Ys = np.vstack(Ys)

    # ignition delay (max dT/dt)
    delay = None
    if detect_ignition:
        dTdt = np.gradient(temps, times, edge_order=2)
        kmax = int(np.argmax(dTdt))
        delay = float(times[kmax])

    # inert-run guard
    dT = float(temps[-1] - temps[0])
    dY = float(np.max(np.abs(Ys[-1] - Ys[0])))
    if dT < 1.0 or dY < 1e-5:
        raise RuntimeError(f"no_reaction: dT={dT:.3g}, max|ΔY|={dY:.3g}")

    return BatchResult(times, temps, Ys, ignition_delay=delay)


def run_constant_volume(
    gas: ct.Solution,
    T0: float,
    P0: float,
    X_or_Y0: dict,
    tf: float,
    nsteps: int = 1000,
    detect_ignition: bool = True,
    use_mole: bool = True,
    log_times: bool = False,
) -> BatchResult:
    """Integrate an adiabatic constant-volume reactor."""

    if use_mole:
        gas.TPX = T0, P0, X_or_Y0
    else:
        gas.TPY = T0, P0, X_or_Y0

    r = ct.IdealGasReactor(gas, energy="on")
    net = ct.ReactorNet([r])
    net.atol, net.rtol = 1e-18, 1e-8

    times = [0.0]
    temps = [r.T]
    Ys = [r.thermo.Y.copy()]

    if log_times:
        ts = np.geomspace(1e-12, tf, nsteps)
        ts = np.insert(ts, 0, 0.0)
    else:
        dt = tf / nsteps
        ts = np.linspace(0.0, tf, nsteps + 1)

    for t in ts[1:]:
        net.advance(t)
        times.append(t)
        temps.append(r.T)
        Ys.append(r.thermo.Y.copy())

    times = np.array(times)
    temps = np.array(temps)
    Ys = np.vstack(Ys)

    delay = None
    if detect_ignition:
        dTdt = np.gradient(temps, times, edge_order=2)
        kmax = int(np.argmax(dTdt))
        delay = float(times[kmax])

    dT = float(temps[-1] - temps[0])
    dY = float(np.max(np.abs(Ys[-1] - Ys[0])))
    if dT < 1.0 or dY < 1e-5:
        raise RuntimeError(f"no_reaction: dT={dT:.3g}, max|ΔY|={dY:.3g}")

    return BatchResult(times, temps, Ys, ignition_delay=delay)
