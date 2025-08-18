import cantera as ct
import numpy as np
from mechanism.mix import methane_air_mole_fractions, mole_to_mass_fractions
from reactor.batch import run_constant_pressure

def test_hot_methane_run():
    gas = ct.Solution('gri30.yaml')
    x0 = methane_air_mole_fractions(1.0)
    Y0 = mole_to_mass_fractions(gas, x0)
    res = run_constant_pressure(gas, 1500.0, ct.one_atm, Y0, 0.2, nsteps=200, use_mole=False)
    assert (res.temperature[-1] - res.temperature[0]) > 50.0
    assert res.ignition_delay is not None and res.ignition_delay < 0.2
    assert np.max(np.abs(res.mass_fractions[-1] - res.mass_fractions[0])) > 1e-4
