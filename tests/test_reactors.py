import cantera as ct
from mechanism.mix import methane_air_mole_fractions, mole_to_mass_fractions
from reactor.batch import run_constant_pressure


def test_ch4_consumed():
    gas = ct.Solution('gri30.yaml')
    x0 = methane_air_mole_fractions(1.0)
    Y0 = mole_to_mass_fractions(gas, x0)
    res = run_constant_pressure(gas, 1500.0, ct.one_atm, Y0, 0.05, nsteps=20, use_mole=False)
    idx = gas.species_index('CH4')
    assert res.mass_fractions[-1, idx] < res.mass_fractions[0, idx]
