import cantera as ct
import numpy as np
from typing import Dict


def methane_air_mole_fractions(phi: float, diluent="N2") -> Dict[str, float]:
    # Stoich CH4 + 2(O2 + 3.76N2) -> φ = (F/O)_actual / (F/O)_stoich
    # O2 needed per CH4 at stoich = 2
    # For a given φ, equivalence implies O2_actual = 2 / φ
    x = {
        "CH4": 1.0,
        "O2":  2.0 / phi,
    }
    x[diluent] = 3.76 * x["O2"]
    # normalize to mole fractions
    s = sum(x.values())
    return {k: v / s for k, v in x.items()}


def mole_to_mass_fractions(sol: ct.Solution, x: Dict[str, float]) -> Dict[str, float]:
    # Convert to Cantera’s mass fraction dict
    sol.TPX = 300.0, ct.one_atm, x
    return {sp: float(y) for sp, y in zip(sol.species_names, sol.Y) if y > 0}


# (Optional) HR presets from paper (fill exact numbers later if you have them)
def preset_HR2(sol: ct.Solution) -> Dict[str, float]:
    # Placeholder: plug in values extracted from the paper/table
    # Return a **mass-fraction** dict aligned to sol.species_names
    raise NotImplementedError
