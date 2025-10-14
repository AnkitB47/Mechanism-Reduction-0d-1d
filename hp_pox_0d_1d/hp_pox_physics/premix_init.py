# hp_pox_physics/premix_init.py
"""
Premix ignition initialization for HP-POX PSR.
Creates HP-equilibrated initial states for PSR ignition kernel.
"""

import cantera as ct
from typing import Dict, Tuple


def _q_from(T: float, p: float, X_or_Y: Dict[str, float], mech: str, mass: float, basis: str = "X") -> ct.Quantity:
    """Create a Cantera Quantity from T, p, composition, mechanism, and mass."""
    gas = ct.Solution(mech)
    if basis.upper() == "X":
        gas.TPX = T, p, X_or_Y
    else:
        gas.TPY = T, p, X_or_Y
    q = ct.Quantity(gas, constant="HP")
    q.mass = mass
    return q


def premix_equilibrated_initial_state(
    mech: str,
    p_pa: float,
    # fuel premix: natural gas + N2 (diluent)
    T_fuel_K: float, fuel_comp_X: Dict[str, float], m_fuel_kg_s: float,
    T_N2_K: float, m_N2_kg_s: float,
    # oxidizer premix: O2 + secondary steam
    T_O2_K: float, m_O2_kg_s: float,
    T_secsteam_K: float, m_secsteam_kg_s: float,
    # primary steam (added after ignition kernel is formed)
    T_primsteam_K: float, m_primsteam_kg_s: float,
    # fraction of (fuel+N2) that actually forms the flame kernel initially
    y_feed_to_flame: float = 1.0
) -> Tuple[ct.Solution, ct.Solution]:
    """Return (gas_init_flame, gas_psr_in). All states set on `mech` with HP mixing."""
    # 1) Build premixed fuel (NG+N2) and premixed oxidizer (O2+H2Osec)
    q_fuel = _q_from(T_fuel_K, p_pa, fuel_comp_X, mech, m_fuel_kg_s, basis="X")
    q_n2   = _q_from(T_N2_K,   p_pa, {"N2": 1.0},    mech, m_N2_kg_s,  basis="X")
    q_fuel_premix = q_fuel + q_n2

    q_o2   = _q_from(T_O2_K,       p_pa, {"O2": 1.0},  mech, m_O2_kg_s,       basis="X")
    q_h2oS = _q_from(T_secsteam_K, p_pa, {"H2O": 1.0}, mech, m_secsteam_kg_s, basis="X")
    q_oxid = q_o2 + q_h2oS

    # 2) Portion of fuel premix participating in ignition kernel
    # Create new gas state with scaled mass
    gas_feed = ct.Solution(mech)
    gas_feed.TPY = q_fuel_premix.T, q_fuel_premix.P, q_fuel_premix.Y
    q_feed_to_flame = ct.Quantity(gas_feed, constant="HP")
    q_feed_to_flame.mass = q_fuel_premix.mass * y_feed_to_flame

    # 3) Initial flame mixture = (fuel premix portion) + (oxidizer premix)
    q_init = q_feed_to_flame + q_oxid
    gas_init = ct.Solution(mech)
    gas_init.TPY = q_init.TPY
    gas_init.equilibrate("HP")  # HP-equilibrated ignition seed

    # Cap temperature to prevent numerical issues
    if gas_init.T > 3000.0:
        gas_init.TP = 3000.0, gas_init.P

    # 4) Add primary steam to the equilibrated kernel to form the PSR inlet state
    q_prim = _q_from(T_primsteam_K, p_pa, {"H2O": 1.0}, mech, m_primsteam_kg_s, basis="X")
    q_total_to_psr = ct.Quantity(gas_init, constant="HP") + q_prim
    gas_psr_in = ct.Solution(mech)
    gas_psr_in.TPY = q_total_to_psr.TPY

    return gas_init, gas_psr_in
