#!/usr/bin/env python3
"""
Test Cantera wall heat flux sign convention.
"""

import cantera as ct
import numpy as np

def test_heat_flux_sign():
    print("=== CANTERA HEAT FLUX SIGN TEST ===")
    
    # Create simple gas
    gas = ct.Solution('hp_pox_0d_1d/gri30.yaml')
    gas.TPX = 1000.0, 101325.0, 'CH4:0.1, O2:0.2, N2:0.7'
    
    # Create reactor
    reactor = ct.IdealGasReactor(gas)
    reactor.volume = 0.001  # 1 L
    
    # Create reservoir (ambient)
    reservoir = ct.Reservoir(gas)
    
    # Create wall
    wall = ct.Wall(reactor, reservoir, A=1.0)  # 1 m² area
    
    # Create network
    net = ct.ReactorNet([reactor])
    
    print(f"Initial T = {reactor.T:.1f} K")
    
    # Test positive heat flux (should heat reactor)
    wall.heat_flux = +1000.0  # +1000 W/m²
    print(f"Applied +1000 W/m² heat flux")
    
    # Integrate for a short time
    net.advance(0.01)  # 10 ms
    T_plus = reactor.T
    print(f"After +1000 W/m²: T = {T_plus:.1f} K")
    
    # Reset and test negative heat flux (should cool reactor)
    reactor.thermo.TPX = 1000.0, 101325.0, 'CH4:0.1, O2:0.2, N2:0.7'
    wall.heat_flux = -1000.0  # -1000 W/m²
    print(f"Applied -1000 W/m² heat flux")
    
    net.advance(0.01)  # 10 ms
    T_minus = reactor.T
    print(f"After -1000 W/m²: T = {T_minus:.1f} K")
    
    # Analysis
    dT_plus = T_plus - 1000.0
    dT_minus = T_minus - 1000.0
    
    print(f"\n=== RESULTS ===")
    print(f"+1000 W/m² -> dT = {dT_plus:+.2f} K")
    print(f"-1000 W/m² -> dT = {dT_minus:+.2f} K")
    
    if dT_plus > 0 and dT_minus < 0:
        print("PASS: SIGN CONVENTION CORRECT: +flux heats, -flux cools")
        return True
    elif dT_plus < 0 and dT_minus > 0:
        print("FAIL: SIGN CONVENTION REVERSED: +flux cools, -flux heats")
        print("   Need to flip the sign in PSR implementation!")
        return False
    else:
        print("WARNING: UNCLEAR SIGN CONVENTION")
        return False

if __name__ == '__main__':
    test_heat_flux_sign()
