#!/usr/bin/env python3
"""
Debug Cantera wall setup to understand why heat flux is not working.
"""

import cantera as ct
import numpy as np

def debug_wall_setup():
    print("=== CANTERA WALL SETUP DEBUG ===")
    
    # Create simple gas
    gas = ct.Solution('hp_pox_0d_1d/gri30.yaml')
    gas.TPX = 1000.0, 101325.0, 'CH4:0.1, O2:0.2, N2:0.7'
    
    print(f"Initial gas state: T={gas.T:.1f}K, P={gas.P:.0f}Pa")
    
    # Create reactor
    reactor = ct.IdealGasReactor(gas)
    reactor.volume = 0.001  # 1 L
    
    print(f"Reactor volume: {reactor.volume:.3f} m^3")
    print(f"Reactor T: {reactor.T:.1f} K")
    
    # Create reservoir (ambient) - SAME gas object!
    reservoir = ct.Reservoir(gas)
    
    print(f"Reservoir T: {reservoir.thermo.T:.1f} K")
    
    # Create wall
    wall = ct.Wall(reactor, reservoir, A=1.0)  # 1 m² area
    
    print(f"Wall area: {wall.area} m^2")
    print(f"Wall heat flux: {wall.heat_flux} W/m^2")
    
    # Check wall properties
    print(f"\nWall properties:")
    print(f"  area: {wall.area}")
    print(f"  heat_flux: {wall.heat_flux}")
    print(f"  Available attributes: {[attr for attr in dir(wall) if not attr.startswith('_')]}")
    
    # Try different wall setups
    print(f"\n=== TESTING DIFFERENT WALL SETUPS ===")
    
    # Test 1: Wall with different ambient temperature
    print(f"\nTest 1: Different ambient T")
    reservoir2 = ct.Reservoir(gas)
    reservoir2.thermo.TPX = 300.0, 101325.0, 'CH4:0.1, O2:0.2, N2:0.7'  # Cold ambient
    
    wall2 = ct.Wall(reactor, reservoir2, A=1.0)
    wall2.heat_flux = 1000.0  # +1000 W/m²
    
    net = ct.ReactorNet([reactor])
    print(f"Before integration: T = {reactor.T:.1f} K")
    net.advance(0.01)
    print(f"After +1000 W/m² with cold ambient: T = {reactor.T:.1f} K")
    
    # Test 2: Check if wall is actually connected
    print(f"\nTest 2: Wall connection check")
    print(f"Wall created successfully: {wall2 is not None}")
    
    # Test 3: Try heat transfer coefficient instead
    print(f"\nTest 3: Heat transfer coefficient")
    wall3 = ct.Wall(reactor, reservoir2, A=1.0)
    wall3.heat_transfer_coeff = 1000.0  # W/m²/K
    
    reactor.thermo.TPX = 1000.0, 101325.0, 'CH4:0.1, O2:0.2, N2:0.7'
    net = ct.ReactorNet([reactor])
    print(f"Before integration: T = {reactor.T:.1f} K")
    net.advance(0.01)
    print(f"After htc=1000 W/m²/K: T = {reactor.T:.1f} K")
    
    # Test 4: Check Cantera version and available methods
    print(f"\n=== CANTERA INFO ===")
    print(f"Cantera version: {ct.__version__}")
    print(f"Wall methods: {[m for m in dir(wall) if not m.startswith('_')]}")

if __name__ == '__main__':
    debug_wall_setup()
