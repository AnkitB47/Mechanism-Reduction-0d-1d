"""
Case constants and configuration for HP-POX model.
Defines Richter benchmark cases and geometry parameters.
"""

import json
import math
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass


# Typed configuration schema
@dataclass
class NozzleConfig:
    """Nozzle configuration for jet velocity calculations."""
    area_m2: Optional[float] = None  # m²
    diameter_m: Optional[float] = None  # m
    C_d: float = 0.8  # Discharge coefficient (-)

@dataclass
class MixingConfig:
    """Mixing duct configuration."""
    L_ign_m: float = 0.10  # Mixing duct length (m)
    q_pre_W: float = 0.0   # Pre-mix wall loss (W)

@dataclass
class KernelConfig:
    """Ignition kernel configuration."""
    V_kernel_frac_psr: float = 0.01  # Fraction of PSR volume for kernel
    tau_kernel_s: float = 0.015      # Kernel residence time (s)
    y_feed_to_flame: float = 0.4     # Fraction of fuel+N2 that forms kernel

@dataclass
class StreamConfig:
    """Inlet stream configuration."""
    mass_kg_per_s: float  # kg/s
    T_K: float            # K
    composition_X: Dict[str, float]  # Mole fractions (-)

@dataclass
class CaseConfig:
    """Complete case configuration."""
    # Basic properties
    case_id: str
    pressure_pa: float  # Pa

    # Inlet streams (by name)
    primary_steam: StreamConfig
    secondary_steam: StreamConfig
    oxygen: StreamConfig
    nitrogen: StreamConfig
    natural_gas: StreamConfig

    # Losses
    burner_heat_loss_W: float  # W
    wall_heat_loss_W: float    # W

    # Geometry
    psr_volume_m3: float       # m³
    pfr_length_m: float        # m
    pfr_diameter_m: float      # m
    pfr_area_m2: float         # m²

    # New velocity and ignition parameters
    nozzles: Dict[str, NozzleConfig]
    mixing: MixingConfig
    kernel: KernelConfig

    # Validation targets (dry basis unless noted)
    target_T_K: float
    target_H2_pct: float
    target_N2_pct: float
    target_CO_pct: float
    target_CO2_pct: float
    target_CH4_pct: float
    target_H2O_wet_pct: float


# Baseline diagnostics loader -------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
BASELINE_RESULTS_ROOT = REPO_ROOT / 'hp_pox_results_full'


def _diagnostics_path_candidates(case_id: str) -> List[Path]:
    """Enumerate plausible diagnostics.json locations for a given case."""
    return [
        BASELINE_RESULTS_ROOT / case_id / case_id / 'diagnostics.json',
        BASELINE_RESULTS_ROOT / case_id / 'diagnostics.json',
    ]


def _load_target_temperature(case_id: str, default_K: float) -> float:
    """
    Load baseline outlet temperature (K) from diagnostics if available.

    Falls back to provided default when diagnostics are missing or invalid.
    """
    for candidate in _diagnostics_path_candidates(case_id):
        try:
            data = json.loads(candidate.read_text())
        except FileNotFoundError:
            continue
        except json.JSONDecodeError:
            continue

        value = data.get('T_out')
        try:
            value_f = float(value)
        except (TypeError, ValueError):
            continue

        if math.isfinite(value_f) and value_f > 0.0:
            return value_f
    return default_K


class CaseConstants:
    """Case constants for HP-POX model."""
    
    # Case A (Richter geometry)
    CASE_A_LENGTH = 2.32  # m
    CASE_A_DIAMETER = 0.49  # m
    CASE_A_AREA = 3.14159 * (0.49 / 2)**2  # m²
    
    # Case B (134L geometry)
    CASE_B_LENGTH = 1.00  # m
    CASE_B_DIAMETER = 0.387  # m (exact from REF)
    CASE_B_AREA = 3.14159 * (0.387 / 2)**2  # m²
    
    # Burner heat loss (kW) - Fixed values from Richter benchmark
    BURNER_HEAT_LOSS = {
        'A_case1_richter': 34.29,    # kW
        'A_case4_richter': 36.36,    # kW
        'B_case1_134L': 34.29,       # kW (mirror A_case1)
        'B_case4_134L': 36.36,       # kW (mirror A_case4)
    }
    
    # Wall heat loss per length (W/m)
    # Case A1: 18.50 kW / 2.32 m = 7974.1 W/m
    # Case A4: 22.92 kW / 2.32 m = 9876.7 W/m
    # Case B1: 18.50 kW / 1.00 m = 18500.0 W/m
    # Case B4: 22.92 kW / 1.00 m = 22920.0 W/m
    
    WALL_HEAT_LOSS = {
        'A_case1_richter': 18500.0 * (1.0/2.32),   # 18.50 kW / 2.32 m ≈ 7,974.1 W/m
        'A_case4_richter': 22920.0 * (1.0/2.32),   # 22.92 kW / 2.32 m ≈ 9,887.9 W/m
        'B_case1_134L': 18500.0,                   # W/m (18.50 kW over 1.00 m)
        'B_case4_134L': 22920.0,                   # W/m (22.92 kW over 1.00 m)
    }
    
    # Target temperatures (K) sourced from latest baseline diagnostics when available
    TARGET_TEMPERATURES = {
        'A_case1_richter': _load_target_temperature('A_case1_richter', 1201 + 273.15),
        'A_case4_richter': _load_target_temperature('A_case4_richter', 1401 + 273.15),
        'B_case1_134L': _load_target_temperature('B_case1_134L', 1201 + 273.15),
        'B_case4_134L': _load_target_temperature('B_case4_134L', 1401 + 273.15),
    }
    
    # Target KPIs
    TARGET_KPIS = {
        'h2_co_ratio': 2.0,      # Target H2/CO ratio
        'ch4_conversion': 0.95,  # 95% CH4 conversion
        'co2_conversion': 0.05,  # 5% CO2 conversion
    }

    # Pressure (Pa)
    PRESSURE = 51.013e5  # 51.013 bar abs

    # Mass flow rate (kg/s) - from Richter benchmark
    MASS_FLOW_RATE = 0.164  # kg/s

    # Legacy helper methods
    @classmethod
    def get_geometry(cls, case_id: str) -> Dict[str, float]:
        """Get geometry parameters for a case."""
        key = case_id.lower()
        if 'case_a' in key or 'richter' in key:
            return {
                'length_m': cls.CASE_A_LENGTH,
                'diameter_m': cls.CASE_A_DIAMETER,
                'area_m2': cls.CASE_A_AREA,
            }
        if 'case_b' in key or '134l' in key:
            return {
                'length_m': cls.CASE_B_LENGTH,
                'diameter_m': cls.CASE_B_DIAMETER,
                'area_m2': cls.CASE_B_AREA,
            }
        raise ValueError(f"Unknown case ID: {case_id}")

    @classmethod
    def get_burner_heat_loss(cls, case_id: str) -> float:
        """Get burner heat loss (W)."""
        if case_id in cls.BURNER_HEAT_LOSS:
            return cls.BURNER_HEAT_LOSS[case_id] * 1000.0
        raise ValueError(f"Unknown case ID: {case_id}")

    @classmethod
    def get_wall_heat_loss(cls, case_id: str) -> float:
        """Get wall heat loss per length (W/m)."""
        if case_id in cls.WALL_HEAT_LOSS:
            return cls.WALL_HEAT_LOSS[case_id]
        raise ValueError(f"Unknown case ID: {case_id}")

    @classmethod
    def get_target_temperature(cls, case_id: str) -> float:
        """Get target outlet temperature (K)."""
        if case_id in cls.TARGET_TEMPERATURES:
            return cls.TARGET_TEMPERATURES[case_id]
        raise ValueError(f"Unknown case ID: {case_id}")

    @classmethod
    def get_case_config(cls, case_id: str) -> Dict[str, Any]:
        """Get complete case configuration bundle."""
        return {
            'case_id': case_id,
            'geometry': cls.get_geometry(case_id),
            'burner_heat_loss': cls.get_burner_heat_loss(case_id),
            'wall_heat_loss': cls.get_wall_heat_loss(case_id),
            'target_temperature': cls.get_target_temperature(case_id),
            'pressure': cls.PRESSURE,
            'mass_flow_rate': cls.MASS_FLOW_RATE,
            'target_kpis': cls.TARGET_KPIS,
        }

# Richter benchmark case configurations
def _create_richter_case_1() -> CaseConfig:
    """Create Case 1 configuration."""
    # PSR volume: 0.016136 m³
    psr_volume = 0.016136

    # PFR geometry: L = 2.32 m, D = 0.49 m
    pfr_length = 2.32
    pfr_diameter = 0.49
    pfr_area = 3.14159 * (pfr_diameter / 2)**2

    # Pressure: 51 bar abs
    p_abs = 5.1e6

    wall_loss_total_W = 18_500.0
    wall_loss_per_m = wall_loss_total_W / pfr_length

    return CaseConfig(
        case_id='A_case1_richter',
        pressure_pa=p_abs,
        primary_steam=StreamConfig(
            mass_kg_per_s=0.02167,
            T_K=562.75,
            composition_X={'H2O': 1.0}
        ),
        secondary_steam=StreamConfig(
            mass_kg_per_s=0.00643,
            T_K=513.35,
            composition_X={'H2O': 1.0}
        ),
        oxygen=StreamConfig(
            mass_kg_per_s=0.07371,
            T_K=513.35,
            composition_X={'O2': 1.0}
        ),
        nitrogen=StreamConfig(
            mass_kg_per_s=0.00131,
            T_K=298.15,
            composition_X={'N2': 1.0}
        ),
        natural_gas=StreamConfig(
            mass_kg_per_s=0.06224,
            T_K=632.25,
            composition_X={
                'CH4': 0.9597,
                'CO2': 0.0034,
                'C2H6': 0.0226,
                'C3H8': 0.0046,
                'N2': 0.0086
            }
        ),
        burner_heat_loss_W=34290.0,  # W
        wall_heat_loss_W=wall_loss_per_m,  # W/m
        psr_volume_m3=psr_volume,
        pfr_length_m=pfr_length,
        pfr_diameter_m=pfr_diameter,
        pfr_area_m2=pfr_area,
        nozzles={
            'fuel_premix': NozzleConfig(area_m2=None, diameter_m=None, C_d=0.8),
            'oxid_premix': NozzleConfig(area_m2=None, diameter_m=None, C_d=0.8),
            'steam_primary': NozzleConfig(area_m2=None, diameter_m=None, C_d=0.8)
        },
        mixing=MixingConfig(L_ign_m=0.10, q_pre_W=0.0),
        kernel=KernelConfig(V_kernel_frac_psr=0.01, tau_kernel_s=0.015, y_feed_to_flame=0.4),
        target_T_K=CaseConstants.TARGET_TEMPERATURES['A_case1_richter'],
        target_H2_pct=48.27,
        target_N2_pct=0.65,
        target_CO_pct=23.79,
        target_CO2_pct=4.19,
        target_CH4_pct=3.76,
        target_H2O_wet_pct=19.33
    )


def _create_richter_case_4() -> CaseConfig:
    """Create Case 4 configuration."""
    # PSR volume: 0.016136 m³
    psr_volume = 0.016136

    # PFR geometry: L = 2.32 m, D = 0.49 m
    pfr_length = 2.32
    pfr_diameter = 0.49
    pfr_area = 3.14159 * (pfr_diameter / 2)**2

    # Pressure: 51 bar abs
    p_abs = 5.1e6

    wall_loss_total_W = 22_920.0
    wall_loss_per_m = wall_loss_total_W / pfr_length

    return CaseConfig(
        case_id='A_case4_richter',
        pressure_pa=p_abs,
        primary_steam=StreamConfig(
            mass_kg_per_s=0.01789,
            T_K=550.65,
            composition_X={'H2O': 1.0}
        ),
        secondary_steam=StreamConfig(
            mass_kg_per_s=0.00633,
            T_K=513.05,
            composition_X={'H2O': 1.0}
        ),
        oxygen=StreamConfig(
            mass_kg_per_s=0.07785,
            T_K=513.05,
            composition_X={'O2': 1.0}
        ),
        nitrogen=StreamConfig(
            mass_kg_per_s=0.00133,
            T_K=298.15,
            composition_X={'N2': 1.0}
        ),
        natural_gas=StreamConfig(
            mass_kg_per_s=0.05442,
            T_K=628.35,
            composition_X={
                'CH4': 0.9620,
                'CO2': 0.0026,
                'C2H6': 0.0209,
                'C3H8': 0.0047,
                'N2': 0.0087
            }
        ),
        burner_heat_loss_W=36360.0,  # W
        wall_heat_loss_W=wall_loss_per_m,  # W/m
        psr_volume_m3=psr_volume,
        pfr_length_m=pfr_length,
        pfr_diameter_m=pfr_diameter,
        pfr_area_m2=pfr_area,
        nozzles={
            'fuel_premix': NozzleConfig(area_m2=None, diameter_m=None, C_d=0.8),
            'oxid_premix': NozzleConfig(area_m2=None, diameter_m=None, C_d=0.8),
            'steam_primary': NozzleConfig(area_m2=None, diameter_m=None, C_d=0.8)
        },
        mixing=MixingConfig(L_ign_m=0.10, q_pre_W=0.0),
        kernel=KernelConfig(V_kernel_frac_psr=0.01, tau_kernel_s=0.015, y_feed_to_flame=0.4),
        target_T_K=CaseConstants.TARGET_TEMPERATURES['A_case4_richter'],
        target_H2_pct=48.06,
        target_N2_pct=0.67,
        target_CO_pct=25.61,
        target_CO2_pct=3.89,
        target_CH4_pct=0.06,
        target_H2O_wet_pct=21.71
    )


# Predefined cases
richter_case_1 = _create_richter_case_1()
richter_case_4 = _create_richter_case_4()




def get_available_cases() -> list:
    """Get list of available cases.
    
    Returns:
        List of case identifiers
    """
    return list(CaseConstants.WALL_HEAT_LOSS.keys())


def print_case_summary(case_id: str):
    """Print case summary.
    
    Args:
        case_id: Case identifier
    """
    config = CaseConstants.get_case_config(case_id)
    
    print(f"\nCase: {case_id}")
    print(f"  Geometry: L={config['geometry']['length_m']:.2f}m, D={config['geometry']['diameter_m']:.2f}m")
    print(f"  Burner heat loss: {config['burner_heat_loss']:.2f} kW")
    print(f"  Wall heat loss: {config['wall_heat_loss']:.1f} W/m")
    print(f"  Target temperature: {config['target_temperature']-273.15:.1f}°C")
    print(f"  Pressure: {config['pressure']/1e5:.1f} bar")
    print(f"  Mass flow rate: {config['mass_flow_rate']:.3f} kg/s")


if __name__ == "__main__":
    # Print all available cases
    print("Available cases:")
    for case_id in get_available_cases():
        print_case_summary(case_id)


