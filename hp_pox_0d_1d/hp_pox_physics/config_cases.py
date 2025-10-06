"""
Case constants and configuration for HP-POX model.
Defines Richter benchmark cases and geometry parameters.
"""

from typing import Dict, Any


class CaseConstants:
    """Case constants for HP-POX model."""
    
    # Case A (Richter geometry)
    CASE_A_LENGTH = 2.32  # m
    CASE_A_DIAMETER = 0.49  # m
    CASE_A_AREA = 3.14159 * (0.49 / 2)**2  # m²
    
    # Case B (134L geometry)
    CASE_B_LENGTH = 1.00  # m
    CASE_B_DIAMETER = 0.39  # m
    CASE_B_AREA = 3.14159 * (0.39 / 2)**2  # m²
    
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
        'A_case1_richter': 7974.1,   # W/m
        'A_case4_richter': 9876.7,   # W/m
        'B_case1_134L': 18500.0,     # W/m
        'B_case4_134L': 22920.0,     # W/m
    }
    
    # Target temperatures (K)
    TARGET_TEMPERATURES = {
        'A_case1_richter': 1201 + 273.15,  # 1474 K
        'A_case4_richter': 1401 + 273.15,  # 1674 K
        'B_case1_134L': 1201 + 273.15,     # 1474 K
        'B_case4_134L': 1401 + 273.15,     # 1674 K
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
    
    @classmethod
    def get_geometry(cls, case_id: str) -> Dict[str, float]:
        """Get geometry parameters for case.
        
        Args:
            case_id: Case identifier
            
        Returns:
            Geometry parameters
        """
        if 'case_a' in case_id.lower() or 'richter' in case_id.lower():
            return {
                'length_m': cls.CASE_A_LENGTH,
                'diameter_m': cls.CASE_A_DIAMETER,
                'area_m2': cls.CASE_A_AREA
            }
        elif 'case_b' in case_id.lower() or '134l' in case_id.lower():
            return {
                'length_m': cls.CASE_B_LENGTH,
                'diameter_m': cls.CASE_B_DIAMETER,
                'area_m2': cls.CASE_B_AREA
            }
        else:
            raise ValueError(f"Unknown case ID: {case_id}")
    
    @classmethod
    def get_burner_heat_loss(cls, case_id: str) -> float:
        """Get burner heat loss for case.
        
        Args:
            case_id: Case identifier
            
        Returns:
            Burner heat loss (kW)
        """
        if case_id in cls.BURNER_HEAT_LOSS:
            return cls.BURNER_HEAT_LOSS[case_id]
        else:
            raise ValueError(f"Unknown case ID: {case_id}")
    
    @classmethod
    def get_wall_heat_loss(cls, case_id: str) -> float:
        """Get wall heat loss per length for case.
        
        Args:
            case_id: Case identifier
            
        Returns:
            Wall heat loss per length (W/m)
        """
        if case_id in cls.WALL_HEAT_LOSS:
            return cls.WALL_HEAT_LOSS[case_id]
        else:
            raise ValueError(f"Unknown case ID: {case_id}")
    
    @classmethod
    def get_target_temperature(cls, case_id: str) -> float:
        """Get target temperature for case.
        
        Args:
            case_id: Case identifier
            
        Returns:
            Target temperature (K)
        """
        if case_id in cls.TARGET_TEMPERATURES:
            return cls.TARGET_TEMPERATURES[case_id]
        else:
            raise ValueError(f"Unknown case ID: {case_id}")
    
    @classmethod
    def get_case_config(cls, case_id: str) -> Dict[str, Any]:
        """Get complete case configuration.
        
        Args:
            case_id: Case identifier
            
        Returns:
            Complete case configuration
        """
        return {
            'case_id': case_id,
            'geometry': cls.get_geometry(case_id),
            'burner_heat_loss': cls.get_burner_heat_loss(case_id),
            'wall_heat_loss': cls.get_wall_heat_loss(case_id),
            'target_temperature': cls.get_target_temperature(case_id),
            'pressure': cls.PRESSURE,
            'mass_flow_rate': cls.MASS_FLOW_RATE,
            'target_kpis': cls.TARGET_KPIS
        }


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
