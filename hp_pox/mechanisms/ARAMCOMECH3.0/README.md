# ARAMCOMECH3.0 - AramcoMech 3.0 Chemical Mechanism

## Source Files
- **mechanism.yaml**: CHEMKIN mechanism file (reactions) - originally from 17_04_irreversable_changes_CLEAN_RED_LTC.MECH
- **thermodynamics.yaml**: CHEMKIN thermo file (NASA polynomials) - originally from 16_26.THERM  
- **transport.yaml**: CHEMKIN transport file (LJ parameters) - originally from 16_26.TRAN

## Conversion Process
These files were converted from CHEMKIN format to Cantera YAML using:

```bash
python -m cantera.ck2yaml \
  --input mechanism.yaml \
  --thermo thermodynamics.yaml \
  --transport transport.yaml \
  --output aramco3.yaml \
  --permissive
```

## Final Output
- **aramco3.yaml**: Complete Cantera YAML file with mechanism + thermo + transport

## Species Count
- Total species: ~200+ (C1-C4 hydrocarbons, H2/CO chemistry)
- Elements: C, H, N, O, AR, HE

## Transport Notes
- Uses mixture-averaged transport
- All species have transport parameters
- Some species may have estimated transport data

## Usage
```python
import cantera as ct
gas = ct.Solution('aramco3.yaml')
```

## References
- AramcoMech 3.0: A detailed chemical kinetic model for natural gas combustion
- Original files cleaned and sorted on 07/07/2016
