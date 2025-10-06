# USC_MECH_II - USC Mech Version II

## Source Files
- **mechanism.yaml**: CHEMKIN mechanism file (reactions) - USC Mech ver 2.0, May 2007
- **thermdat.yaml**: CHEMKIN thermo file (NASA polynomials) 
- **transdat.yaml**: CHEMKIN transport file (LJ parameters)

## Conversion Process
These files were converted from CHEMKIN format to Cantera YAML using:

```bash
python -m cantera.ck2yaml \
  --input mechanism.yaml \
  --thermo thermdat.yaml \
  --transport transdat.yaml \
  --output uscmech2.yaml \
  --permissive
```

## Final Output
- **uscmech2.yaml**: Complete Cantera YAML file with mechanism + thermo + transport

## Species Count
- Total species: ~100+ (H2/CO/C1-C4 compounds)
- Elements: C, H, N, O, AR, HE, NE

## Transport Notes
- Uses mixture-averaged transport
- All species have transport parameters
- Some species may have estimated transport data

## Usage
```python
import cantera as ct
gas = ct.Solution('uscmech2.yaml')
```

## References
- USC Mech Version II. High-Temperature Combustion Reaction Model of H2/CO/C1-C4 Compounds
- Hai Wang, Xiaoqing You, Ameya V. Joshi, Scott G. Davis, Alexander Laskin, Fokion Egolfopoulos, and Chung K. Law
- http://ignis.usc.edu/USC_Mech_II.htm, May 2007
