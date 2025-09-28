# Mechanism Data Directory

Place mechanism files in this directory.

## Required Files

- `gri30.yaml`: GRI-Mech 3.0 mechanism (53 species, 325 reactions)
- `polimi.yaml`: POLIMI C1-C3 mechanism (78 species, 1200 reactions)
- `aramco.yaml`: AramcoMech 2.0 mechanism (295 species, 1840 reactions)
- `usc_ii.yaml`: USC-II mechanism (111 species, 784 reactions)

## Download Sources

- GRI-Mech 3.0: http://www.me.berkeley.edu/gri_mech/
- POLIMI: https://www.chem.polimi.it/combustion/mechanism.php
- AramcoMech: https://www.chem.ucla.edu/~combustion/cermech/
- USC-II: https://ignis.usc.edu/Mechanisms/

## File Format

All mechanism files should be in Cantera YAML format (.yaml) or CTI format (.cti).

## Validation

The framework will automatically validate mechanism compatibility with HP-POX conditions.
