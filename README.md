# HP-POX Benchmark Framework

Industry-grade 1-D reactor framework for HP-POX validation and plant applications. Built on the Richter et al. (Fuel 2015) benchmark dataset with comprehensive validation, mechanism reduction, and plasma surrogate capabilities.

## Features

- **HP-POX Validation**: Complete validation against Richter et al. Cases 1-4
- **Plant Applications**: Support for industrial geometries A & B with operating envelope sweeps
- **Mechanism Reduction**: GA-GNN pipeline for automated mechanism reduction
- **Plasma Surrogates**: Thermal and radical plasma models for enhanced combustion
- **Comprehensive CLI**: Full command-line interface for all operations

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

### 1. HP-POX Validation

Validate against the HP-POX benchmark:

```bash
python main.py validate --case case_1,case_4 --mechanism data/gri30.yaml
```

### 2. Plant Application

Run operating envelope sweep for plant geometry:

```bash
python main.py plant --geometry plant_a --mechanism data/gri30.yaml --sweep
```

### 3. Mechanism Reduction

Reduce mechanism using GA-GNN pipeline:

```bash
python main.py reduce --mechanism data/gri30.yaml --target-species 30 --cases case_1,case_4
```

### 4. Plasma Analysis

Test plasma surrogate models:

```bash
python main.py plasma --type thermal --mechanism data/gri30.yaml --power-range 50,100,150,200
```

### 5. Complete Pipeline

Run full validation pipeline:

```bash
python main.py pipeline --mechanism data/gri30.yaml --output results/full_validation
```

## Configuration

### HP-POX Cases

Configuration files are located in `config/`:

- `hp_pox_cases.yaml`: Complete HP-POX benchmark data (Cases 1-4)
- `plant_geometries.yaml`: Plant geometries A & B
- `mechanisms.yaml`: Supported mechanisms and settings

### Mechanism Files

Place mechanism files in `data/`:

- `gri30.yaml`: GRI-Mech 3.0 (recommended)
- `polimi.yaml`: POLIMI C1-C3
- `aramco.yaml`: AramcoMech 2.0
- `usc_ii.yaml`: USC-II

## Architecture

```
├── main.py                     # CLI interface
├── config/                     # Configuration files
│   ├── hp_pox_cases.yaml      # HP-POX benchmark data
│   ├── plant_geometries.yaml  # Plant geometries
│   └── mechanisms.yaml        # Mechanism settings
├── reactor/                    # Reactor models
│   ├── plug_flow.py           # 1-D plug-flow reactor
│   └── batch.py               # Batch reactor (legacy)
├── validation/                 # Validation framework
│   └── hp_pox_validator.py    # HP-POX validator
├── applications/               # Plant applications
│   └── plant_application.py   # Plant application framework
├── surrogates/                 # Surrogate models
│   └── plasma_models.py       # Plasma surrogate models
├── reduction/                  # Mechanism reduction
│   └── ga_gnn_pipeline.py     # GA-GNN reduction pipeline
├── mechanism/                  # Mechanism handling
│   ├── loader.py              # Mechanism loader
│   └── editor.py              # Mechanism editor
└── results/                    # Output directory
```

## Key Performance Indicators

The framework calculates and validates:

- **CH4 Conversion**: Methane conversion percentage
- **CO2 Conversion**: Carbon dioxide conversion percentage
- **H2/CO Ratio**: Hydrogen to carbon monoxide ratio
- **Ignition Length**: Flame stabilization length
- **Pressure Drop**: Axial pressure drop
- **Temperature Profiles**: Axial temperature distribution
- **Species Profiles**: Axial species concentration profiles

## Validation Tolerances

Default validation tolerances:

- Temperature: ±50°C
- Composition: ±5%
- Flame length: ±20 mm
- Pressure drop: ±10%

## Output Files

Results are saved in the specified output directory:

- `*_profiles.csv`: Axial profiles (T, P, species)
- `*_validation.yaml`: Validation results
- `*_composition.png`: Composition comparison plots
- `*_profiles.png`: Axial profile plots
- `convergence.csv`: GA convergence data
- `reduced_mechanism.yaml`: Reduced mechanism file

## Examples

### Example 1: Basic Validation

```bash
# Validate Case 1 against GRI-3.0
python main.py validate --case case_1 --mechanism data/gri30.yaml --output results/case1
```

### Example 2: Plant Sweep

```bash
# Run operating envelope sweep for Plant A
python main.py plant --geometry plant_a --mechanism data/gri30.yaml --sweep --output results/plant_a_sweep
```

### Example 3: Mechanism Reduction

```bash
# Reduce GRI-3.0 to 30 species
python main.py reduce --mechanism data/gri30.yaml --target-species 30 --cases case_1,case_4 --generations 50
```

### Example 4: Plasma Analysis

```bash
# Test thermal plasma surrogate
python main.py plasma --type thermal --mechanism data/gri30.yaml --power-range 50,100,150,200,250,300
```

## Troubleshooting

### Common Issues

1. **Mechanism file not found**: Ensure mechanism files are in `data/` directory
2. **Validation failures**: Check mechanism compatibility with HP-POX conditions
3. **Memory issues**: Reduce `--n-points` for large mechanisms
4. **Convergence issues**: Increase GA generations or population size

### Debug Mode

Use `--verbose` flag for detailed output:

```bash
python main.py validate --case case_1 --mechanism data/gri30.yaml --verbose
```

## References

- Richter, H., et al. "High-pressure partial oxidation of natural gas: A validation dataset for CFD and 1-D model development." Fuel 2015.
- GRI-Mech 3.0: http://www.me.berkeley.edu/gri_mech/
- Cantera: https://cantera.org/

## License

This project is part of a Master's thesis research. Please cite appropriately if used in academic work.

## Contact

For questions or issues, please refer to the thesis documentation or contact the author.
