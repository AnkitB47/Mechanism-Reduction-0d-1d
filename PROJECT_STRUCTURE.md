# HP-POX Benchmark Framework - Project Structure

## Overview

This is a complete, industry-grade 1-D reactor framework built on the Richter et al. (Fuel 2015) HP-POX benchmark dataset. The framework provides comprehensive validation, plant applications, mechanism reduction, and plasma surrogate capabilities.

## Directory Structure

```
Master_Thesis/
├── main.py                          # Main CLI interface
├── setup.py                         # Package setup
├── requirements.txt                 # Python dependencies
├── README.md                        # Main documentation
├── PROJECT_STRUCTURE.md            # This file
│
├── config/                          # Configuration files
│   ├── hp_pox_cases.yaml           # HP-POX benchmark data (Cases 1-4)
│   ├── plant_geometries.yaml       # Plant geometries A & B
│   ├── mechanisms.yaml             # Mechanism configurations
│   └── plasma_config.yaml          # Plasma surrogate settings
│
├── data/                            # Mechanism data
│   ├── README.md                   # Data directory instructions
│   ├── gri30.yaml                  # GRI-Mech 3.0 (placeholder)
│   ├── polimi.yaml                 # POLIMI C1-C3 (placeholder)
│   ├── aramco.yaml                 # AramcoMech 2.0 (placeholder)
│   └── usc_ii.yaml                 # USC-II (placeholder)
│
├── reactor/                         # Reactor models
│   ├── __init__.py                 # Reactor module init
│   ├── plug_flow.py                # 1-D plug-flow reactor (NEW)
│   ├── batch.py                    # Batch reactor (existing)
│   ├── piston.py                   # Piston reactor (existing)
│   ├── flame.py                    # Flame reactor (existing)
│   └── network.py                  # Reactor network (existing)
│
├── validation/                      # Validation framework
│   └── hp_pox_validator.py         # HP-POX validator (NEW)
│
├── applications/                    # Plant applications
│   └── plant_application.py        # Plant application framework (NEW)
│
├── surrogates/                      # Surrogate models
│   └── plasma_models.py            # Plasma surrogate models (NEW)
│
├── reduction/                       # Mechanism reduction
│   └── ga_gnn_pipeline.py          # GA-GNN reduction pipeline (NEW)
│
├── mechanism/                       # Mechanism handling
│   ├── __init__.py                 # Mechanism module init
│   ├── loader.py                   # Mechanism loader (existing)
│   ├── editor.py                   # Mechanism editor (existing)
│   └── mix.py                      # Mixture utilities (existing)
│
├── examples/                        # Example demonstrations
│   ├── README.md                   # Example documentation
│   └── run_examples.py             # Example runner script
│
├── results/                         # Output directory (created at runtime)
│   ├── validation/                 # HP-POX validation results
│   ├── plant_sweeps/               # Plant application results
│   ├── mechanism_reduction/        # Mechanism reduction results
│   └── plasma_sweeps/              # Plasma analysis results
│
└── [legacy files]                  # Existing thesis files (preserved)
    ├── gnn/                        # Graph neural network models
    ├── graph/                      # Graph construction
    ├── metaheuristics/             # Genetic algorithm
    ├── testing/                    # Test framework
    ├── tests/                      # Unit tests
    ├── visualizations/             # Plotting utilities
    └── ...                         # Other existing files
```

## Key Components

### 1. Core Reactor System (`reactor/plug_flow.py`)

- **1-D Plug-Flow Reactor**: Complete implementation with mass, energy, and momentum conservation
- **Heat Transfer Models**: Fixed wall temperature and effective heat transfer coefficient
- **Chemistry Integration**: Full Cantera integration with detailed kinetics
- **KPIs**: CH4/CO2 conversion, H2/CO ratio, ignition length, pressure drop

### 2. Validation Framework (`validation/hp_pox_validator.py`)

- **HP-POX Validation**: Complete validation against Richter et al. Cases 1-4
- **Temperature Validation**: Axial temperature profile comparison
- **Composition Validation**: Outlet syngas composition validation
- **Flame Validation**: Ignition length and flame characteristics
- **Automated Reporting**: CSV outputs and validation plots

### 3. Plant Application (`applications/plant_application.py`)

- **Plant Geometries**: Support for Plant A (450L) and Plant B (134L)
- **Operating Envelope Sweeps**: Comprehensive parameter sweeps
- **Plant Validation**: Validation against plant measurement data
- **Response Surfaces**: Visualization of operating parameter effects

### 4. Plasma Surrogates (`surrogates/plasma_models.py`)

- **Thermal Model**: Enthalpy rise and temperature enhancement
- **Radical Model**: Radical injection (H, O, OH) for enhanced kinetics
- **Parameter Sweeps**: Power range analysis
- **Performance Metrics**: Conversion improvement and efficiency gains

### 5. Mechanism Reduction (`reduction/ga_gnn_pipeline.py`)

- **Genetic Algorithm**: Species selection optimization
- **GNN Integration**: Graph neural network scoring (placeholder)
- **Validation Integration**: Automatic validation against test cases
- **Convergence Analysis**: GA convergence tracking and visualization

### 6. CLI Interface (`main.py`)

- **Comprehensive Commands**: validate, plant, reduce, plasma, pipeline
- **Flexible Configuration**: Command-line parameter control
- **Batch Processing**: Support for parameter sweeps
- **Error Handling**: Robust error handling and reporting

## Configuration System

### HP-POX Cases (`config/hp_pox_cases.yaml`)

Complete Richter et al. benchmark data:
- **Case 1**: Baseline (50 bar, 1200°C)
- **Case 2**: Pressure effect (60 bar, 1200°C)
- **Case 3**: Pressure effect (70 bar, 1200°C)
- **Case 4**: High temperature (50 bar, 1400°C)

Each case includes:
- Operating conditions and inlet streams
- Natural gas composition
- Expected outlet composition
- Flame characteristics
- Axial temperature measurements
- Wall temperature profiles

### Plant Geometries (`config/plant_geometries.yaml`)

- **Plant A**: Large scale (450L, 4.5m length)
- **Plant B**: Medium scale (134L, 3.2m length)

Each geometry includes:
- Dimensions and segments
- Operating envelopes
- Heat transfer configurations
- Natural gas compositions

## Usage Patterns

### 1. Research Validation

```bash
# Validate against HP-POX benchmark
python main.py validate --case case_1,case_4 --mechanism data/gri30.yaml
```

### 2. Industrial Application

```bash
# Plant application with operating envelope sweep
python main.py plant --geometry plant_a --mechanism data/gri30.yaml --sweep
```

### 3. Mechanism Development

```bash
# Reduce mechanism using GA-GNN
python main.py reduce --mechanism data/gri30.yaml --target-species 30 --cases case_1,case_4
```

### 4. Enhanced Combustion

```bash
# Test plasma surrogates
python main.py plasma --type thermal --mechanism data/gri30.yaml --power-range 50,100,150,200
```

### 5. Complete Pipeline

```bash
# Run full validation pipeline
python main.py pipeline --mechanism data/gri30.yaml --output results/full_validation
```

## Output Structure

### Validation Results

- `*_profiles.csv`: Axial profiles (T, P, species)
- `*_validation.yaml`: Validation summary
- `*_composition.png`: Composition comparison plots
- `*_profiles.png`: Axial profile plots

### Plant Results

- `*_sweep_results.csv`: Operating envelope results
- `*_response_surfaces.png`: Response surface plots
- `*_validation.csv`: Plant validation results

### Reduction Results

- `reduced_mechanism.yaml`: Reduced mechanism file
- `convergence.csv`: GA convergence data
- `best_individual.yaml`: Best individual parameters
- `reduction_analysis.png`: Convergence plots
- `species_importance.png`: Species importance analysis

### Plasma Results

- `*_sweep_results.csv`: Plasma parameter sweep results
- `*_analysis.png`: Plasma performance plots

## Integration Points

### 1. Cantera Integration

- Full Cantera solution integration
- Support for .yaml and .cti mechanism formats
- Mixture-averaged and multicomponent transport
- Detailed chemistry with configurable tolerances

### 2. Existing Codebase Integration

- Preserves existing reactor models (batch, piston, flame, network)
- Integrates with existing mechanism handling
- Compatible with existing GNN and GA implementations
- Maintains existing test framework

### 3. Configuration Integration

- YAML-based configuration system
- Command-line parameter override
- Flexible validation tolerances
- Extensible geometry and mechanism support

## Quality Assurance

### 1. Validation

- Comprehensive HP-POX benchmark validation
- Plant measurement validation
- Mechanism reduction validation
- Plasma surrogate validation

### 2. Error Handling

- Robust error handling throughout
- Graceful failure modes
- Detailed error reporting
- Validation tolerance checking

### 3. Documentation

- Comprehensive inline documentation
- Example scripts and demonstrations
- Configuration file documentation
- Usage pattern examples

## Future Extensions

### 1. Enhanced Chemistry

- More detailed chemistry models
- Advanced transport models
- Multi-phase considerations
- Catalytic effects

### 2. Advanced Surrogates

- Machine learning surrogates
- Neural network models
- Reduced-order models
- Uncertainty quantification

### 3. Industrial Integration

- Real-time plant data integration
- Advanced control strategies
- Optimization algorithms
- Digital twin capabilities

## Conclusion

This framework provides a complete, industry-grade solution for HP-POX validation and plant applications. It combines rigorous validation against the Richter et al. benchmark with practical industrial applications, mechanism reduction capabilities, and plasma surrogate models. The modular design allows for easy extension and customization while maintaining high code quality and comprehensive documentation.
