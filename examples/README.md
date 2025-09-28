# Example Demonstrations

This directory contains example scripts demonstrating the HP-POX Benchmark Framework capabilities.

## Quick Start

Run all examples:

```bash
python examples/run_examples.py
```

## Individual Examples

### 1. HP-POX Validation

```bash
python main.py validate --case case_1 --mechanism data/gri30.yaml --output examples/results/validation
```

### 2. Plant Application

```bash
python main.py plant --geometry plant_a --mechanism data/gri30.yaml --sweep --output examples/results/plant
```

### 3. Mechanism Reduction

```bash
python main.py reduce --mechanism data/gri30.yaml --target-species 20 --cases case_1 --generations 10 --output examples/results/reduction
```

### 4. Plasma Analysis

```bash
python main.py plasma --type thermal --mechanism data/gri30.yaml --power-range 50,100,150 --output examples/results/plasma
```

## Expected Outputs

After running examples, you should see:

- `examples/results/validation/`: HP-POX validation results
- `examples/results/plant/`: Plant application results
- `examples/results/reduction/`: Mechanism reduction results
- `examples/results/plasma/`: Plasma analysis results

## Notes

- Examples use reduced parameters for faster execution
- For production runs, increase generations, population size, and discretization points
- Ensure mechanism files are available in `data/` directory
