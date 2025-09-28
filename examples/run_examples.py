#!/usr/bin/env python3
"""Example scripts demonstrating the HP-POX Benchmark Framework.

This script runs example cases to demonstrate the framework capabilities.
"""

import subprocess
import sys
from pathlib import Path


def run_command(cmd, description):
    """Run a command and print results."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print('='*60)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("SUCCESS!")
        if result.stdout:
            print("Output:")
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("ERROR!")
        print(f"Return code: {e.returncode}")
        if e.stdout:
            print("Output:")
            print(e.stdout)
        if e.stderr:
            print("Error:")
            print(e.stderr)
        return False
    except FileNotFoundError:
        print("ERROR: Command not found. Make sure you're in the correct directory.")
        return False
    
    return True


def main():
    """Run example demonstrations."""
    
    print("HP-POX Benchmark Framework - Example Demonstrations")
    print("=" * 60)
    
    # Check if we're in the right directory
    if not Path("main.py").exists():
        print("ERROR: main.py not found. Please run from the project root directory.")
        sys.exit(1)
    
    # Example 1: HP-POX Validation
    print("\n1. HP-POX Validation Example")
    success = run_command([
        "python", "main.py", "validate",
        "--case", "case_1",
        "--mechanism", "data/gri30.yaml",
        "--output", "examples/results/validation"
    ], "Validating Case 1 against HP-POX benchmark")
    
    if not success:
        print("Skipping remaining examples due to error.")
        return
    
    # Example 2: Plant Application
    print("\n2. Plant Application Example")
    run_command([
        "python", "main.py", "plant",
        "--geometry", "plant_a",
        "--mechanism", "data/gri30.yaml",
        "--sweep",
        "--output", "examples/results/plant"
    ], "Running plant application sweep")
    
    # Example 3: Mechanism Reduction (small example)
    print("\n3. Mechanism Reduction Example")
    run_command([
        "python", "main.py", "reduce",
        "--mechanism", "data/gri30.yaml",
        "--target-species", "20",
        "--cases", "case_1",
        "--generations", "10",
        "--population-size", "20",
        "--output", "examples/results/reduction"
    ], "Running mechanism reduction (small example)")
    
    # Example 4: Plasma Analysis
    print("\n4. Plasma Analysis Example")
    run_command([
        "python", "main.py", "plasma",
        "--type", "thermal",
        "--mechanism", "data/gri30.yaml",
        "--power-range", "50,100,150",
        "--output", "examples/results/plasma"
    ], "Running plasma surrogate analysis")
    
    print("\n" + "="*60)
    print("Example demonstrations completed!")
    print("Check the examples/results/ directory for output files.")
    print("="*60)


if __name__ == "__main__":
    main()
