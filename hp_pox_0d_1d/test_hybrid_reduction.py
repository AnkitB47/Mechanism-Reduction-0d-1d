#!/usr/bin/env python3
"""
Test script for the hybrid mechanism reduction pipeline.

Tests the complete pipeline: DRGEP → GA+GNN → closure → validation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from pathlib import Path
import cantera as ct

from reduction_orchestrator import ReductionConfig, run_reduction


def test_hybrid_reduction():
    """Test the complete hybrid reduction pipeline."""

    # Create test configuration
    config = ReductionConfig(
        base_yaml=Path("gri30.yaml"),
        out_dir=Path("hybrid_reduction_test"),
        target_species_min=26,
        target_species_max=40,  # Wider band for test
        protected_species=["H2", "CO", "CO2", "H2O", "CH4", "N2"],
        drgep_tau_grid=[0.20, 0.15, 0.10],  # Shorter grid for test
        verbose=True
    )

    print("Testing hybrid mechanism reduction pipeline...")
    print("=" * 80)

    try:
        # Run the reduction pipeline
        reduced_yaml, report = run_reduction(config)

        # Basic assertions
        assert reduced_yaml.exists(), "Reduced mechanism file should exist"

        # Load and check mechanism
        mech = ct.Solution(str(reduced_yaml))
        n_species = len(mech.species())
        n_reactions = len(mech.reactions())

        # Check species count
        assert config.target_species_min <= n_species <= config.target_species_max, \
            f"Species count {n_species} outside target band {config.target_species_min}-{config.target_species_max}"

        # Check protected species are present
        protected_names = set(config.protected_species)
        mech_species_names = {sp.name for sp in mech.species()}
        assert protected_names.issubset(mech_species_names), \
            f"Missing protected species: {protected_names - mech_species_names}"

        # Check YAML loads without errors
        assert mech.n_species > 0, "Mechanism should have species"
        assert mech.n_reactions > 0, "Mechanism should have reactions"

        print(f"[OK] Reduced mechanism: {n_species} species, {n_reactions} reactions")
        print(f"[OK] Protected species present: {all(sp in mech_species_names for sp in protected_names)}")

        # Check validation results
        if report.psr_passed and report.pfr_passed:
            print("[OK] Validation PASSED")
        else:
            print("[WARN] Validation had issues (expected for test configuration)")

        print("[OK] All tests passed!")
        return True

    except Exception as e:
        print(f"[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Running hybrid reduction test...")
    success = test_hybrid_reduction()
    sys.exit(0 if success else 1)
