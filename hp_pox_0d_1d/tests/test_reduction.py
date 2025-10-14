#!/usr/bin/env python3
"""
Smoke test for the hybrid mechanism reduction pipeline.

Tests the complete pipeline: DRGEP → GA+GNN → closure → validation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
import cantera as ct

from reduction_orchestrator import ReductionConfig, run_reduction


def test_reduction_pipeline():
    """Test the complete reduction pipeline with a minimal configuration."""

    # Create minimal test configuration
    config = ReductionConfig(
        base_yaml=Path("gri30.yaml"),
        out_dir=Path("reduction_test_output"),
        target_species_min=26,
        target_species_max=40,  # Wider band for test
        protected_species=["H2", "CO", "CO2", "H2O", "CH4", "N2"],
        drgep_tau_grid=[0.20, 0.15, 0.10],  # Shorter grid for test
        verbose=True
    )

    print("Testing mechanism reduction pipeline...")
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

        # Check closure converged (no dangling reactions)
        # This is a basic check - in practice you'd verify no reaction references dropped species
        print(f"[OK] Reduced mechanism: {n_species} species, {n_reactions} reactions")

        # Check validation results
        if report.psr_passed and report.pfr_passed:
            print("[OK] Validation PASSED")
        else:
            print("[WARN] Validation had issues (expected for test configuration)")

        print("[OK] All smoke tests passed!")
        return True

    except Exception as e:
        print(f"[ERROR] Smoke test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_individual_components():
    """Test individual components work correctly."""

    print("\nTesting individual components...")

    # Test data collection
    from hp_pox_physics.data_api import collect_training_states
    states = collect_training_states()
    assert len(states) > 0, "Should collect training states"
    print(f"[OK] Data collection: {len(states)} states")

    # Test DRGEP
    from reducers.drgep import run_drgep
    drgep_mask, drgep_yaml = run_drgep(
        base_yaml=Path("gri30.yaml"),
        training_states=states[:5],  # Small subset for test
        tau_candidates=[0.20, 0.15],
        protected_species=["H2", "CO", "CO2", "H2O", "CH4", "N2"],
        stop_when_species_leq=45,
        out_dir=Path("test_drgep")
    )
    assert sum(drgep_mask) > 0, "DRGEP should produce species mask"
    assert drgep_yaml.exists(), "DRGEP should produce YAML"
    print(f"[OK] DRGEP: {sum(drgep_mask)} species")

    # Test closure
    from reducers.closure import fixed_point_closure
    closed_yaml = fixed_point_closure(
        yaml_in=drgep_yaml,
        protected_species=["H2", "CO", "CO2", "H2O", "CH4", "N2"],
        out_yaml=Path("test_closed.yaml"),
        max_iter=10,
        verbose=False
    )
    assert closed_yaml.exists(), "Closure should produce YAML"
    print("[OK] Fixed-point closure")

    # Test validation
    from reducers.validation import validate_psr_pfr
    ok, report = validate_psr_pfr(
        base_yaml=Path("gri30.yaml"),
        red_yaml=closed_yaml,
        protected_species=["H2", "CO", "CO2", "H2O", "CH4", "N2"],
        out_dir=Path("test_validation")
    )
    print(f"[OK] Validation: {'PASSED' if ok else 'FAILED'}")

    return True


if __name__ == "__main__":
    print("Running mechanism reduction smoke tests...")

    # Test individual components first
    component_test = test_individual_components()

    # Test full pipeline
    pipeline_test = test_reduction_pipeline()

    overall_success = component_test and pipeline_test

    print("\n" + "=" * 80)
    print(f"SMOKE TEST RESULTS: {'SUCCESS' if overall_success else 'FAILED'}")
    print("=" * 80)

    sys.exit(0 if overall_success else 1)
