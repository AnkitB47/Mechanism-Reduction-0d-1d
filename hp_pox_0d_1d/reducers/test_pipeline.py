#!/usr/bin/env python3
"""
Test script for the advanced mechanism reduction pipeline.

Demonstrates the full pipeline: data collection → DRGEP → DRGASA → CSP/QSS → closure → validation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pathlib import Path
from . import ReductionConfig
from .pipeline import reduce_mechanism_grimech


def test_reduction_pipeline():
    """Test the complete reduction pipeline."""

    # Create a test configuration
    config = ReductionConfig(
        base_mechanism=Path("gri30.yaml"),
        output_dir=Path("reduction_test_output"),
        verbose=True,
        # Reduce sampling for quick test
        psr_conditions=[
            {"T": 900, "phi": 2.5, "P": 1.0, "tau": 0.010},
            {"T": 1000, "phi": 2.5, "P": 1.0, "tau": 0.010},
        ],
        pfr_cases=["A_case1_richter"],
        sample_points=[0.0, 0.5, 1.0],
        tau_s_grid=[0.01, 0.05, 0.10]  # Reduced grid for quick test
    )

    print("Testing mechanism reduction pipeline...")
    print("=" * 80)

    try:
        # Run the reduction pipeline
        result = reduce_mechanism_grimech(config)

        print("\n" + "=" * 80)
        print("REDUCTION PIPELINE TEST COMPLETE")
        print("=" * 80)
        print(f"Original: {result.original_species} species, {result.original_reactions} reactions")
        print(f"Reduced: {result.reduced_species} species, {result.reduced_reactions} reactions")
        print(f"Reduction: {100*(1-result.reduced_species/result.original_species):.1f}% species")
        print(f"Output files in: {config.output_dir}")

        # Check if validation passed
        if result.psr_convergence:
            passed_cases = sum(result.psr_convergence.values())
            total_cases = len(result.psr_convergence)
            print(f"Validation: {passed_cases}/{total_cases} cases passed")

        return True

    except Exception as e:
        print(f"Pipeline test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_reduction_pipeline()
    sys.exit(0 if success else 1)
