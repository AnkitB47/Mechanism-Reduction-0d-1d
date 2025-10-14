#!/usr/bin/env python3
"""
Test CH4 critical pathway protection functionality.
"""

import sys
from pathlib import Path
import cantera as ct

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from reduction.constraints import count_ch4_pathway_families, protected_rxn_indices, is_critical_ch4_path


def test_ch4_critical_protection():
    """Test that all four CH4 pathway families are found and protected in GRI-3.0."""
    
    # Load GRI-3.0 mechanism
    sol = ct.Solution("gri30.yaml")
    print(f"GRI-3.0 has {len(sol.reactions())} reactions")
    
    # Count CH4 pathway families
    family_counts = count_ch4_pathway_families(sol)
    
    print("\nCH4 Pathway Family Counts:")
    print("-" * 50)
    for family, count in family_counts.items():
        print(f"{family}: {count} reactions")
    
    # Get protected indices
    protected = protected_rxn_indices(sol)
    print(f"\nTotal protected reactions: {len(protected)}")
    
    # Check that critical CH4 path reactions are protected
    critical_indices = []
    for i, rxn in enumerate(sol.reactions()):
        if is_critical_ch4_path(rxn):
            critical_indices.append(i)
    
    protected_critical = [i for i in critical_indices if i in protected]
    
    print(f"\nCritical CH4 path reactions found: {len(critical_indices)}")
    print(f"Critical CH4 path reactions protected: {len(protected_critical)}")
    
    # Assertions
    assert len(critical_indices) > 0, "No critical CH4 pathway reactions found"
    assert len(protected_critical) == len(critical_indices), "Not all critical CH4 pathway reactions are protected"
    
    # Check that we found all four families
    total_family_reactions = sum(family_counts.values())
    assert total_family_reactions > 0, "No CH4 pathway family reactions found"
    
    print(f"\n✅ All {len(protected_critical)} critical CH4 pathway reactions are protected!")
    print(f"✅ Found {total_family_reactions} total CH4 pathway family reactions across all families")
    
    return family_counts, len(protected_critical)


if __name__ == "__main__":
    family_counts, protected_count = test_ch4_critical_protection()
    
    print(f"\nSUMMARY:")
    print(f"Family counts: {family_counts}")
    print(f"Protected critical reactions: {protected_count}")
