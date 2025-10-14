#!/usr/bin/env python3
"""
Quick test of the CH4 pathway constraints.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import cantera as ct
from reduction.constraints import count_ch4_pathway_families, protected_rxn_indices

def test_constraints():
    """Test the CH4 oxidation pathway constraints on GRI-3.0."""
    print("Testing CH4 oxidation pathway constraints...")
    
    # Load GRI-3.0
    sol = ct.Solution("gri30.yaml")
    print(f"GRI-3.0 has {len(sol.reactions())} reactions")
    
    # Get protected indices
    protected = protected_rxn_indices(sol)
    print(f"Protected reactions: {len(protected)}")
    
    # Count CH4 pathway families
    family_counts = count_ch4_pathway_families(sol)
    
    print(f"\nCH4 Pathway Family Counts:")
    for family, count in family_counts.items():
        print(f"  {family}: {count} reactions")
    
    total_family_reactions = sum(family_counts.values())
    print(f"\nTotal CH4 pathway family reactions: {total_family_reactions}")
    
    # Check if they're all protected
    critical_count = 0
    for i, rxn in enumerate(sol.reactions()):
        from reduction.constraints import is_critical_ch4_path
        if is_critical_ch4_path(rxn):
            critical_count += 1
    
    protected_critical = 0
    for i, rxn in enumerate(sol.reactions()):
        from reduction.constraints import is_critical_ch4_path
        if is_critical_ch4_path(rxn) and i in protected:
            protected_critical += 1
    
    print(f"Critical reactions protected: {protected_critical}/{critical_count}")
    
    if protected_critical == critical_count and critical_count > 0:
        print("✅ All critical CH4 pathway reactions are protected!")
    else:
        print("❌ Some critical CH4 pathway reactions are NOT protected!")

if __name__ == "__main__":
    test_constraints()
