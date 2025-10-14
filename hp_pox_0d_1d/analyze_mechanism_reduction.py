#!/usr/bin/env python3
"""
Analyze the reduced GRI-3.0 mechanism to identify missing CH4→CO oxidation reactions.
"""

import cantera as ct
import yaml
from pathlib import Path
import numpy as np
from typing import List, Dict, Set, Tuple

def load_reactions_from_yaml(yaml_path: Path) -> List[Dict]:
    """Load reactions from YAML file."""
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data.get('reactions', [])

def analyze_ch4_oxidation_pathway(full_mech_path: Path, reduced_mech_path: Path) -> Dict:
    """Analyze which CH4→CO oxidation reactions were pruned."""
    
    # Load full mechanism
    full_sol = ct.Solution(str(full_mech_path))
    reduced_sol = ct.Solution(str(reduced_mech_path))
    
    # Load reaction data from YAML
    full_reactions = load_reactions_from_yaml(full_mech_path)
    reduced_reactions = load_reactions_from_yaml(reduced_mech_path)
    
    print(f"Full mechanism: {len(full_sol.reactions())} reactions")
    print(f"Reduced mechanism: {len(reduced_sol.reactions())} reactions")
    print(f"Reactions pruned: {len(full_sol.reactions()) - len(reduced_sol.reactions())}")
    
    # Key CH4 oxidation species
    ch4_pathway_species = {'CH4', 'CH3', 'CH2', 'CH2(S)', 'CH', 'C', 'CH2O', 'HCO', 'CO', 'CO2'}
    
    # Find reactions involving CH4 pathway species
    ch4_reactions_full = []
    ch4_reactions_reduced = []
    
    for i, rxn in enumerate(full_sol.reactions()):
        species_involved = set(rxn.reactants.keys()) | set(rxn.products.keys())
        if species_involved & ch4_pathway_species:
            ch4_reactions_full.append({
                'index': i,
                'equation': rxn.equation,
                'species': species_involved & ch4_pathway_species,
                'reactants': list(rxn.reactants.keys()),
                'products': list(rxn.products.keys())
            })
    
    for i, rxn in enumerate(reduced_sol.reactions()):
        species_involved = set(rxn.reactants.keys()) | set(rxn.products.keys())
        if species_involved & ch4_pathway_species:
            ch4_reactions_reduced.append({
                'index': i,
                'equation': rxn.equation,
                'species': species_involved & ch4_pathway_species,
                'reactants': list(rxn.reactants.keys()),
                'products': list(rxn.products.keys())
            })
    
    print(f"\nCH4 pathway reactions in full mechanism: {len(ch4_reactions_full)}")
    print(f"CH4 pathway reactions in reduced mechanism: {len(ch4_reactions_reduced)}")
    
    # Find pruned reactions
    reduced_equations = {r['equation'] for r in ch4_reactions_reduced}
    pruned_reactions = []
    
    for rxn in ch4_reactions_full:
        if rxn['equation'] not in reduced_equations:
            pruned_reactions.append(rxn)
    
    print(f"Pruned CH4 pathway reactions: {len(pruned_reactions)}")
    
    # Categorize pruned reactions by type
    reaction_categories = {
        'CH4_abstraction': [],
        'CH3_oxidation': [],
        'CH2_chemistry': [],
        'CH2O_chemistry': [],
        'HCO_chemistry': [],
        'CO_chemistry': [],
        'other': []
    }
    
    for rxn in pruned_reactions:
        equation = rxn['equation']
        reactants = set(rxn['reactants'])
        products = set(rxn['products'])
        
        # CH4 abstraction reactions
        if 'CH4' in reactants and 'CH3' in products:
            reaction_categories['CH4_abstraction'].append(rxn)
        # CH3 oxidation
        elif 'CH3' in reactants and any(sp in products for sp in ['CH2', 'CH2O', 'CH2(S)']):
            reaction_categories['CH3_oxidation'].append(rxn)
        # CH2 chemistry
        elif any(sp in reactants for sp in ['CH2', 'CH2(S)']):
            reaction_categories['CH2_chemistry'].append(rxn)
        # CH2O chemistry
        elif 'CH2O' in reactants or 'CH2O' in products:
            reaction_categories['CH2O_chemistry'].append(rxn)
        # HCO chemistry
        elif 'HCO' in reactants or 'HCO' in products:
            reaction_categories['HCO_chemistry'].append(rxn)
        # CO chemistry
        elif 'CO' in reactants or 'CO' in products:
            reaction_categories['CO_chemistry'].append(rxn)
        else:
            reaction_categories['other'].append(rxn)
    
    # Print analysis
    print("\n" + "="*80)
    print("PRUNED CH4→CO OXIDATION REACTIONS ANALYSIS")
    print("="*80)
    
    for category, reactions in reaction_categories.items():
        if reactions:
            print(f"\n{category.upper().replace('_', ' ')} ({len(reactions)} reactions):")
            print("-" * 50)
            for rxn in reactions:
                print(f"  {rxn['equation']}")
    
    # Analyze critical missing reactions
    print("\n" + "="*80)
    print("CRITICAL MISSING REACTIONS FOR CH4 OXIDATION")
    print("="*80)
    
    critical_missing = []
    
    # Check for CH4 + H abstraction
    ch4_h_abstraction = [r for r in pruned_reactions 
                        if 'CH4' in r['reactants'] and 'H' in r['reactants'] and 'CH3' in r['products']]
    if ch4_h_abstraction:
        critical_missing.extend(ch4_h_abstraction)
        print("\nCH4 + H abstraction reactions (CRITICAL):")
        for rxn in ch4_h_abstraction:
            print(f"  {rxn['equation']}")
    
    # Check for CH4 + OH abstraction
    ch4_oh_abstraction = [r for r in pruned_reactions 
                         if 'CH4' in r['reactants'] and 'OH' in r['reactants'] and 'CH3' in r['products']]
    if ch4_oh_abstraction:
        critical_missing.extend(ch4_oh_abstraction)
        print("\nCH4 + OH abstraction reactions (CRITICAL):")
        for rxn in ch4_oh_abstraction:
            print(f"  {rxn['equation']}")
    
    # Check for CH3 + O2 oxidation
    ch3_o2_oxidation = [r for r in pruned_reactions 
                       if 'CH3' in r['reactants'] and 'O2' in r['reactants']]
    if ch3_o2_oxidation:
        critical_missing.extend(ch3_o2_oxidation)
        print("\nCH3 + O2 oxidation reactions (CRITICAL):")
        for rxn in ch3_o2_oxidation:
            print(f"  {rxn['equation']}")
    
    # Check for CH2O decomposition
    ch2o_decomposition = [r for r in pruned_reactions 
                         if 'CH2O' in r['reactants'] and 'HCO' in r['products']]
    if ch2o_decomposition:
        critical_missing.extend(ch2o_decomposition)
        print("\nCH2O decomposition reactions (CRITICAL):")
        for rxn in ch2o_decomposition:
            print(f"  {rxn['equation']}")
    
    # Check for HCO decomposition
    hco_decomposition = [r for r in pruned_reactions 
                        if 'HCO' in r['reactants'] and 'CO' in r['products']]
    if hco_decomposition:
        critical_missing.extend(hco_decomposition)
        print("\nHCO decomposition reactions (CRITICAL):")
        for rxn in hco_decomposition:
            print(f"  {rxn['equation']}")
    
    return {
        'total_pruned': len(pruned_reactions),
        'critical_missing': critical_missing,
        'categories': reaction_categories,
        'ch4_pathway_reactions_full': len(ch4_reactions_full),
        'ch4_pathway_reactions_reduced': len(ch4_reactions_reduced)
    }

if __name__ == "__main__":
    full_mech = Path("gri30.yaml")
    reduced_mech = Path("mechanisms_reduced/gri30_reduced.yaml")
    
    if not full_mech.exists():
        print(f"Error: {full_mech} not found")
        exit(1)
    if not reduced_mech.exists():
        print(f"Error: {reduced_mech} not found")
        exit(1)
    
    results = analyze_ch4_oxidation_pathway(full_mech, reduced_mech)
    
    print(f"\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total CH4 pathway reactions pruned: {results['total_pruned']}")
    print(f"Critical missing reactions: {len(results['critical_missing'])}")
    print(f"CH4 pathway retention rate: {results['ch4_pathway_reactions_reduced']/results['ch4_pathway_reactions_full']*100:.1f}%")
