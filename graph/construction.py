import networkx as nx
from typing import Dict, Optional
from cantera import Solution


def build_species_graph(solution: Solution, weights: Dict[str, float] = None, node_scores: Optional[Dict[str, float]] = None) -> nx.DiGraph:
    """Create a species graph with chemically meaningful node features and reaction edges."""
    G = nx.DiGraph()

    for sp in solution.species_names:
        score = 1.0
        if node_scores and sp in node_scores:
            score = node_scores[sp]
        G.add_node(sp, score=score)

    for i, rxn in enumerate(solution.reactions()):
        for reactant in rxn.reactants.keys():
            for product in rxn.products.keys():
                weight = 1.0
                if node_scores:
                    w_r = G.nodes[reactant]["score"]
                    w_p = G.nodes[product]["score"]
                    weight = (w_r + w_p) / 2
                if weights and rxn.reaction_type in weights:
                    weight *= weights[rxn.reaction_type]
                G.add_edge(reactant, product, reaction=i, weight=weight)

    return G

def save_graphml(G: nx.DiGraph, path: str):
    nx.write_graphml(G, path)
