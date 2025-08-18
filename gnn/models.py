"""GNN models for ranking chemical species.

This module now includes richer logging and safeguards to avoid degenerate
training targets, alongside slightly more informative node features.
"""

import csv
import logging

import cantera as ct
import networkx as nx
import numpy as np
import torch
from torch import nn
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv
from cantera import Solution

logger = logging.getLogger(__name__)

class SpeciesGCN(nn.Module):
    def __init__(self, in_dim: int, hidden: int):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hidden)
        self.conv2 = GCNConv(hidden, hidden)
        self.conv3 = GCNConv(hidden, 1)

    def forward(self, x, edge_index):
        h = torch.relu(self.conv1(x, edge_index))
        h = torch.relu(self.conv2(h, edge_index))
        return torch.sigmoid(self.conv3(h, edge_index))

def graph_to_data(G: nx.DiGraph, solution: Solution) -> Data:
    """Convert a NetworkX graph to a PyG ``Data`` instance.

    The graph nodes **must** correspond to species names in ``solution``.
    Node features include normalised in/out degrees and a few bounded chemical
    proxies such as molecular weight. Features are kept roughly in ``[0, 1]`` to
    aid optimisation stability.
    """

    mapping = {n: i for i, n in enumerate(G.nodes())}
    edges = [[mapping[u], mapping[v]] for u, v in G.edges()]
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()

    # degree normalisation safeguards against constant features
    out_deg = dict(G.out_degree())
    in_deg = dict(G.in_degree())
    max_out = max(out_deg.values()) or 1
    max_in = max(in_deg.values()) or 1

    features = []
    for name in G.nodes():
        if name not in solution.species_names:
            raise ValueError(f"Graph node '{name}' not in mechanism species list")
        species = solution.species(name)
        features.append([
            out_deg[name] / max_out,
            in_deg[name] / max_in,
            species.molecular_weight / 100.0,  # bounded proxy
            species.composition.get("C", 0) / 10.0,
            species.composition.get("H", 0) / 20.0,
            species.composition.get("O", 0) / 10.0,
        ])

    x = torch.tensor(features, dtype=torch.float)
    return Data(x=x, edge_index=edge_index)


def train_gnn(
    G: nx.DiGraph,
    solution: Solution,
    labels: dict[str, float] | None = None,
    epochs: int = 200,
) -> SpeciesGCN:
    """Train a small GCN to predict species importance.

    The training targets are rescaled to ``[0, 1]`` and a warning is raised if
    too many zeros are supplied which would lead to a degenerate model.
    """

    data = graph_to_data(G, solution)
    model = SpeciesGCN(in_dim=data.x.size(1), hidden=32)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    if labels is None:
        y = torch.ones(len(G.nodes), 1)
    else:
        y = torch.tensor([labels.get(n, 0.0) for n in G.nodes], dtype=torch.float).unsqueeze(1)

    pairs = list(zip(G.nodes, y.squeeze().tolist()))
    for sp, target in pairs[:10]:
        logger.info("train pair %s -> %.3f", sp, target)

    zero_frac = float((y == 0.0).sum()) / y.numel()
    if zero_frac > 0.05:
        raise AssertionError(f"Too many zero labels ({zero_frac*100:.1f}% > 5%)")

    y_min, y_max = y.min(), y.max()
    if torch.isclose(y_max, y_min):
        y = torch.full_like(y, 0.5)
    else:
        y = (y - y_min) / (y_max - y_min + 1e-8)

    for epoch in range(epochs):
        optimizer.zero_grad()
        out = model(data.x, data.edge_index)
        loss = nn.functional.mse_loss(out, y)
        loss.backward()
        optimizer.step()
        if epoch % 10 == 0:
            logger.info("epoch %03d loss %.4f", epoch, loss.item())

    return model

def predict_scores(
    model: nn.Module,
    G: nx.DiGraph,
    solution: ct.Solution,
    save_path: str = "gnn_scores.csv",
) -> dict:
    """Predict species scores and persist them to ``save_path``.

    Logged statistics include min/mean/max as well as a small histogram to help
    diagnose training collapse. An assertion guards against an all-zero output.
    """

    data = graph_to_data(G, solution)
    with torch.no_grad():
        scores = model(data.x, data.edge_index).squeeze().cpu().numpy()

    min_s, mean_s, max_s = float(scores.min()), float(scores.mean()), float(scores.max())
    logger.info("score stats min=%.3f mean=%.3f max=%.3f", min_s, mean_s, max_s)
    hist, _ = np.histogram(scores, bins=5, range=(0, 1))
    logger.info("score histogram %s", hist.tolist())

    if float(np.abs(scores).sum()) == 0.0:
        raise AssertionError("All predicted scores are zero")

    score_dict = {n: float(s) for n, s in zip(G.nodes, scores)}

    with open(save_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "score"])
        for species, score in score_dict.items():
            writer.writerow([species, score])

    return score_dict

