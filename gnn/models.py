import torch
from torch import nn
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.data import Data
import networkx as nx
import csv
from cantera import Solution
import cantera as ct

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
    mapping = {n: i for i, n in enumerate(G.nodes())}
    edges = [[mapping[u], mapping[v]] for u, v in G.edges()]
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()

    features = []
    for name in G.nodes():
        species = solution.species(name)
        Tmin = species.thermo.min_temp
        Tmax = species.thermo.max_temp
        features.append([
            G.out_degree(name),
            G.in_degree(name),
            species.composition.get("C", 0),
            species.composition.get("H", 0),
            species.composition.get("O", 0),
            Tmin,
            Tmax
        ])
    x = torch.tensor(features, dtype=torch.float)
    return Data(x=x, edge_index=edge_index)


def train_gnn(G: nx.DiGraph, solution: Solution, labels: dict[str, float] | None = None, epochs: int = 10) -> SpeciesGCN:
    data = graph_to_data(G, solution)
    model = SpeciesGCN(in_dim=data.x.size(1), hidden=16)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    if labels is None:
        y = torch.ones(len(G.nodes), 1)
    else:
        y = torch.tensor([labels.get(n, 0.0) for n in G.nodes], dtype=torch.float).unsqueeze(1)

    for _ in range(epochs):
        optimizer.zero_grad()
        out = model(data.x, data.edge_index)
        loss = nn.functional.mse_loss(out, y)
        loss.backward()
        optimizer.step()
    return model

def predict_scores(model: nn.Module, G: nx.DiGraph, solution: ct.Solution, save_path: str = "gnn_scores.csv") -> dict:
    data = graph_to_data(G, solution)
    with torch.no_grad():
        scores = model(data.x, data.edge_index).squeeze().cpu().numpy()
    score_dict = {n: float(s) for n, s in zip(G.nodes, scores)}

    # Save scores to CSV
    with open(save_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "score"])
        for species, score in score_dict.items():
            writer.writerow([species, score])

    return score_dict

