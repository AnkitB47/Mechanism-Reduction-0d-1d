from pathlib import Path
import json
import numpy as np

import cantera as ct

try:
    from graph.construction import build_species_graph
    from gnn.models import train_gnn, predict_scores
except Exception:
    build_species_graph = None
    train_gnn = None
    predict_scores = None


def train_or_load_gnn(data_dir: Path, config_path: Path, outdir: Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)

    # Load GNN hyperparameters from config
    import yaml
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    gnn_params = config.get('gnn_hyperparams', {})
    epochs = gnn_params.get('epochs', 120)
    patience = gnn_params.get('patience', 20)
    lr = gnn_params.get('lr', 3e-4)

    # In this integration, we train once per-mechanism on demand during scoring; here we just persist config
    with open(outdir / 'metrics.json', 'w') as f:
        json.dump({
            'status': 'ready',
            'config': str(config_path.resolve()),
            'gnn_params': gnn_params
        }, f, indent=2)
    return 0


def _reaction_scores_from_species(solution: ct.Solution, species_scores: dict[str, float]) -> np.ndarray:
    rxns = solution.reactions()
    scores = np.zeros(len(rxns), dtype=float)
    for i, r in enumerate(rxns):
        involved = set(r.reactants.keys()) | set(r.products.keys())
        if not involved:
            scores[i] = 0.0
            continue
        vals = [species_scores.get(s, 0.0) for s in involved]
        scores[i] = float(np.max(vals)) if len(vals) else 0.0
    return scores


def score_reactions(data_dir: Path, workdir: Path, mech_path: Path) -> int:
    workdir.mkdir(parents=True, exist_ok=True)
    # Load mechanism and build species graph
    sol = ct.Solution(str(mech_path))
    if build_species_graph is None or train_gnn is None or predict_scores is None:
        # Fallback: simple degree-centrality-like heuristic using counts across reactions
        involvement = {sp: 0 for sp in sol.species_names}
        for r in sol.reactions():
            for s in set(r.reactants.keys()) | set(r.products.keys()):
                if s in involvement:
                    involvement[s] += 1
        # Normalize to [0,1]
        max_cnt = max(involvement.values()) or 1
        species_scores = {s: involvement[s] / max_cnt for s in involvement}
    else:
        G = build_species_graph(sol)
        # Train without explicit labels as default; labels can be added by extending dataset parsing
        model = train_gnn(G, sol, labels=None, epochs=200)
        species_scores = predict_scores(model, G, sol, save_path=str((workdir / 'gnn_species_scores.csv').resolve()))

    rxn_scores = _reaction_scores_from_species(sol, species_scores)
    out_npy = workdir / 'rxn_scores.npy'
    np.save(out_npy, rxn_scores)
    print(f"Reaction scores saved: {out_npy.resolve().as_posix()} (n={len(rxn_scores)})")
    return 0


