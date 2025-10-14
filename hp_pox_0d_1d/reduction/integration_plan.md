# Integration Plan: GA+GNN → hp_pox_0d_1d/reduction

This plan maps the legacy modules discovered under `D:/Thesis/Master_Thesis` to thin adapters in `hp_pox_0d_1d/reduction`.

## Legacy entry points

- GA (binary genomes): `metaheuristics/ga.py`
  - `run_ga(genome_length, eval_fn, options, fixed_indices=None, mask=None)`
- Mechanism wrapper: `mechanism/loader.py`
  - `Mechanism(file).species_names()`, `.reactions()`, `.remove_species(names)`, `.save(out_path)`
- Species graph + GNN: `graph/construction.py`, `gnn/models.py`
  - `build_species_graph(solution) -> nx.DiGraph`
  - `train_gnn(G, solution, labels, epochs) -> model`
  - `predict_scores(model, G, solution, save_path) -> dict species->score`

## Adapters to implement/use

- `reduction/dataset.py`
  - Build manifest from axial CSVs (z_m, T_C, X_*) + KPIs.
  - Record mechanism path and ordered reaction indices for alignment.
  - For GNN features, construct species graph via Cantera `Solution` loaded from mechanism.

- `reduction/gnn_bridge.py`
  - Import: `from gnn.models import train_gnn, predict_scores`
  - Import: `from graph.construction import build_species_graph`
  - Steps:
    1) Load mechanism with `cantera.Solution(mech_path)`
    2) Build `G = build_species_graph(solution)`
    3) Train: `model = train_gnn(G, solution, labels=None, epochs=config.train.epochs)`
    4) Score: `scores = predict_scores(model, G, solution, save_path)`
    5) Convert species scores → reaction scores by aggregating edge endpoints per-reaction id stored in Graph edges.
    6) Save `rxn_scores.npy` aligned to `solution.reactions()` ordering.

- `reduction/ga_bridge.py`
  - Import: `from metaheuristics.ga import run_ga, GAOptions`
  - Define eval_fn(mask):
    - Construct reaction keep-mask (respect constraints), write pruned YAML via `mechanism_io.apply_reaction_mask`
    - Call `evaluator.evaluate(mech_tmp, cases_root, fast_cases, outdir_tmp)` to obtain fitness (negative blended error + size bonus)
  - Run GA with parsed budget; on `eval_every`, persist topK candidates and run full evaluator.

- `reduction/mechanism_io.py`
  - Load YAML, prune reactions by indices to KEEP, preserve species blocks and metadata.
  - Ensure no orphan species; if created, re-add species or reject candidate.

- `reduction/constraints.py`
  - Protect reaction types: `falloff`, `three-body`, `duplicate`.
  - Protect mandatory species: `H2, CO, CO2, CH4, H2O, O2, N2` and any produced-only species.
  - Provide `filtered_mask(mask)->mask` and `fixed_indices()` for GA.

- `reduction/evaluator.py`
  - Fast mode: fewer PFR cells; full mode: original resolution.
  - KPIs: H2_out, CO_out, CO2_out, CH4_out, H2/CO, T_out_K, CH4_slip_pct, deltas vs targets.

## Data mapping

- From `hp_pox_results_{mech}/{case}/axial_profiles.csv`:
  - Use: `z_m`, `T_C`, `X_H2`, `X_CO`, `X_CO2`, `X_CH4` at taps {5%, 20%, 50%, 100%} for optional labels.
- Mechanism reaction ordering: `solution.reactions()` index defines `rxn_scores.npy` alignment.

## CLI mapping

- `prepare-data`: builds dataset directory with `manifest.json`, feature parquet(s), and mechanism path.
- `train-gnn`: trains/loads GNN; writes `metrics.json`, `gnn.ckpt`.
- `score`: predicts reaction scores; writes `rxn_scores.npy`.
- `reduce`: runs GA using `rxn_scores.npy` and writes `reduced_mech.yaml` + `scorecard.json`.
- `validate`: runs full PSR→PFR for 4 cases; writes `metrics.json` and overlay plots.

## Open questions

1) Confirm where the reaction index is exposed in `graph/construction.py` (edge attribute `reaction`).
2) Provide expected label schema for `train_gnn` (currently optional).
3) Confirm any additional constraints (e.g., pressure-dependent reactions) to hard-protect.
4) Confirm preferred fitness blending weights (surrogate vs hard evaluation).
