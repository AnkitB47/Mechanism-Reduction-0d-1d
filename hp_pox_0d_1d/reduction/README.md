# Reduction (GA+GNN) for HP-POX Mechanisms

This package wires a GA+GNN reducer (hosted in ../master_thesis) into the HP-POX 0D-1D pipeline. For now, the GA/GNN hooks are stubbed; the CLI and evaluators run end-to-end so you can integrate legacy components incrementally.

## Quickstart

1) Prepare dataset
```
python -m hp_pox_0d_1d.reduction.cli prepare-data \
  --cases hp_pox_results_sandiego A_case1_richter A_case4_richter B_case1_134L B_case4_134L \
  --out hp_pox_0d_1d/reduction_runs/exp1/dataset
```

2) Train/load GNN (stub)
```
python -m hp_pox_0d_1d.reduction.cli train-gnn \
  --data hp_pox_0d_1d/reduction_runs/exp1/dataset \
  --config hp_pox_0d_1d/reduction/configs/gnn.yaml \
  --out hp_pox_0d_1d/reduction_runs/exp1/gnn
```

3) Reduce (demo wiring; copies mechanism as reduced)
```
python -m hp_pox_0d_1d.reduction.cli reduce \
  --mech SAN_DIEGO_MECH/c0_mech.yaml \
  --cases hp_pox_results_sandiego A_case1_richter A_case4_richter \
  --budget "generations=3,pop=8,topK=1,eval_every=2" \
  --constraints hp_pox_0d_1d/reduction/configs/constraints.yaml \
  --out hp_pox_0d_1d/reduction_runs/exp1/runs/sandiego
```

4) Validate reduced mechanism
```
python -m hp_pox_0d_1d.reduction.cli validate \
  --mech hp_pox_0d_1d/reduction_runs/exp1/runs/sandiego/reduced_mech.yaml \
  --cases hp_pox_results_sandiego A_case1_richter A_case4_richter \
  --out hp_pox_0d_1d/reduction_runs/exp1/validate_sandiego
```

Artifacts are written to hp_pox_0d_1d/reduction_runs/<exp>/.

## Notes
- If ../master_thesis is not present, GNN and GA steps run in stub mode. Integrate your legacy modules by replacing the stubs in gnn_bridge.py and wiring GA sampling in runner.py.
- This does not modify the existing PSR/PFR APIs and plotting.


