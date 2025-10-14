## Reduction Pipeline Changes

- Added multiprocessing GA evaluation with `--workers` flag; BLAS threads capped to 1.
- Introduced constraints-only GA search (no pruning by GNN scores).
- Implemented Cantera-based YAML writer with duplicate auto-repair; preflight validator.
- Added checkpointing per generation under `ga_runs/<mech>/gen_XX.json`.
- Added CLI flags: `--workers`, `--fast-cells`, `--pop`, `--gens`, `--elite`, `--mut`, `--cx`, `--resume`.
- Kept physics identical between fast and full modes; only axial cells differ.

