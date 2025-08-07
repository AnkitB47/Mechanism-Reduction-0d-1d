# === testing/run_tests.py ===

import argparse
import os
from testing.pipeline import full_pipeline


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", default="data/gri30.yaml", help="Path to mechanism YAML file")
    parser.add_argument("--out", default="results", help="Output directory")
    parser.add_argument("--steps", type=int, default=200, help="Number of time steps")
    parser.add_argument("--tf", type=float, default=1.0, help="Final simulation time [s]")
    args = parser.parse_args()

    full_pipeline(args.mechanism, args.out, steps=args.steps, tf=args.tf)
    print(f"\n✅ Pipeline complete. Results written to '{args.out}'\n")


if __name__ == "__main__":
    main()
