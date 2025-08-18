# === testing/run_tests.py ===

import argparse
import os
import cantera as ct
from testing.pipeline import full_pipeline


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", default="data/gri30.yaml", help="Path to mechanism YAML file")
    parser.add_argument("--out", default="results", help="Output directory")
    parser.add_argument("--steps", type=int, default=1000, help="Number of time steps")
    parser.add_argument("--tf", type=float, default=0.5, help="Final simulation time [s]")
    parser.add_argument("--phi", type=float, default=1.0, help="Equivalence ratio")
    parser.add_argument(
        "--preset",
        default="methane_air",
        choices=["methane_air", "HR1", "HR2", "HR3"],
        help="Mixture preset",
    )
    parser.add_argument("--T0", type=float, default=1500.0, help="Initial temperature [K]")
    parser.add_argument("--p0", type=float, default=float(ct.one_atm), help="Initial pressure [Pa]")
    parser.add_argument("--log-times", action="store_true", help="Use log-spaced time sampling")
    parser.add_argument("--isothermal", action="store_true", help="Use isothermal reactor for HR presets")

    # === NEW ARGUMENTS for label building and aggressiveness ===
    parser.add_argument("--alpha", type=float, default=0.8, help="Weight between LOO and flux importance")
    parser.add_argument("--tf-short", type=float, default=None, help="Shortened tf for LOO runs")
    parser.add_argument("--steps-short", type=int, default=None, help="Steps for LOO runs")
    parser.add_argument("--no-cache-labels", action="store_true", help="Recompute labels instead of loading cache")
    parser.add_argument("--min-species", type=int, default=None, help="Minimum species allowed in GA")
    parser.add_argument("--max-species", type=int, default=None, help="Maximum species allowed in GA")
    parser.add_argument("--target-species", type=int, default=None, help="Target species count")
    parser.add_argument("--size-weight", type=float, default=3.0, help="Penalty weight for size term")
    parser.add_argument("--species-weight", type=float, default=2.0, help="Penalty weight for species residual")
    parser.add_argument("--zoom-around-tau", type=float, default=None, help="Decades half-span around τ")
    parser.add_argument("--ymax", type=float, default=None, help="Optional upper y-limit for species plot")

    args = parser.parse_args()

    full_pipeline(
        args.mechanism,
        args.out,
        steps=args.steps,
        tf=args.tf,
        phi=args.phi,
        preset=args.preset,
        T0=args.T0,
        p0=args.p0,
        log_times=args.log_times,
        alpha=args.alpha,
        tf_short=args.tf_short,
        steps_short=args.steps_short,
        cache_labels=not args.no_cache_labels,
        isothermal=args.isothermal,
        min_species=args.min_species,
        max_species=args.max_species,
        target_species=args.target_species,
        size_weight=args.size_weight,
        species_weight=args.species_weight,
        zoom_around_tau=args.zoom_around_tau,
        ymax=args.ymax,
    )
    print(f"\n✅ Pipeline complete. Results written to '{args.out}'\n")


if __name__ == "__main__":
    main()
