"""CLI entry point for running the reduction testing pipeline."""

import argparse

from testing.pipeline import full_pipeline


def main():
    """Parse command line arguments and run the pipeline."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", required=True, help="Path to mechanism YAML file")
    parser.add_argument("--out", default="results", help="Output directory")
    args = parser.parse_args()
    full_pipeline(args.mechanism, args.out, steps=50)

if __name__ == '__main__':
    main()
