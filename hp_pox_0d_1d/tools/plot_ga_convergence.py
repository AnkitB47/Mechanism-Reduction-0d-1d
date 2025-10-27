#!/usr/bin/env python3
"""
Plot GA convergence from ga_trace.csv.

Usage example:
  python hp_pox_0d_1d/tools/plot_ga_convergence.py \
      --trace runs/20251013_234429/ga_trace.csv \
      --out runs/20251013_234429/ga_convergence.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd


def load_trace(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Trace file not found: {path}")
    df = pd.read_csv(path)
    # older runs may have unexpected whitespace in headers
    df.columns = [c.strip() for c in df.columns]
    if "generation" not in df.columns or "penalty" not in df.columns:
        raise KeyError(f"Trace file {path} must contain 'generation' and 'penalty' columns.")
    df = df.sort_values("generation")
    df["penalty"] = df["penalty"].astype(float)
    return df


def plot_convergence(df: pd.DataFrame, title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(df["generation"], df["penalty"], marker="o", lw=2)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness penalty")
    ax.set_title(title)
    ax.grid(True, which="both", ls="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot GA convergence from ga_trace.csv.")
    parser.add_argument(
        "--trace",
        nargs="+",
        required=True,
        help="One or more ga_trace.csv paths. If multiple are provided, curves are overlaid.",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output plot path (PNG). Defaults to <trace_dir>/ga_convergence.png.",
    )
    parser.add_argument(
        "--y-scale",
        choices=["linear", "log", "symlog"],
        default="log",
        help="Y-axis scaling (default: log).",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom plot title. Default uses name of first trace directory.",
    )
    args = parser.parse_args()

    trace_paths = [Path(p) for p in args.trace]
    data_frames = [load_trace(p) for p in trace_paths]

    if args.out:
        out_path = Path(args.out)
    else:
        out_path = trace_paths[0].parent / "ga_convergence.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    title = args.title or f"GA convergence ({trace_paths[0].parent.name})"

    fig, ax = plt.subplots(figsize=(9, 4.8))
    for path, df in zip(trace_paths, data_frames):
        label = path.parent.name
        ax.plot(df["generation"], df["penalty"], marker="o", lw=2, label=label)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness penalty")
    ax.set_title(title)
    ax.set_yscale(args.y_scale)
    ax.grid(True, which="both", ls="--", alpha=0.3)
    ax.legend(title="Run", fontsize="small")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    print(f"[plot_ga_convergence] wrote {out_path.as_posix()}")


if __name__ == "__main__":
    main()
