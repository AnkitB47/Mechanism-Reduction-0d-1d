#!/usr/bin/env python3
"""
Smoke tests for the CLI orchestrator (dry-run mode for speed).
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CLI_PATH = REPO_ROOT / "hp_pox_0d_1d" / "tools" / "run_full_vs_reduced.py"


def test_cli_dry_run(tmp_path):
    """Ensure the CLI wiring works in dry-run mode."""
    comparison_dir = tmp_path / "comparisons"
    full_results = tmp_path / "full_results"
    reduced_results = tmp_path / "reduced_results"
    full_config = tmp_path / "cfg" / "gri_full.yaml"

    cmd = [
        sys.executable,
        str(CLI_PATH),
        "--mech",
        "gri30",
        "--cases",
        "A_case1_richter",
        "--out",
        str(comparison_dir),
        "--full-config",
        str(full_config),
        "--full-results-dir",
        str(full_results),
        "--reduced-results-dir",
        str(reduced_results),
        "--skip-full",
        "--resume",
        "--dry-run",
    ]

    result = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    reduced_config = full_config.with_name("gri30_reduced.yaml")
    assert full_config.exists(), "Full configuration should be generated in dry-run"
    assert reduced_config.exists(), "Reduced configuration should be generated in dry-run"

    dry_run_dirs = sorted((REPO_ROOT / "runs").glob("dry_run_*"))
    assert dry_run_dirs, "Dry-run should create a stub runs/ directory"
    summary_json = dry_run_dirs[-1] / "summary.json"
    assert summary_json.exists(), "Dry-run summary.json missing"

    summary_csv = comparison_dir / "summary_metrics.csv"
    assert summary_csv.exists(), "Dry-run comparison summary missing"

    # Clean up stub artefacts to keep the repository tidy
    try:
        full_config.unlink(missing_ok=True)
        reduced_config.unlink(missing_ok=True)
        shutil.rmtree(dry_run_dirs[-1], ignore_errors=True)
    finally:
        pass
