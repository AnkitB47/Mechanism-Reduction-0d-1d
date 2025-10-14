#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Thin wrapper to run the HP-POX model with package-relative imports")
    parser.add_argument("--config", required=True, help="Path to configuration YAML")
    return parser.parse_known_args()


def main() -> None:
    args, rest = parse_args()
    module = "hp_pox_0d_1d.hp_pox_physics.hp_pox_model"
    cmd = [sys.executable, "-m", module, "--config", args.config, *rest]
    print("[wrapper]", " ".join(cmd))
    sys.exit(subprocess.call(cmd, cwd=REPO_ROOT))


if __name__ == "__main__":
    main()
