"""HP-POX Validator.

Validates HP-POX cases against Richter et al. (Fuel 2015) benchmark specifications.
Enforces oxygen-blown HP-POX feeds with proper normalization and error analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import yaml
from datetime import datetime
import cantera as ct

from reactor.plug_flow import PlugFlowReactor, ReactorGeometry, HeatTransferConfig, InletConditions


class HP_POXValidator:
    """HP-POX validator with oxygen-blown feeds per Richter benchmark."""
    
    def __init__(self, config_file: str = "config/hp_pox_benchmark.yaml"):
        """Initialize HP-POX validator."""
        self.config_file = config_file
        self.config = self._load_config()
        
        # Create output directory
        self.output_dir = Path("results_1d_hp_pox")
        self.output_dir.mkdir(exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)
        
    def _load_config(self) -> Dict:
        """Load HP-POX benchmark configuration."""
        with open(self.config_file, 'r') as f:
            return yaml.safe_load(f)
