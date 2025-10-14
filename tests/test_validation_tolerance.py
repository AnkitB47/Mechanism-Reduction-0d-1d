from pathlib import Path

import yaml

from hp_pox_0d_1d.tools.config_templates import write_full_config


def test_generated_config_has_relaxed_h2co_tol(tmp_path: Path) -> None:
    out_cfg = tmp_path / "full.yaml"
    dummy_mech = tmp_path / "dummy_mech.yaml"

    write_full_config(
        template_src=None,
        out_path=out_cfg,
        mechanism_path=dummy_mech,
        cases_root=tmp_path / "cases",
        output_root=tmp_path / "out",
    )

    data = yaml.safe_load(out_cfg.read_text())
    assert data["validation"]["h2co_target"] == 2.0
    assert data["validation"]["h2co_rel_tol"] >= 0.05
