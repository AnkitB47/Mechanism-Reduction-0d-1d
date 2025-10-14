from pathlib import Path

from hp_pox_0d_1d.hp_pox_physics import hp_pox_model


def test_load_inlet_streams_accepts_str_and_path(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("inlet: {}\n")

    # Should not raise for str or Path
    hp_pox_model.load_inlet_streams(str(cfg))
    hp_pox_model.load_inlet_streams(cfg)
