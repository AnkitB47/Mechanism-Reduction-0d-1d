import pytest

from hp_pox_0d_1d.hp_pox_physics import hp_pox_model
from hp_pox_0d_1d.hp_pox_physics.config_cases import richter_case_1


def test_load_inlet_streams_uses_config_cases() -> None:
    streams = hp_pox_model.load_inlet_streams()
    assert "case1" in streams

    primary = streams["case1"]["primary_steam"]
    expected_mass = richter_case_1.primary_steam.mass_kg_per_s * 3600.0
    assert primary["mass_kg_per_h"] == pytest.approx(expected_mass)

    natural_gas = streams["case1"]["natural_gas"]
    assert "composition_volpct" in natural_gas
    assert pytest.approx(sum(natural_gas["composition_volpct"].values()), rel=1e-6, abs=1e-6) == 100.0


def test_create_inlet_state_accepts_caseconfig() -> None:
    model = hp_pox_model.HPPOXModel()
    gas_state, m_dot = model.create_inlet_state(richter_case_1, richter_case_1.case_id)
    total_mass_expected = (
        richter_case_1.primary_steam.mass_kg_per_s
        + richter_case_1.secondary_steam.mass_kg_per_s
        + richter_case_1.oxygen.mass_kg_per_s
        + richter_case_1.nitrogen.mass_kg_per_s
        + richter_case_1.natural_gas.mass_kg_per_s
    )
    assert m_dot == pytest.approx(total_mass_expected)
    assert gas_state.temperature > 0.0
