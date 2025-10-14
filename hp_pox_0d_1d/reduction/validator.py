import cantera as ct


def assert_loadable(mech_path: str, transport_model: str = 'mixture-averaged') -> None:
    try:
        _ = ct.Solution(mech_path, transport_model=transport_model)
    except Exception as e:
        raise RuntimeError(f"Mechanism not loadable: {mech_path}\n{e}") from e


