from testing.pipeline import full_pipeline


def test_pipeline_runs(tmp_path):
    out = tmp_path / "out"
    full_pipeline("data/gri30.yaml", str(out), steps=10)
    assert (out / "metrics.csv").exists()
