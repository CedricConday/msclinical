from msclinical.edssnoise.core import estimate, reversal_rate
from msclinical.sim.core import Model, simulate


def test_reversals_rise_with_noise():
    a = reversal_rate(simulate(Model(n_per_arm=100, rating_noise_p=0.0, seed=5))["visits"])[0]
    b = reversal_rate(simulate(Model(n_per_arm=100, rating_noise_p=0.3, seed=5))["visits"])[0]
    assert a < 0.02 < b


def test_recovers_known_noise():
    v = simulate(Model(n_per_arm=150, rating_noise_p=0.2, seed=11))["visits"]
    r = estimate(v, boot=40)
    assert abs(r["misrating_p"] - 0.2) < 0.06 and r["ci95"][0] <= r["misrating_p"] <= r["ci95"][1]


def test_sparse_data_raises():
    import pandas as pd
    import pytest

    v = pd.DataFrame([("a", "2020-01-01", 2.0), ("a", "2021-01-01", 2.0)], columns=["patient", "date", "edss"])
    with pytest.raises(ValueError):
        estimate(v, boot=5)
