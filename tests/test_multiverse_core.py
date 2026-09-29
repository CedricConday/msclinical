import pandas as pd

from msclinical.multiverse.core import influence, run


def test_grid_and_influence():
    vis = pd.DataFrame([(p, d, e) for p, seq in (("a", (2, 3, 3, 3)), ("b", (2, 2, 3, 3)))
                        for d, e in zip(("2020-01-01", "2020-04-01", "2020-08-01", "2021-01-01"), seq)], columns=["patient", "date", "edss"])
    rel = pd.DataFrame([("a", "2020-02-10")], columns=["patient", "date"])
    fs = pd.DataFrame([("a", "2019-01-01"), ("b", "2019-01-01")], columns=["patient", "date"])
    t = run(vis, rel, None, fs)
    assert len(t) == 32
    assert t[(t.raw_days == 90) & (t.confirm_days == 84) & (t.reference == "roving")]["raw"].max() == 1
    assert t[(t.raw_days == 30)]["raw"].max() == 0  # relapse 51 days before onset: RAW under 90 d, PIRA under 30 d
    assert set(influence(t)["choice"]) >= {"confirm_days", "raw_days"}
