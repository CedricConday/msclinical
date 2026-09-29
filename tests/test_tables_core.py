import pandas as pd

from msclinical.tables.core import validate


def test_checks():
    t = {"patients": pd.DataFrame([("a", "2020-01-01", "2022-01-01")], columns=["patient", "start", "end"]),
         "visits": pd.DataFrame([("a", "2020-01-01", 0.5), ("a", "2020-02-10", 5.0), ("b", "2020-03-01", 2.0)], columns=["patient", "date", "edss"]),
         "relapses": pd.DataFrame([("a", "2020-02-01")], columns=["patient", "date"]),
         "mri": pd.DataFrame([("a", "2020-06-01", -1)], columns=["patient", "date", "new_t2"])}
    msgs = " | ".join(i["message"] for i in validate(t))
    assert "0.5 does not exist" in msgs and "3 points" in msgs and "after a relapse onset" in msgs
    assert "non-negative integers" in msgs and "not in patients table" in msgs


def test_trialtruth_output_is_clean():
    from msclinical.sim.core import Model, simulate

    s = simulate(Model(n_per_arm=10, seed=0))
    t = {k: s[k] for k in ("patients", "visits", "relapses", "mri")}
    assert not [i for i in validate(t) if i["level"] == "error"]



def test_review_regressions():
    v = pd.DataFrame([("a", "2020-01-01", "3,5"), ("a", "2020-06-01", "x")], columns=["patient", "date", "edss"])
    msgs = [i["message"] for i in validate({"visits": v})]
    assert any("not a number" in m for m in msgs)
    v2 = pd.DataFrame([("a", "2020-01-01", 2.0), ("a", "2020-03-01", 2.0), ("a", "2020-05-01", 2.0), ("a", "2020-07-10", 2.0)],
                      columns=["patient", "date", "edss"])
    r = pd.DataFrame([("a", "2019-01-01"), ("a", "2019-06-01"), ("a", "2020-07-01")], columns=["patient", "date"])
    rows = next(i["rows"] for i in validate({"visits": v2, "relapses": r}) if "after a relapse onset" in i["message"])
    assert rows == [3]
