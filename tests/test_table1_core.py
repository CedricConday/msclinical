from msclinical.sim.core import Model, simulate
from msclinical.table1.core import table1


def test_table_from_simulated_trial():
    s = simulate(Model(n_per_arm=50, seed=0))
    p = s["patients"].assign(birth_year=1985, sex=["F", "M"] * 50, course="RRMS")
    t = table1({"patients": p, "visits": s["visits"], "relapses": s["relapses"], "mri": s["mri"]})
    assert list(t.columns) == ["All", "active", "placebo"]
    assert t.loc["Patients", "All"] == "100" and t.loc["Female", "All"].startswith("50 (50 %)")
    assert "patient-years" in t.loc["Relapses on study, ARR", "placebo"]



def test_review_regressions():
    import pandas as pd

    p = pd.DataFrame([("a", "arm", "2020-01-01", "2021-01-01")], columns=["patient", "arm", "start", "end"])
    v = pd.DataFrame([("a", "2020-01-01", 2.0), ("a", "2020-02-01", 3.5), ("a", "2020-03-15", 4.0)], columns=["patient", "date", "edss"])
    m = pd.DataFrame([("a", "2019-06-01", 3)], columns=["patient", "date", "new_t2"])
    t = table1({"patients": p, "visits": v, "mri": m})
    assert t.loc["Baseline EDSS", "All"].startswith("2.0") and "[n=0 with MRI]" in t.loc["Any new T2 lesion on study", "All"]
