import pandas as pd

from msclinical.neda.core import Params, cohort, neda

VIS = pd.DataFrame([(p, d, 2.0) for p in "abcd" for d in ("2020-01-01", "2020-07-01", "2021-01-01")],
                   columns=["patient", "date", "edss"])
VIS.loc[(VIS.patient == "b") & (VIS.date != "2020-01-01"), "edss"] = 3.0  # b: confirmed worsening onset 2020-07
MRI = pd.DataFrame([("a", "2020-12-01", 0, 0, 0, -0.3), ("b", "2020-12-01", 0, 0, 0, -0.3),
                    ("c", "2020-12-01", 1, 0, 0, -0.3)], columns=["patient", "date", "new_t2", "enlarging_t2", "gd", "pbvc"])
MRI = pd.concat([MRI, pd.DataFrame([(p, "2020-01-01", 0, 0, 0, None) for p in "abc"], columns=MRI.columns)])
REL = pd.DataFrame([("a", "2021-06-01")], columns=["patient", "date"])


def test_components():
    t = neda(VIS, REL, MRI)
    t = t[t.window == 1].set_index("patient")
    assert t.loc["a", "neda3"] is True and t.loc["a", "neda4"] is True
    assert t.loc["b", "neda3"] is False and "CDW" in t.loc["b", "activity"]
    assert t.loc["c", "neda3"] is False and "MRI" in t.loc["c", "activity"]
    assert t.loc["d", "neda3"] is None  # no MRI in window: not assessable, never NEDA


def test_neda4_atrophy():
    m = MRI.copy()
    m.loc[(m.patient == "a") & (m.date == "2020-12-01"), "pbvc"] = -0.9
    t = neda(VIS, REL, m)
    t = t[t.window == 1].set_index("patient")
    assert t.loc["a", "neda3"] is True and t.loc["a", "neda4"] is False


def test_rebaseline_drops_early_activity():
    rel = pd.DataFrame([("a", "2020-02-01")], columns=["patient", "date"])
    assert neda(VIS, rel, MRI).set_index("patient").loc["a", "neda3"] is False
    t = neda(VIS, rel, MRI, p=Params(rebaseline_months=3))
    assert (t[t.patient == "a"]["relapses"] == 0).all()


def test_cohort_counts_unknown_separately():
    c = cohort(neda(VIS, REL, MRI))
    assert c["neda3"]["patients"] == 4 and c["neda3"]["assessable"] == 3 and c["neda3"]["neda"] == 1


def test_review_regressions():
    import numpy as np

    vis = pd.DataFrame([("a", "2019-01-01", 2.0), ("a", "2020-01-01", 2.5), ("a", "2020-06-01", 3.0), ("a", "2020-10-01", 3.0),
                        ("a", "2021-01-01", 3.0)], columns=["patient", "date", "edss"])
    mri = pd.DataFrame([("a", "2020-12-01", 0, 0, 0)], columns=["patient", "date", "new_t2", "enlarging_t2", "gd"])
    st = pd.DataFrame([("a", "2020-01-01")], columns=["patient", "date"])
    t = neda(vis, pd.DataFrame(columns=["patient", "date"]), mri, st)
    assert t.iloc[0]["cdw_events"] == 0 and t.iloc[0]["neda3"] == True
    m2 = mri.assign(new_t2=np.nan, enlarging_t2=np.nan, gd=np.nan)
    assert neda(vis, pd.DataFrame(columns=["patient", "date"]), m2, st).iloc[0]["neda3"] is None   # unread MRI
    rel = pd.DataFrame([("a", "2021-01-01")], columns=["patient", "date"])
    assert neda(vis, rel, mri, st).iloc[-1]["relapses"] == 1                  # relapse on the last date counts
    short = neda(vis[vis.date <= "2020-10-01"], pd.DataFrame(columns=["patient", "date"]), mri.assign(date="2020-09-01"), st)
    assert cohort(short, 1)["neda3"]["assessable"] == 0                       # 9 months is not one year of NEDA
