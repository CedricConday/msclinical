import pandas as pd

from msclinical.epochs.core import Rules, arr_by_drug, epochs


def test_epochs_and_arr():
    tx = pd.DataFrame([("a", "X", "2018-01-01", "2019-06-01"), ("a", "X", "2019-06-20", "2020-01-01"),
                       ("a", "Y", "2020-03-01", None)], columns=["patient", "drug", "start", "stop"])
    fu = pd.DataFrame([("a", "2017-01-01", "2022-03-01")], columns=["patient", "start", "end"])
    ep = epochs(tx, fu, Rules())
    assert ep["drug"].tolist() == ["none", "X", "none", "Y"]               # gap of 19 days bridged; untreated gap before Y
    x = ep[ep.drug == "X"].iloc[0]
    assert str(x["window_start"]) == "2018-07-01" and str(x["window_end"]) == "2020-03-01" and x["cut"] == "stop + carry-over"
    rel = pd.DataFrame([("a", "2018-03-01"), ("a", "2019-01-01"), ("a", "2021-01-01")], columns=["patient", "date"])
    arr = arr_by_drug(rel, ep)
    assert arr.loc["X", "relapses"] == 1 and arr.loc["Y", "relapses"] == 1   # March 2018 falls in X's re-baselining



def test_review_regressions():
    fu = pd.DataFrame([("a", "2020-01-01", "2024-12-31")], columns=["patient", "start", "end"])
    tx = pd.DataFrame([("a", "X", "2020-01-01", "2020-12-31")], columns=["patient", "drug", "start", "stop"])
    ep = epochs(tx, fu)
    assert ep["drug"].tolist() == ["X", "none"] and str(ep.iloc[-1]["window_end"]) == "2024-12-31"   # trailing time kept
    tx2 = pd.DataFrame([("a", "X", "2020-01-01", "2022-01-01"), ("a", "Y", "2020-06-01", "2020-09-01")],
                       columns=["patient", "drug", "start", "stop"])
    ep2 = epochs(tx2, fu)
    assert ep2["drug"].tolist() == ["X", "Y", "X", "none"] and ep2.iloc[1]["cut"] == "X resumes"       # nested course
