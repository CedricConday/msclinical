import pandas as pd

from msclinical.edssplus.core import events
from msclinical.pirma.core import attribute, summary

VIS = pd.DataFrame([("a", "2020-01-01", 2.0), ("a", "2020-06-01", 3.0), ("a", "2020-10-01", 3.0),
                    ("b", "2020-01-01", 2.0), ("b", "2020-06-01", 3.0), ("b", "2020-10-01", 3.0),
                    ("c", "2020-01-01", 2.0), ("c", "2020-06-01", 3.0), ("c", "2020-10-01", 3.0)],
                   columns=["patient", "date", "edss"])
EV = events(VIS)


def test_raw_pira_pirma():
    rel = pd.DataFrame([("a", "2020-04-15")], columns=["patient", "date"])
    mri = pd.DataFrame([("b", "2020-05-01", 0, 0, 0), ("c", "2020-05-01", 1, 0, 0)],
                       columns=["patient", "date", "new_t2", "enlarging_t2", "gd"])
    fs = pd.DataFrame([(p, "2019-01-01") for p in "abc"], columns=["patient", "date"])
    att = attribute(EV, rel, mri, fs).set_index("patient")
    assert att.loc["a", "class"] == "RAW"
    assert att.loc["b", "class"] == "PIRMA"
    assert att.loc["c", "class"] == "PIRA_MRI_active"
    s = summary(att.reset_index())
    assert s["pira_share"] == 2 / 3 and s["pirma_share_of_pira_with_mri"] == 0.5


def test_no_scan_is_unknown_not_pirma():
    fs = pd.DataFrame([(p, "2019-01-01") for p in "abc"], columns=["patient", "date"])
    att = attribute(EV, pd.DataFrame(columns=["patient", "date"]), pd.DataFrame(columns=["patient", "date", "new_t2"]), fs)
    assert set(att["class"]) == {"PIRA_MRI_unknown"}


def test_relapse_coverage_required():
    late = pd.DataFrame([(p, "2020-05-01") for p in "abc"], columns=["patient", "date"])
    att = attribute(EV, pd.DataFrame(columns=["patient", "date"]), None, late)
    assert set(att["class"]) == {"undetermined"}  # ascertainment starts after onset - 90 d


def test_relapse_during_confirmation_is_raw():
    rel = pd.DataFrame([("b", "2020-08-01")], columns=["patient", "date"])
    fs = pd.DataFrame([(p, "2019-01-01") for p in "abc"], columns=["patient", "date"])
    att = attribute(EV, rel, None, fs).set_index("patient")
    assert att.loc["b", "class"] == "RAW" and att.loc["a", "class"] == "PIRA_MRI_unknown"


def test_unread_mri_is_unknown_not_pirma():
    fs = pd.DataFrame([(p, "2019-01-01") for p in "abc"], columns=["patient", "date"])
    mri = pd.DataFrame([("b", "2020-05-01", None, None, None)], columns=["patient", "date", "new_t2", "enlarging_t2", "gd"])
    att = attribute(EV, pd.DataFrame(columns=["patient", "date"]), mri, fs).set_index("patient")
    assert att.loc["b", "class"] == "PIRA_MRI_unknown"


def test_misnamed_mri_columns_raise():
    import pytest

    with pytest.raises(ValueError):
        attribute(EV, pd.DataFrame(columns=["patient", "date"]), pd.DataFrame([("b", "2020-05-01", 3)], columns=["patient", "date", "new_lesions"]))


def test_default_coverage_is_first_visit_and_duplicates_ok():
    vis = pd.DataFrame([("a", d, e) for d, e in (("2018-01-01", 2), ("2018-07-01", 3), ("2018-10-01", 3), ("2018-12-01", 4), ("2019-03-15", 4))],
                       columns=["patient", "date", "edss"])
    ev = events(vis)
    first = pd.to_datetime(vis["date"]).groupby(vis["patient"]).min()
    att = attribute(ev, pd.DataFrame(columns=["patient", "date"]), None, None, visits_first=first)
    assert set(att["class"]) == {"PIRA_MRI_unknown"}
    fs = pd.DataFrame([("a", "2017-01-01"), ("a", "2016-01-01")], columns=["patient", "date"])
    assert len(attribute(ev, pd.DataFrame(columns=["patient", "date"]), None, fs)) == len(att)
