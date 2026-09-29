import pandas as pd

from msclinical.relapses.core import check


def test_rules():
    r = pd.DataFrame([("a", "2020-01-01", False, 72), ("a", "2020-01-20", False, 72), ("a", "2020-03-01", True, 72),
                      ("a", "2020-05-01", False, 12), ("a", "2020-07-01", False, 48), ("a", "2020-07-01", False, 48),
                      ("b", "2019-01-01", False, 48)], columns=["patient", "date", "fever", "duration_hours"])
    f = pd.DataFrame([("a", "2019-06-01", "2021-01-01"), ("b", "2019-06-01", "2021-01-01")], columns=["patient", "start", "end"])
    clean, removed, s = check(r, f)
    assert [str(d) for d in clean["date"]] == ["2020-01-01", "2020-07-01"]
    assert s["removed"] == 5 and set(removed["reason"].str.split(" ").str[0]) == {"within", "pseudo", "shorter", "duplicate", "outside"}


def test_review_regressions():
    import pytest

    r = pd.DataFrame([("p", "2020-01-01", True), ("p", "2020-01-01", False)], columns=["patient", "date", "fever"])
    assert len(check(r)[0]) == 1 and len(check(r.iloc[::-1])[0]) == 1       # order no longer decides
    r2 = pd.DataFrame([("p", "2020-01-01", "no"), ("p", "2020-06-01", "no")], columns=["patient", "date", "fever"])
    assert len(check(r2)[0]) == 2                                             # "no" is not truthy
    with pytest.raises(ValueError):
        check(r2.assign(fever="maybe"))
    r3 = pd.DataFrame([("p", "05.01.2020"), ("p", "01.02.2020")], columns=["patient", "date"])
    assert len(check(r3, dayfirst=True)[0]) == 1 and len(check(r3)[0]) == 2   # 5 Jan and 1 Feb are 27 days apart: one attack
    f = pd.DataFrame([("q", "2019-01-01", "2021-01-01")], columns=["patient", "start", "end"])
    assert check(r2, f)[1]["reason"].str.contains("not in follow-up").all()
