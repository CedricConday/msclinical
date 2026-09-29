import pandas as pd

from msclinical.cogworse.core import Rules, events


def V(seq, start="2020-01-01", months=6):
    d = pd.date_range(start, periods=len(seq), freq=f"{months}MS")
    return pd.DataFrame({"patient": "a", "date": d.date, "sdmt": seq})


def test_rules():
    assert events(V([50, 54, 49, 49, 49]))["status"].tolist() == ["confirmed"]       # 54 (practice) -> 49: 5 points
    assert events(V([50, 50, 47, 47]))["status"].tolist() == []                       # 3 points, 6 %: no
    assert events(V([30, 30, 27, 27]))["status"].tolist() == ["confirmed"]            # 3 points but 10 % of 30
    assert events(V([30, 30, 27, 27]), Rules(combine="and")).empty
    assert events(V([50, 50, 45, 50, 45]))["status"].tolist() == ["unconfirmed"]      # recovers, then open at the end


def test_first_visit_reference_hides_decline():
    v = V([46, 52, 47, 47])
    assert len(events(v)) == 1 and events(v, Rules(baseline_visits=1)).empty



def test_decline_at_second_visit():
    v = pd.DataFrame({"patient": "a", "date": ["2020-01-01", "2020-04-01", "2020-07-01"], "sdmt": [60, 50, 50]})
    assert events(v)["status"].tolist() == ["confirmed"]
