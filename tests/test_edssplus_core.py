import pandas as pd

from msclinical.edssplus.core import Rules, edss_plus, edss_threshold, events


def V(rows):
    return pd.DataFrame(rows, columns=["patient", "date", "edss", "t25fw"])


def test_thresholds():
    assert edss_threshold(0) == 1.5 and edss_threshold(3.0) == 1.0 and edss_threshold(5.5) == 1.0 and edss_threshold(6.0) == 0.5


def test_confirmed_edss_event():
    v = V([("a", "2020-01-01", 2.0, 5), ("a", "2020-06-01", 3.0, 5), ("a", "2020-09-01", 3.0, 5)])
    e = events(v)
    assert len(e) == 1 and e.iloc[0]["status"] == "confirmed" and str(e.iloc[0]["confirmation_date"]) == "2020-09-01"


def test_not_sustained_is_no_event():
    v = V([("a", "2020-01-01", 2.0, 5), ("a", "2020-06-01", 3.0, 5), ("a", "2020-07-01", 2.0, 5), ("a", "2020-12-01", 2.0, 5)])
    assert events(v).empty


def test_unconfirmed_at_end():
    v = V([("a", "2020-01-01", 2.0, 5), ("a", "2020-06-01", 3.0, 5)])
    assert events(v).iloc[0]["status"] == "unconfirmed"


def test_roving_allows_second_event_fixed_does_not_double_count():
    rows = [("a", "2020-01-01", 1.0, 5), ("a", "2020-04-01", 2.0, 5), ("a", "2020-08-01", 2.0, 5),
            ("a", "2021-01-01", 3.0, 5), ("a", "2021-06-01", 3.0, 5)]
    assert (events(V(rows))["status"] == "confirmed").sum() == 2
    assert (events(V(rows), Rules(reference="fixed"))["status"] == "confirmed").sum() == 1


def test_t25fw_20pct_and_composite():
    v = V([("a", "2020-01-01", 2.0, 5.0), ("a", "2020-06-01", 2.0, 6.1), ("a", "2020-10-01", 2.0, 6.2),
           ("b", "2020-01-01", 2.0, 5.0), ("b", "2021-01-01", 2.0, 5.0)])
    c = edss_plus(v)
    assert c.set_index("patient").loc["a", "component"] == "t25fw" and not c.set_index("patient").loc["b", "event"]


def test_24_week_needs_later_visit():
    v = V([("a", "2020-01-01", 2.0, 5), ("a", "2020-06-01", 3.0, 5), ("a", "2020-09-01", 3.0, 5)])
    assert events(v, Rules(confirm_days=168)).iloc[0]["status"] == "unconfirmed"


def test_non_iso_dates_sorted_after_parsing():
    v = V([("a", "12/01/2019", 3.0, 5), ("a", "01/15/2020", 2.0, 5), ("a", "02/10/2020", 3.0, 5), ("a", "06/01/2020", 3.0, 5)])
    assert events(v).empty  # 3.0 -> 2.0 -> 3.0: no confirmed worsening from the first reference


def test_roving_improvement_needs_every_visit():
    v = V([("a", "2020-01-01", 4.0, 5), ("a", "2020-02-01", 3.0, 5), ("a", "2020-03-01", 4.5, 5), ("a", "2020-05-01", 3.0, 5),
           ("a", "2020-06-01", 4.0, 5), ("a", "2020-09-01", 4.0, 5)])
    assert (events(v)["status"] == "confirmed").sum() == 0


def test_custom_steps_are_used():
    r = Rules(edss_steps=((0.5, 1.5), (5.5, 1.0), (10.0, 0.5)))
    assert edss_threshold(0.5, r) == 1.5 and edss_threshold(1.0, r) == 1.0
