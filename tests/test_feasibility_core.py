import pandas as pd

from msclinical.feasibility.core import evaluate


def test_funnel():
    t = {"patients": pd.DataFrame([("a", 1980, "RRMS"), ("b", 1950, "RRMS"), ("c", 1990, "SPMS"), ("d", 1985, "RRMS")],
                                  columns=["patient", "birth_year", "course"]),
         "visits": pd.DataFrame([("a", "2024-06-01", 3.0), ("b", "2024-06-01", 3.0), ("c", "2024-06-01", 6.5), ("d", "2024-06-01", 2.0)],
                                columns=["patient", "date", "edss"]),
         "relapses": pd.DataFrame([("a", "2024-01-01"), ("d", "2024-05-20")], columns=["patient", "date"]),
         "mri": pd.DataFrame([("c", "2024-03-01", 2)], columns=["patient", "date", "gd"])}
    crit = {"age": [18, 55], "course": ["RRMS"], "edss": [0, 5.5], "activity": {"relapses": 1, "gd": 1}, "no_relapse_days": 30}
    df, f = evaluate(t, crit)
    assert f["remaining"].tolist() == [4, 3, 2, 2, 2, 1] and df.set_index("patient").loc["a", "eligible"]



def test_unknowns():
    t = {"patients": pd.DataFrame([("a", None, None)], columns=["patient", "birth_year", "course"]),
         "visits": pd.DataFrame([("a", "2024-06-01", 3.0)], columns=["patient", "date", "edss"]),
         "relapses": pd.DataFrame(columns=["patient", "date"]), "mri": pd.DataFrame(columns=["patient", "date", "gd"])}
    df, f = evaluate(t, {"age": [18, 55], "course": ["RRMS"], "activity": {"relapses": 1, "gd": 1}})
    assert f["unknown"].tolist()[1:] == [1, 0, 0]          # the funnel counts unknowns among patients still in it
    assert df[["age", "course", "activity"]].isna().all(axis=None)
