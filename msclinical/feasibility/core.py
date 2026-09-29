"""Trial feasibility on registry tables (mstables schema).

Criteria are a small dict (YAML or JSON on the command line). Each is evaluated per patient at an index
date (default: the last visit), and the funnel reports how many remain after each criterion in the given
order, so the binding criterion is visible. Supported keys:

    age: [min, max]                      needs patients.birth_year
    edss: [min, max]                     latest EDSS in the ``edss_lookback_days`` before the index
    disease_years: [min, max]            needs patients.onset_date
    course: [RRMS, SPMS]                 patients.course
    activity: {relapses: n, relapse_days: 365, gd: n, gd_days: 365, any: true}
                                         relapses and/or enhancing lesions in the look-back; any = either suffices
    no_relapse_days: 30                  no relapse in the last N days (stability before screening)
"""

from __future__ import annotations

import pandas as pd


def evaluate(tables: dict, criteria: dict, index_date=None, edss_lookback_days: int = 365) -> tuple[pd.DataFrame, pd.DataFrame]:
    p = tables["patients"].copy()
    v = tables.get("visits", pd.DataFrame(columns=["patient", "date", "edss"])).assign(date=lambda d: pd.to_datetime(d["date"]))
    r = tables.get("relapses", pd.DataFrame(columns=["patient", "date"])).assign(date=lambda d: pd.to_datetime(d["date"]))
    m = tables.get("mri", pd.DataFrame(columns=["patient", "date", "gd"])).assign(date=lambda d: pd.to_datetime(d["date"]))
    idx = {pid: pd.Timestamp(index_date) if index_date else (g["date"].max() if len(g) else pd.NaT) for pid, g in v.groupby("patient")}
    rows = []
    for _, pt in p.iterrows():
        pid = pt["patient"]
        t = idx.get(pid, pd.Timestamp(index_date) if index_date else pd.NaT)
        rec = {"patient": pid, "index_date": t.date() if pd.notna(t) else None}
        for key, rule in criteria.items():
            ok = None
            if pd.isna(t):
                ok = None
            elif key == "age" and "birth_year" in p:
                ok = None if pd.isna(pt["birth_year"]) else bool(rule[0] <= t.year - pt["birth_year"] <= rule[1])
            elif key == "edss":
                e = v[(v["patient"] == pid) & (v["date"] <= t) & (v["date"] >= t - pd.Timedelta(days=edss_lookback_days))].dropna(subset=["edss"])
                ok = None if not len(e) else bool(rule[0] <= float(e.sort_values("date")["edss"].iloc[-1]) <= rule[1])
            elif key == "disease_years" and "onset_date" in p and pd.notna(pt.get("onset_date")):
                y = (t - pd.Timestamp(pt["onset_date"])).days / 365.25
                ok = rule[0] <= y <= rule[1]
            elif key == "course" and "course" in p:
                ok = None if pd.isna(pt["course"]) else str(pt["course"]) in set(rule)
            elif key == "activity":
                nr = ((r["patient"] == pid) & (r["date"] <= t) & (r["date"] > t - pd.Timedelta(days=rule.get("relapse_days", 365)))).sum()
                mm = m[(m["patient"] == pid) & (m["date"] <= t) & (m["date"] > t - pd.Timedelta(days=rule.get("gd_days", 365)))]
                read = mm.dropna(subset=["gd"]) if "gd" in mm else mm.iloc[0:0]
                a = bool(nr >= rule.get("relapses", 1))
                b = None if not len(read) else bool(read["gd"].sum() >= rule.get("gd", 1))  # no scan: unknown, not zero
                if rule.get("any", True):
                    ok = True if (a or b) else (None if b is None else False)
                else:
                    ok = False if (not a or b is False) else (None if b is None else True)
            elif key == "no_relapse_days":
                ok = not ((r["patient"] == pid) & (r["date"] <= t) & (r["date"] > t - pd.Timedelta(days=rule))).any()
            rec[key] = ok
        rows.append(rec)
    df = pd.DataFrame(rows)
    alive = pd.Series(True, index=df.index)
    funnel = [{"criterion": "all patients", "remaining": len(df), "unknown": 0}]
    for key in criteria:
        unknown = int((df[key].isna() & alive).sum())
        alive &= df[key].fillna(False).astype(bool)
        funnel.append({"criterion": key, "remaining": int(alive.sum()), "unknown": unknown})
    df["eligible"] = alive
    return df, pd.DataFrame(funnel)
