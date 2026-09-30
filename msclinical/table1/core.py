"""Descriptive table of an MS cohort, per group.

Continuous variables: median (IQR), since age, disease duration and EDSS are skewed; EDSS is ordinal and
is never summarised as a mean. Counts: n (%). ARR: total relapses over total patient-years (the cohort
rate), not the mean of per-patient rates, which small follow-ups inflate. Every row states how many
patients contributed, so missingness is visible instead of silently shrinking the denominator.
"""

from __future__ import annotations

import pandas as pd


def _med(s: pd.Series, d: int = 1) -> str:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if not len(s):
        return "–"
    q1, q2, q3 = s.quantile([0.25, 0.5, 0.75])
    return f"{q2:.{d}f} ({q1:.{d}f}–{q3:.{d}f})"


def _pct(n: int, total: int) -> str:
    return f"{n} ({100 * n / total:.0f} %)" if total else "–"


def table1(tables: dict, group: str | None = "arm", index_date: str | None = None) -> pd.DataFrame:
    p = tables["patients"].copy()
    p["start"] = pd.to_datetime(p["start"])
    p["end"] = pd.to_datetime(p["end"])
    ref = pd.Timestamp(index_date) if index_date else None
    groups = {"All": p} if (group is None or group not in p) else {"All": p, **{str(k): g for k, g in p.groupby(group)}}
    v = tables.get("visits")
    r = tables.get("relapses")
    m = tables.get("mri")
    rows = []
    for name, g in groups.items():
        col = {}
        n = len(g)
        col["Patients"] = str(n)
        at = g["start"] if ref is None else pd.Series(ref, index=g.index)
        if "birth_year" in g:
            col["Age at start, years"] = _med(at.dt.year - g["birth_year"], 0) + f"  [n={g['birth_year'].notna().sum()}]"
        if "sex" in g:
            f = g["sex"].astype(str).str.upper().str[0]
            col["Female"] = _pct(int((f == "F").sum()), int(g["sex"].notna().sum()))
        if "course" in g:
            for c, k in g["course"].value_counts().items():
                col[f"Course: {c}"] = _pct(int(k), n)
        if "onset_date" in g:
            dur = (at - pd.to_datetime(g["onset_date"])).dt.days / 365.25
            col["Disease duration, years"] = _med(dur) + f"  [n={dur.notna().sum()}]"
        fu = (g["end"] - g["start"]).dt.days / 365.25
        col["Follow-up, years"] = _med(fu)
        if v is not None and "edss" in v:
            vv = v.assign(date=pd.to_datetime(v["date"]))
            vv = vv[vv["patient"].isin(g["patient"])].merge(g[["patient", "start"]], on="patient")
            # the last EDSS on or before the start; when none exists, the first within 90 days after it
            pre = vv[vv["date"] <= vv["start"]].sort_values("date").groupby("patient")["edss"].last()
            post = vv[(vv["date"] > vv["start"]) & (vv["date"] <= vv["start"] + pd.Timedelta(days=90))].sort_values("date").groupby("patient")["edss"].first()
            base = pre.combine_first(post)
            col["Baseline EDSS"] = _med(base) + f"  [n={len(base)}]"
        if r is not None:
            rr = r.assign(date=pd.to_datetime(r["date"])).merge(g[["patient", "start", "end"]], on="patient")
            inside = rr[(rr["date"] > rr["start"]) & (rr["date"] <= rr["end"])]
            py = fu.sum()
            col["Relapses on study, ARR"] = f"{len(inside)} / {py:.0f} patient-years = {len(inside) / py:.2f}" if py else "–"
        if m is not None and "new_t2" in m:
            mm = m.assign(date=pd.to_datetime(m["date"])).merge(g[["patient", "start", "end"]], on="patient")
            mm = mm[(mm["date"] > mm["start"]) & (mm["date"] <= mm["end"])]  # on study only
            read = mm.dropna(subset=["new_t2"])
            any_new = read.groupby("patient")["new_t2"].sum() > 0
            col["Any new T2 lesion on study"] = _pct(int(any_new.sum()), len(any_new)) + f"  [n={len(any_new)} with MRI]"
        rows.append(pd.Series(col, name=name))
    return pd.DataFrame(rows).T.fillna("")
