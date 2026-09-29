"""No evidence of disease activity (NEDA), per patient and window.

NEDA-3 (Havrdova et al. 2009): no relapse, no confirmed disability worsening, no MRI activity
(new or enlarging T2, or gadolinium-enhancing lesions). NEDA-4 (Kappos et al. 2016) adds annualised
brain volume loss no worse than 0.4 %/yr. Each component is True (no activity), False (activity) or
None (not assessable: no MRI in the window, no brain-volume data). A window is NEDA only when every
component is True; any False makes it evident activity; otherwise the status is None, never True.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from msclinical.edssplus.core import Rules, events


@dataclass(frozen=True)
class Params:
    window_months: int = 12
    rebaseline_months: int = 0  # ignore activity in the first N months after start (e.g. 6 after a new therapy)
    atrophy_pct_per_year: float = 0.4
    cdw: Rules = field(default_factory=Rules)
    mri_columns: tuple = ("new_t2", "enlarging_t2", "gd")


def _windows(start: pd.Timestamp, end: pd.Timestamp, p: Params) -> list[tuple]:
    s = start + pd.DateOffset(months=p.rebaseline_months)
    out = []
    while s < end:
        e = s + pd.DateOffset(months=p.window_months)
        out.append((s, min(e, end), e <= end))
        s = e
    return out


def neda(visits: pd.DataFrame, relapses: pd.DataFrame, mri: pd.DataFrame, starts: pd.DataFrame | None = None,
         p: Params | None = None) -> pd.DataFrame:
    p = p or Params()
    vis = visits.assign(date=pd.to_datetime(visits["date"]))
    rel = relapses.assign(date=pd.to_datetime(relapses["date"]))
    mr = mri.assign(date=pd.to_datetime(mri["date"]))
    st = starts.assign(date=pd.to_datetime(starts["date"])).groupby("patient")["date"].min() if starts is not None else None
    cols = [x for x in p.mri_columns if x in mr]
    if len(mr) and not cols:
        raise ValueError(f"MRI table has none of the activity columns {p.mri_columns}")
    rows = []
    for pid, g in vis.groupby("patient"):
        start = st.get(pid, g["date"].min()) if st is not None else g["date"].min()
        base = start + pd.DateOffset(months=p.rebaseline_months)
        # worsening is measured from the (re)baseline: the last visit on or before it is the reference
        pre = g[g["date"] <= base]
        g_cdw = g[g["date"] >= (pre["date"].max() if len(pre) else base)]
        ev = events(g_cdw.assign(date=g_cdw["date"].dt.date), p.cdw, components=("edss",))
        ev = ev.assign(onset_date=pd.to_datetime(ev["onset_date"]))
        conf = ev[ev["status"] == "confirmed"]
        unconf = ev[ev["status"] == "unconfirmed"]
        pm = mr[mr["patient"] == pid]
        pr = rel[rel["patient"] == pid]
        end = max([g["date"].max()] + ([pm["date"].max()] if len(pm) else []) + ([pr["date"].max()] if len(pr) else []))
        wins = _windows(start, end, p)
        for k, (s, e, full) in enumerate(wins, start=1):
            last = k == len(wins)
            r = pr[(pr["date"] >= s) & ((pr["date"] <= e) if last else (pr["date"] < e))]
            c = conf[(conf["onset_date"] >= s) & (conf["onset_date"] < e)]
            u = unconf[(unconf["onset_date"] >= s) & (unconf["onset_date"] < e)]
            has_visit = ((g["date"] >= s) & (g["date"] <= e)).any()
            cdw_ok = False if len(c) else (None if (len(u) or not has_visit) else True)
            m = pm[(pm["date"] > s) & (pm["date"] <= e)]
            m = m[m[cols].notna().all(axis=1)] if cols else m.iloc[0:0]  # unread counts are not "no activity"
            mri_ok = None if not len(m) else bool((m[cols].sum(axis=1) == 0).all())
            atro, rate = None, None
            if "pbvc" in m and m["pbvc"].notna().any():
                allp = mr[mr["patient"] == pid].sort_values("date").reset_index(drop=True)
                allp["prev"] = allp["date"].shift()
                w = allp[(allp["date"] > s) & (allp["date"] <= e) & allp["pbvc"].notna() & allp["prev"].notna()]
                yrs = ((w["date"] - w["prev"]).dt.days / 365.25).sum()
                if yrs > 0:
                    rate = float(w["pbvc"].sum() / yrs)
                    atro = rate >= -p.atrophy_pct_per_year
            comp3 = [len(r) == 0, cdw_ok, mri_ok]
            n3 = False if False in comp3 else (None if None in comp3 else True)
            n4 = False if (n3 is False or atro is False) else (None if (n3 is None or atro is None) else True)
            reasons = [x for x, bad in (("relapse", len(r) > 0), ("CDW", len(c) > 0), ("MRI activity", mri_ok is False),
                                        ("brain volume loss", atro is False)) if bad]
            rows.append({"patient": pid, "window": k, "start": s.date(), "end": e.date(), "complete_window": full,
                         "relapses": len(r), "cdw_events": len(c), "mri_scans": len(m), "no_mri_activity": mri_ok,
                         "pbvc_per_year": rate, "neda3": n3, "neda4": n4, "activity": ", ".join(reasons)})
    return pd.DataFrame(rows)


def cohort(table: pd.DataFrame, window: int = 1) -> dict:
    t = table[(table["window"] <= window)]
    out = {}
    for k in ("neda3", "neda4"):
        def status(g, k=k):
            if (g[k] == False).any():
                return False
            complete = g[g["complete_window"].astype(bool)]
            if len(complete) < window or complete[k].isna().any():
                return None  # follow-up shorter than the requested period, or a window not assessable
            return True
        per = pd.Series({pid: status(g) for pid, g in t.groupby("patient")}, dtype=object)
        known = per.dropna()
        out[k] = {"patients": len(per), "assessable": len(known), "neda": int(known.sum()),
                  "proportion": float(known.mean()) if len(known) else None}
    return out
