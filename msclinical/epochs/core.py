"""Treatment epochs from a treatments table (patient, drug, start[, stop]).

An epoch is one uninterrupted period on one drug. Overlapping courses of the same drug merge; a gap of at
most ``gap_days`` between two courses of the same drug is bridged. An epoch's analysis window starts
``rebaseline_months`` after its start (the drug's onset of action; 3 to 6 months are common choices) and
ends at the next drug's start, the stop date plus ``carryover_days`` (washout and carry-over), or the end
of follow-up, whichever comes first. Untreated time between epochs is its own epoch with drug "none".
Every row keeps the rule that cut it, so a reviewer can see why an event belongs where it does.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise

import pandas as pd


@dataclass(frozen=True)
class Rules:
    gap_days: int = 30
    rebaseline_months: int = 6
    carryover_days: int = 90
    min_window_days: int = 90


def epochs(treatments: pd.DataFrame, followup: pd.DataFrame, r: Rules | None = None) -> pd.DataFrame:
    """Timeline segmentation: at every moment the most recently started active course holds the time, so a
    course nested inside another interrupts it and the outer drug resumes afterwards; time with no course is
    "none". A drug epoch's window starts ``rebaseline_months`` after the course start (not after a resumption),
    and when a drug is followed by untreated time its window extends by ``carryover_days`` into it."""
    r = r or Rules()
    t = treatments.assign(start=pd.to_datetime(treatments["start"]),
                          stop=pd.to_datetime(treatments["stop"]) if "stop" in treatments else pd.NaT)
    f = followup.assign(start=pd.to_datetime(followup["start"]), end=pd.to_datetime(followup["end"])).set_index("patient")
    rows = []
    for pid, fu in f.iterrows():
        g = t[t["patient"] == pid].sort_values("start")
        courses = []
        for _, c in g.iterrows():
            stop = min(c["stop"] if pd.notna(c["stop"]) else fu["end"], fu["end"])
            if courses and courses[-1]["drug"] == c["drug"] and c["start"] <= courses[-1]["stop"] + pd.Timedelta(days=r.gap_days):
                courses[-1]["stop"] = max(courses[-1]["stop"], stop)
            else:
                courses.append({"drug": c["drug"], "start": c["start"], "stop": stop})
        pts = sorted({fu["start"], fu["end"], *[c["start"] for c in courses], *[c["stop"] for c in courses]})
        pts = [x for x in pts if fu["start"] <= x <= fu["end"]]
        segs = []
        for a, b in pairwise(pts):
            active = [c for c in courses if c["start"] <= a < c["stop"]]
            cur = max(active, key=lambda c: c["start"]) if active else None
            key = (cur["drug"], cur["start"]) if cur else ("none", None)
            if segs and segs[-1]["key"] == key and segs[-1]["end"] == a:
                segs[-1]["end"] = b
            else:
                segs.append({"key": key, "start": a, "end": b, "course": cur})
        for i, sg in enumerate(segs):
            drug = sg["key"][0]
            nxt = segs[i + 1] if i + 1 < len(segs) else None
            if drug == "none":
                prev = segs[i - 1] if i else None
                ws = sg["start"] + (pd.Timedelta(days=r.carryover_days) if prev is not None else pd.Timedelta(0))
                rows.append({"patient": pid, "drug": "none", "start": sg["start"], "stop": sg["end"], "window_start": min(ws, sg["end"]),
                             "window_end": sg["end"], "cut": "untreated" + (" (after carry-over)" if prev is not None else "")})
                continue
            c = sg["course"]
            ws = sg["start"] + pd.DateOffset(months=r.rebaseline_months) if sg["start"] == c["start"] else sg["start"]
            if nxt is None:
                we, why = sg["end"], "end of follow-up"
            elif nxt["key"][0] == "none":
                we, why = min(sg["end"] + pd.Timedelta(days=r.carryover_days), nxt["end"]), "stop + carry-over"
            else:
                we, why = sg["end"], "next drug" if nxt["course"]["start"] == nxt["start"] else f"{nxt['key'][0]} resumes"
            rows.append({"patient": pid, "drug": drug, "start": sg["start"], "stop": sg["end"], "window_start": ws,
                         "window_end": we, "cut": why})
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    df["window_days"] = (df["window_end"] - df["window_start"]).dt.days.clip(lower=0)
    df["analysable"] = df["window_days"] >= r.min_window_days
    for c in ("start", "stop", "window_start", "window_end"):
        df[c] = df[c].dt.date
    return df


def slice_table(table: pd.DataFrame, ep: pd.DataFrame, date_col: str = "date") -> pd.DataFrame:
    """Attach each row of a dated table (relapses, visits, MRI) to the epoch whose analysis window holds it."""
    t = table.assign(**{date_col: pd.to_datetime(table[date_col])})
    e = ep.assign(ws=pd.to_datetime(ep["window_start"]), we=pd.to_datetime(ep["window_end"])).reset_index(names="epoch")
    m = t.merge(e[["patient", "epoch", "drug", "ws", "we"]], on="patient")
    m = m[(m[date_col] >= m["ws"]) & (m[date_col] < m["we"])]
    return m.drop(columns=["ws", "we"])


def arr_by_drug(relapses: pd.DataFrame, ep: pd.DataFrame) -> pd.DataFrame:
    e = ep[ep["analysable"]]
    s = slice_table(relapses, e)
    years = e.groupby("drug")["window_days"].sum() / 365.25
    n = s.groupby("drug").size().reindex(years.index, fill_value=0)
    return pd.DataFrame({"patient_years": years.round(2), "relapses": n, "arr": (n / years).round(3)})
