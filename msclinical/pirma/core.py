"""Attribute confirmed worsening events to relapses and MRI activity.

For each confirmed event (from ``edssplus``) with onset date o, reference date r and confirmation date c:

* **RAW** if a relapse began within ``raw_days`` (default 90) before o, or between o and c.
  Relapses between o and c also make the event RAW rather than PIRA: the confirmation was not
  relapse-free (Kappos et al. 2020; Müller et al. 2023 harmonisation).
* **PIRA** otherwise, when relapse ascertainment covers [o - raw_days, c] (the window a relapse must be
  absent from; relapses between the reference and onset - raw_days do not change the class).
* **PIRMA** is PIRA with no MRI activity (new or enlarging T2 lesions, or gadolinium-enhancing
  lesions) on any scan dated in [r, c] (``mri_window``: ``reference`` or ``onset``), and at least one
  scan in that window. PIRA with activity is ``PIRA_MRI_active``; PIRA with no scan in the window is
  ``PIRA_MRI_unknown`` (never silently PIRMA).
* An event whose window starts before ``followup_start`` (default: the patient's first visit) has no
  relapse coverage and is ``undetermined``.
* MRI scans count only when every activity column present is read (not missing); a scan row with missing
  counts is not evidence of no activity. An MRI table with none of the activity columns is an error.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Windows:
    raw_days: int = 90
    mri_window: str = "reference"  # or "onset"
    mri_columns: tuple = ("new_t2", "enlarging_t2", "gd")


def attribute(events: pd.DataFrame, relapses: pd.DataFrame, mri: pd.DataFrame | None = None,
              followup_start: pd.DataFrame | None = None, w: Windows | None = None,
              visits_first: pd.Series | None = None) -> pd.DataFrame:
    """events: edssplus.events output. relapses: patient, date. mri: patient, date + mri_columns.
    followup_start: patient, date (start of relapse ascertainment); default is each patient's reference date."""
    w = w or Windows()
    ev = events[events["status"] == "confirmed"].copy()
    for c in ("onset_date", "reference_date", "confirmation_date"):
        ev[c] = pd.to_datetime(ev[c])
    rel = relapses.assign(date=pd.to_datetime(relapses["date"]))
    mr = mri.assign(date=pd.to_datetime(mri["date"])) if mri is not None else None
    if mr is not None and len(mr) and not [c for c in w.mri_columns if c in mr]:
        raise ValueError(f"MRI table has none of the activity columns {w.mri_columns}; rename them or set Windows.mri_columns")
    if followup_start is not None:
        fs = followup_start.assign(date=pd.to_datetime(followup_start["date"])).groupby("patient")["date"].min()
    else:
        fs = None
    first_visit = pd.to_datetime(events["reference_date"]).groupby(events["patient"]).min() if len(events) else pd.Series(dtype="datetime64[ns]")
    if visits_first is not None:
        first_visit = visits_first
    out = []
    for _, e in ev.iterrows():
        pid, o, r, c = e["patient"], e["onset_date"], e["reference_date"], e["confirmation_date"]
        pr = rel[rel["patient"] == pid]["date"]
        start = fs.get(pid, first_visit.get(pid, r)) if fs is not None else first_visit.get(pid, r)
        near = pr[(pr >= o - pd.Timedelta(days=w.raw_days)) & (pr <= c)]
        if len(near):
            cls, why = "RAW", f"relapse {near.min().date()}"
        elif start > o - pd.Timedelta(days=w.raw_days):
            cls, why = "undetermined", f"relapse ascertainment starts {start.date()}, after onset - {w.raw_days} d"
        else:
            cls, why = "PIRA", f"no relapse from onset - {w.raw_days} d to confirmation"
            if mr is None:
                cls, why = "PIRA_MRI_unknown", why + "; no MRI table"
            else:
                lo = r if w.mri_window == "reference" else o
                cols = [x for x in w.mri_columns if x in mr]
                scans = mr[(mr["patient"] == pid) & (mr["date"] >= lo) & (mr["date"] <= c)]
                scans = scans[scans[cols].notna().all(axis=1)]  # unread counts are not "no activity"
                if not len(scans):
                    cls, why = "PIRA_MRI_unknown", why + f"; no scan {lo.date()} to {c.date()}"
                elif (scans[cols].fillna(0).sum(axis=1) > 0).any():
                    cls, why = "PIRA_MRI_active", why + "; MRI activity " + ", ".join(
                        str(d.date()) for d in scans[scans[cols].fillna(0).sum(axis=1) > 0]["date"])
                else:
                    cls, why = "PIRMA", why + f"; {len(scans)} scan(s) without activity"
        out.append({**e.to_dict(), "class": cls, "reason": why})
    df = pd.DataFrame(out)
    for c in ("onset_date", "reference_date", "confirmation_date"):
        if c in df:
            df[c] = df[c].dt.date
    return df


def summary(att: pd.DataFrame) -> dict:
    counts = att["class"].value_counts().to_dict() if len(att) else {}
    n = len(att)
    pira_all = sum(counts.get(k, 0) for k in ("PIRA", "PIRMA", "PIRA_MRI_active", "PIRA_MRI_unknown"))
    return {"events": n, "counts": counts, "pira_share": pira_all / n if n else None,
            "pirma_share_of_pira_with_mri": (counts.get("PIRMA", 0) / (counts.get("PIRMA", 0) + counts.get("PIRA_MRI_active", 0))
                                             if counts.get("PIRMA", 0) + counts.get("PIRA_MRI_active", 0) else None)}
