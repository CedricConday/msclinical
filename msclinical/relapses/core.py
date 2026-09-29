"""Relapse table hygiene.

Clinical definition (McDonald criteria; Thompson et al. 2018): new or worsening symptoms lasting at least
24 hours, in the absence of fever or infection, separated from a previous attack by at least 30 days.
Checks (columns beyond patient and date are optional):
* ``merge``: onsets less than ``min_gap_days`` after the previous accepted onset of the same patient are
  folded into it (one attack reported twice);
* ``pseudo``: ``fever`` or ``infection`` true, flagged as possible pseudo-relapse;
* ``short``: ``duration_hours`` under 24;
* ``outside``: before the patient's follow-up start or after its end (``followup`` table);
* ``duplicate``: same patient and date twice.
The clean table keeps accepted relapses only; every removed row is listed with its reason.
"""

from __future__ import annotations

import pandas as pd

TRUE = {"yes", "true", "1", "y", "ja", "j"}
FALSE = {"no", "false", "0", "n", "", "nein", "nan", "none"}


def _flag(s: pd.Series, name: str) -> pd.Series:
    def one(x):
        if isinstance(x, bool):
            return x
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return False
        t = str(x).strip().lower()
        if t in TRUE:
            return True
        if t in FALSE:
            return False
        raise ValueError(f"{name}: cannot read {x!r} as yes/no")
    return s.map(one).astype(bool)


def check(relapses: pd.DataFrame, followup: pd.DataFrame | None = None, min_gap_days: int = 30,
          drop_pseudo: bool = True, drop_short: bool = True, dayfirst: bool = False,
          date_format: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    r = relapses.copy()
    r["date"] = pd.to_datetime(r["date"], dayfirst=dayfirst, format=date_format)
    r = r.sort_values(["patient", "date"], kind="stable").reset_index(drop=True)
    r["reason"] = ""
    r["pseudo_flag"] = False
    if "fever" in r or "infection" in r:
        pseudo = pd.Series(False, index=r.index)
        for col in ("fever", "infection"):
            if col in r:
                pseudo |= _flag(r[col], col)
        r["pseudo_flag"] = pseudo
        if drop_pseudo:
            r.loc[pseudo & (r["reason"] == ""), "reason"] = "pseudo (fever or infection)"
    if "duration_hours" in r and drop_short:
        r.loc[(r["duration_hours"] < 24) & (r["reason"] == ""), "reason"] = "shorter than 24 h"
    if followup is not None:
        f = followup.assign(start=pd.to_datetime(followup["start"], dayfirst=dayfirst, format=date_format),
                            end=pd.to_datetime(followup["end"], dayfirst=dayfirst, format=date_format))
        f = f.groupby("patient").agg(start=("start", "min"), end=("end", "max"))
        s = r["patient"].map(f["start"])
        e = r["patient"].map(f["end"])
        r.loc[s.isna() & (r["reason"] == ""), "reason"] = "patient not in follow-up table"
        r.loc[((r["date"] < s) | (r["date"] > e)) & (r["reason"] == ""), "reason"] = "outside follow-up"
    last = {}
    for i, row in r.iterrows():
        if row["reason"]:
            continue
        p = row["patient"]
        if p in last and row["date"] == last[p]:
            r.at[i, "reason"] = "duplicate"
        elif p in last and (row["date"] - last[p]).days < min_gap_days:
            r.at[i, "reason"] = f"within {min_gap_days} d of the attack on {last[p].date()} (merged)"
        else:
            last[p] = row["date"]
    clean = r[r["reason"] == ""].drop(columns="reason")
    removed = r[r["reason"] != ""]
    clean = clean.assign(date=clean["date"].dt.date)
    removed = removed.assign(date=removed["date"].dt.date)
    summ = {"input": len(r), "kept": len(clean), "removed": len(removed),
            "by_reason": removed["reason"].str.replace(r" on \d{4}-\d{2}-\d{2}", "", regex=True).value_counts().to_dict()}
    return clean, removed, summ
