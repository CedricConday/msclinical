"""Confirmed SDMT worsening.

A decline of at least ``points`` (default 4) raw SDMT points, or at least ``pct`` (default 10 %) of the
reference, is the change most often treated as clinically meaningful (Benedict et al., Mult Scler 2017;
the rule is configurable because the literature uses both, and "or" versus "and" matters). The decline must
hold at every SDMT visit from onset to the first visit at least ``confirm_days`` later. Practice effects:
within the first ``baseline_visits`` administrations (default 2) a higher score replaces the reference, because
the second administration is typically higher than the first and a first-visit reference makes later declines
look smaller than they are; an early visit that is lower is tested as an onset like any other. Roving reference as in edssplus: after a confirmed event the reference moves to
the confirmed value.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Rules:
    points: float = 4.0
    pct: float = 10.0
    combine: str = "or"  # "or": either criterion; "and": both
    confirm_days: int = 84
    baseline_visits: int = 2
    roving: bool = True


def declined(ref: float, val: float, r: Rules) -> bool:
    a = ref - val >= r.points - 1e-9
    b = ref - val >= ref * r.pct / 100 - 1e-9
    return (a or b) if r.combine == "or" else (a and b)


def events(visits: pd.DataFrame, r: Rules | None = None) -> pd.DataFrame:
    r = r or Rules()
    rows = []
    for pid, g in visits.dropna(subset=["sdmt"]).groupby("patient"):
        v = g.assign(date=pd.to_datetime(g["date"])).sort_values("date").reset_index(drop=True)
        if len(v) < 2:
            continue
        ref_i, ref = 0, float(v.loc[0, "sdmt"])
        i = 1
        while i < len(v):
            val = float(v.loc[i, "sdmt"])
            if i < r.baseline_visits and val > ref:
                ref_i, ref = i, val  # practice effect: a higher early score becomes the reference
                i += 1
                continue
            if declined(ref, val, r):
                onset = v.loc[i, "date"]
                j, held, conf = i, True, None
                while j < len(v):
                    if not declined(ref, float(v.loc[j, "sdmt"]), r):
                        held = False
                        break
                    if (v.loc[j, "date"] - onset).days >= r.confirm_days:
                        conf = j
                        break
                    j += 1
                status = "confirmed" if (held and conf is not None) else ("unconfirmed" if held else None)
                if status:
                    rows.append({"patient": pid, "status": status, "reference_date": v.loc[ref_i, "date"].date(), "reference": ref,
                                 "onset_date": onset.date(), "onset_value": val,
                                 "confirmation_date": v.loc[conf, "date"].date() if conf is not None else None,
                                 "confirmed_value": float(v.loc[conf, "sdmt"]) if conf is not None else None})
                if status == "confirmed":
                    if not r.roving:
                        break
                    ref_i, ref = conf, float(v.loc[conf, "sdmt"])
                    i = conf + 1
                    continue
                if status == "unconfirmed":
                    break
            i += 1
    return pd.DataFrame(rows, columns=["patient", "status", "reference_date", "reference", "onset_date", "onset_value",
                                       "confirmation_date", "confirmed_value"])
