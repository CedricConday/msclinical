"""Confirmed disability worsening (CDW) events, one component at a time or as the EDSS-plus composite.

Rules (all in ``Rules``):
* EDSS worsening from reference r: +1.5 if r == 0, +1.0 if 0 < r <= 5.5, +0.5 if r > 5.5
  (the thresholds used in most phase 3 programmes and in Kappos et al. 2018/2020).
* Timed 25-Foot Walk and 9-Hole Peg Test (either hand): >= 20 % slower than reference
  (EDSS-plus, Cadavid et al. 2017).
* Confirmation: the worsening must hold at every visit of that component from onset to the first
  visit at least ``confirm_days`` after onset, which is the confirmation visit. A worsening that
  has no such visit (end of follow-up) is reported as ``unconfirmed``, never as an event.
* Reference: ``fixed`` keeps the first visit and reports the first event only; ``roving`` moves the reference to the confirmed
  value after each event, and (``rove_on_improvement``) to a lower value once that lower value is
  itself confirmed by the next visit at least ``confirm_days`` later (Kappos et al. 2018).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

COMPONENTS = ("edss", "t25fw", "hpt9_dom", "hpt9_nondom")


@dataclass(frozen=True)
class Rules:
    confirm_days: int = 84  # 12 weeks; 168 for 24-week CDW
    reference: str = "roving"  # or "fixed"
    rove_on_improvement: bool = True
    timed_pct: float = 20.0
    edss_steps: tuple = ((0.0, 1.5), (5.5, 1.0), (10.0, 0.5))  # (reference upper bound, required increase)
    dayfirst: bool = False  # set True for 05.01.2020 = 5 January
    date_format: str | None = None  # e.g. "%d.%m.%Y"; None lets pandas infer (ISO recommended)


def edss_threshold(ref: float, rules: Rules | None = None) -> float:
    rules = rules or Rules()
    for upper, inc in rules.edss_steps:
        if ref <= upper:
            return inc
    return rules.edss_steps[-1][1]


def worse(component: str, ref: float, value: float, rules: Rules | None = None) -> bool:
    rules = rules or Rules()
    if pd.isna(ref) or pd.isna(value):
        return False
    if component == "edss":
        return value - ref >= edss_threshold(ref, rules) - 1e-9
    return value >= ref * (1 + rules.timed_pct / 100.0) - 1e-9


def better(component: str, ref: float, value: float) -> bool:
    return (not pd.isna(value)) and value < ref


def component_events(visits: pd.DataFrame, component: str, rules: Rules | None = None) -> list[dict]:
    """Events for one patient and one component. ``visits`` has columns date and ``component``."""
    rules = rules or Rules()
    v = visits[["date", component]].dropna().copy()
    v["date"] = pd.to_datetime(v["date"], dayfirst=rules.dayfirst, format=rules.date_format)
    v = v.sort_values("date").reset_index(drop=True)
    out: list[dict] = []
    if len(v) < 2:
        return out
    ref_i = 0
    ref = float(v.loc[0, component])
    i = 1
    while i < len(v):
        val = float(v.loc[i, component])
        if worse(component, ref, val, rules):
            onset = v.loc[i, "date"]
            j, held, conf_i = i, True, None
            while j < len(v):
                if not worse(component, ref, float(v.loc[j, component]), rules):
                    held = False
                    break
                if (v.loc[j, "date"] - onset).days >= rules.confirm_days:
                    conf_i = j
                    break
                j += 1
            if held and conf_i is not None:
                out.append({"component": component, "status": "confirmed",
                            "reference_date": v.loc[ref_i, "date"].date(), "reference_value": ref,
                            "onset_date": onset.date(), "onset_value": val,
                            "confirmation_date": v.loc[conf_i, "date"].date(),
                            "confirmed_value": float(v.loc[conf_i, component])})
                if rules.reference != "roving":
                    break  # a fixed reference defines time to first event; later ones would re-count it
                ref_i, ref = conf_i, float(v.loc[conf_i, component])
                i = conf_i + 1
                continue
            if held and conf_i is None:
                out.append({"component": component, "status": "unconfirmed",
                            "reference_date": v.loc[ref_i, "date"].date(), "reference_value": ref,
                            "onset_date": onset.date(), "onset_value": val,
                            "confirmation_date": None, "confirmed_value": None})
                break
        elif rules.reference == "roving" and rules.rove_on_improvement and better(component, ref, val):
            # a lower value becomes the reference once a visit >= confirm_days later is no higher
            # the lower value becomes the reference only when every visit up to and including the first one
            # >= confirm_days later stays at or below it (the same standard the worsening side meets)
            gap = (v["date"] - v.loc[i, "date"]).dt.days
            conf = v.index[(gap >= rules.confirm_days)]
            if len(conf):
                j = conf[0]
                if (v.loc[i:j, component].astype(float) <= val).all():
                    ref_i, ref = i, val
        i += 1
    return out


def events(visits: pd.DataFrame, rules: Rules | None = None, components: tuple = COMPONENTS) -> pd.DataFrame:
    """All components for all patients. ``visits`` columns: patient, date, and any of COMPONENTS."""
    rules = rules or Rules()
    rows = []
    for pid, g in visits.groupby("patient", sort=True):
        for c in components:
            if c in g:
                for e in component_events(g, c, rules):
                    rows.append({"patient": pid, **e})
    cols = ["patient", "component", "status", "reference_date", "reference_value", "onset_date", "onset_value",
            "confirmation_date", "confirmed_value"]
    return pd.DataFrame(rows, columns=cols)


def edss_plus(visits: pd.DataFrame, rules: Rules | None = None) -> pd.DataFrame:
    """First confirmed event per patient on any component (the EDSS-plus composite), with its trigger."""
    ev = events(visits, rules)
    conf = ev[ev["status"] == "confirmed"].sort_values(["patient", "onset_date"])
    first = conf.groupby("patient", as_index=False).first()
    pats = pd.DataFrame({"patient": sorted(visits["patient"].unique())})
    out = pats.merge(first, on="patient", how="left")
    out["event"] = out["component"].notna()
    last = visits.assign(date=pd.to_datetime(visits["date"])).groupby("patient")["date"].max().dt.date
    out["censor_date"] = out["patient"].map(last)
    return out
