"""The definitional multiverse of disability worsening in MS.

Choices (each seen in published cohorts; see Bsteh et al. 2025/2026 for how far they spread):
confirmation 12 or 24 weeks; roving or fixed reference; relapse window 30 or 90 days; MRI window from
the reference or from onset; EDSS alone or the EDSS-plus composite. For every combination: confirmed
events, PIRA share, PIRMA share. Then, per choice, the mean absolute change in PIRA share when only that
choice flips: which decision moves the result most.
"""

from __future__ import annotations

from itertools import product

import pandas as pd

from msclinical.edssplus.core import Rules, events
from msclinical.pirma.core import Windows, attribute, summary

GRID = {"confirm_days": (84, 168), "reference": ("roving", "fixed"), "raw_days": (30, 90),
        "mri_window": ("reference", "onset"), "components": ("edss", "edss_plus")}


def run(visits, relapses, mri=None, followup_start=None, grid=None) -> pd.DataFrame:
    grid = grid or GRID
    rows = []
    keys = list(grid)
    cache = {}
    for combo in product(*(grid[k] for k in keys)):
        c = dict(zip(keys, combo))
        ek = (c["confirm_days"], c["reference"], c["components"])
        if ek not in cache:
            comps = ("edss",) if c["components"] == "edss" else ("edss", "t25fw", "hpt9_dom", "hpt9_nondom")
            ev = events(visits, Rules(confirm_days=c["confirm_days"], reference=c["reference"]), components=comps)
            cache[ek] = ev
        att = attribute(cache[ek], relapses, mri, followup_start, Windows(raw_days=c["raw_days"], mri_window=c["mri_window"]))
        s = summary(att) if len(att) else {"events": 0, "counts": {}, "pira_share": None}
        cnt = s["counts"]
        rows.append({**c, "events": s["events"], "patients_with_event": int(att["patient"].nunique()) if len(att) else 0,
                     "pira_share": s["pira_share"], "pirma": cnt.get("PIRMA", 0), "raw": cnt.get("RAW", 0)})
    return pd.DataFrame(rows)


def influence(table: pd.DataFrame, metric: str = "pira_share") -> pd.DataFrame:
    keys = [k for k in GRID if k in table]
    out = []
    for k in keys:
        others = [x for x in keys if x != k]
        piv = table.pivot_table(index=others, columns=k, values=metric)
        if piv.shape[1] == 2:
            d = (piv.iloc[:, 1] - piv.iloc[:, 0]).abs()
            out.append({"choice": k, "mean_abs_change": float(d.mean()), "max_abs_change": float(d.max())})
    return pd.DataFrame(out).sort_values("mean_abs_change", ascending=False)
