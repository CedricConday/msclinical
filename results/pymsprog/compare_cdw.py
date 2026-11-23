"""Track A: msclinical vs pymsprog on the same simulated visit tables, settings matched as far as the two APIs allow.

Each setting runs both tools and matches confirmed worsening (CDW) events on (patient, onset date).
Output: agreement.tsv (one row per setting) and disagreements_<setting>.tsv next to this script.

Usage: pip install pymsprog==1.0.7; python results/pymsprog/compare_cdw.py [--n 500] [--seeds 0 1 2]
"""

from __future__ import annotations

import argparse
import contextlib
import io
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pymsprog import MSprog

from msclinical.edssplus.core import Rules, edss_threshold, events
from msclinical.sim.core import Model, simulate

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent
INF = float("inf")


def ms_delta(ref: float) -> float:  # msclinical's EDSS threshold handed to pymsprog, to isolate everything else
    return edss_threshold(ref, Rules())


# name -> (msclinical Rules, pymsprog kwargs). Relapse exclusions off in pymsprog: msclinical's CDW ignores relapses.
NO_RELAPSE_RULES = {"relapse_to_bl": 0, "relapse_to_event": 0, "relapse_to_conf": 0}
SETTINGS = {
    # 1. first CDW, fixed baseline, 12 weeks, confirmation at the first visit >= 84 d, all visits in between worse
    "fixed_first_matched": (Rules(reference="fixed"),
                            dict(event="firstCDW", baseline="fixed", conf_days=84, conf_tol_days=[0, INF],
                                 delta_fun=ms_delta, **NO_RELAPSE_RULES)),
    # 2. same, but each tool keeps its own EDSS threshold table (differs only at reference 5.5)
    "fixed_first_own_delta": (Rules(reference="fixed"),
                              dict(event="firstCDW", baseline="fixed", conf_days=84, conf_tol_days=[0, INF],
                                   **NO_RELAPSE_RULES)),
    # 3. pymsprog at its real defaults, relapses passed in: confirmation window 84 d -7/+730, no baseline or
    #    confirmation visit within 30 d after a relapse (relapse_to_event is 0 by default); own EDSS table.
    #    msclinical with a fixed reference (its own default is roving) so the event definition is comparable.
    "fixed_first_pymsprog_defaults": (Rules(reference="fixed"),
                                      {"event": "firstCDW", "baseline": "fixed", "_relapses": True}),
    # 4. all CDW events, roving reference: msclinical re-baselines after each CDW and on a confirmed lower value
    "roving_multiple_matched": (Rules(reference="roving"),
                                dict(event="multiple", baseline="roving", proceed_from="firstconf",
                                     sub_threshold_rebl="improvement", conf_days=84, conf_tol_days=[0, INF],
                                     delta_fun=ms_delta, **NO_RELAPSE_RULES)),
}


def run_pymsprog(visits: pd.DataFrame, kw: dict, relapses: pd.DataFrame | None = None) -> pd.DataFrame:
    kw = dict(kw)
    if kw.pop("_relapses", False):
        kw.update(relapse=relapses, rsubj_col="patient", rdate_col="date")
    with contextlib.redirect_stdout(io.StringIO()):
        _, res = MSprog(visits, subj_col="patient", value_col="edss", date_col="date", outcome="edss",
                        include_dates=True, include_value=True, include_stable=False, verbose=0, **kw)
    res = res[res["event_type"].astype(str).str.contains("CDW|PIRA|RAW", regex=True)] if len(res) else res
    return pd.DataFrame({"patient": res.get("patient", []), "onset": pd.to_datetime(res.get("date", pd.Series(dtype=str))).dt.date,
                         "py_ref": res.get("bl_value", []), "py_value": res.get("value", []),
                         "py_type": res.get("event_type", [])})


def run_msclinical(visits: pd.DataFrame, rules: Rules) -> pd.DataFrame:
    ev = events(visits, rules, components=("edss",))
    ev = ev[ev["status"] == "confirmed"]
    return pd.DataFrame({"patient": ev["patient"], "onset": pd.to_datetime(ev["onset_date"]).dt.date,
                         "ms_ref": ev["reference_value"], "ms_value": ev["onset_value"]})


def perturb(visits: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Irregular visits, so confirmation tolerance and intermediate-visit checks are exercised: every date but the
    first moved by up to +-20 days, 10 % of later visits dropped, and an interim visit (carrying the previous
    value forward) added halfway into 40 % of the gaps. The simulator's own visits are exactly every 91-92 days."""
    rng = np.random.default_rng(1000 + seed)
    out = []
    for _, g in visits.groupby("patient", sort=True):
        g = g.assign(date=pd.to_datetime(g["date"])).sort_values("date").reset_index(drop=True)
        first, rest = g.iloc[:1], g.iloc[1:]
        rest = rest[rng.random(len(rest)) >= 0.10].copy()
        rest["date"] = rest["date"] + pd.to_timedelta(rng.integers(-20, 21, len(rest)), unit="D")
        g = pd.concat([first, rest]).sort_values("date").reset_index(drop=True)
        extra = [{"patient": g.loc[i, "patient"], "date": g.loc[i, "date"] + (g.loc[i + 1, "date"] - g.loc[i, "date"]) / 2,
                  "edss": g.loc[i, "edss"]} for i in range(len(g) - 1) if rng.random() < 0.40]
        g = pd.concat([g, pd.DataFrame(extra)]).sort_values("date")
        g["date"] = g["date"].dt.normalize()
        out.append(g.drop_duplicates(subset=["date"], keep="first"))
    return pd.concat(out, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500, help="patients per arm")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--perturb", action="store_true", help="irregular visits (see perturb())")
    a = ap.parse_args()
    tag = "_perturbed" if a.perturb else ""
    OUT.mkdir(exist_ok=True)
    rows = []
    for name, (rules, kw) in SETTINGS.items():
        dis = []
        tot = {"both": 0, "only_ms": 0, "only_py": 0, "patients": 0}
        for seed in a.seeds:
            t = simulate(Model(n_per_arm=a.n, seed=seed))
            v = perturb(t["visits"], seed) if a.perturb else t["visits"]
            py, ms = run_pymsprog(v, kw, t["relapses"]), run_msclinical(v, rules)
            m = ms.merge(py, on=["patient", "onset"], how="outer", indicator=True)
            tot["both"] += int((m["_merge"] == "both").sum())
            tot["only_ms"] += int((m["_merge"] == "left_only").sum())
            tot["only_py"] += int((m["_merge"] == "right_only").sum())
            tot["patients"] += v["patient"].nunique()
            d = m[m["_merge"] != "both"].assign(seed=seed)
            dis.append(d)
        n = tot["both"] + tot["only_ms"] + tot["only_py"]
        rows.append({"setting": name, **tot, "agreement": round(tot["both"] / n, 4) if n else None})
        pd.concat(dis).to_csv(OUT / f"disagreements_{name}{tag}.tsv", sep="\t", index=False)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"agreement{tag}.tsv", sep="\t", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
