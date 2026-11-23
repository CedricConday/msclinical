"""Track A, part 2: RAW / PIRA attribution, msclinical.pirma vs pymsprog (RAW_PIRA=True), on the same CDW events.

CDW detection is matched as in compare_cdw.py (roving, all events, msclinical's EDSS thresholds, no relapse
exclusions), so any difference here is the attribution rule alone. msclinical: RAW if a relapse falls in
[onset - 90 d, confirmation]; PIRA if none does (MRI not used here: class PIRA_MRI_unknown counts as PIRA).
pymsprog: two settings, its defaults and the closest match its API allows.

Output: attribution.tsv and attr_disagreements_<setting>.tsv next to this script; scored against simulator truth too.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
import warnings
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from compare_cdw import INF, ms_delta
from pymsprog import MSprog

from msclinical.edssplus.core import Rules, events
from msclinical.pirma.core import attribute
from msclinical.sim.core import Model, simulate

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent
CDW = {"event": "multiple", "baseline": "roving", "proceed_from": "firstconf", "sub_threshold_rebl": "improvement",
           "conf_days": 84, "conf_tol_days": [0, INF], "delta_fun": ms_delta, "relapse_to_bl": 0, "relapse_to_event": 0,
           "relapse_to_conf": 0, "RAW_PIRA": True}
PY = {
    # pymsprog's real defaults for everything (CDW detection included: own EDSS table, confirmation -7/+730 d,
    # 30-d relapse exclusions for baseline and confirmation visits; RAW: relapse up to 90 d before onset; PIRA:
    # relapse-free over [onset - 90, onset + 30] and [conf - 90, conf + 30]); only the event scan is set to match.
    "py_defaults": {"_defaults": True},
    # RAW: relapse up to 90 d before onset or up to 84 d after it; PIRA: no relapse from onset - 90 d to confirmation
    "py_matched": {"relapse_assoc": [90, 84], "relapse_indep": {"prec": (0, 0), "event": (90, None),
                                                              "conf": (None, 0), "prec_type": "baseline"}},
}


def py_attr(v, rel, kw):
    with contextlib.redirect_stdout(io.StringIO()):
        _, r = MSprog(v, "patient", "edss", "date", "edss", relapse=rel, rsubj_col="patient", rdate_col="date",
                      include_dates=True, include_stable=False, verbose=0,
                      **({"event": "multiple", "baseline": "roving", "RAW_PIRA": True} if kw.get("_defaults") else {**CDW, **kw}))
    r = r[r["event_type"] == "CDW"]  # CDW_type is PIRA, RAW or undefined
    return pd.DataFrame({"patient": r["patient"], "onset": pd.to_datetime(r["date"]).dt.date,
                         "py_class": r["CDW_type"]})


def ms_attr(v, rel):
    ev = events(v, Rules(reference="roving"), components=("edss",))
    a = attribute(ev, rel)
    if not len(a):
        return pd.DataFrame(columns=["patient", "onset", "ms_class"])
    cls = a["class"].replace({"PIRA_MRI_unknown": "PIRA"})
    return pd.DataFrame({"patient": a["patient"], "onset": pd.to_datetime(a["onset_date"]).dt.date, "ms_class": cls})


def truth_cause(truth, ms):  # cause of the true step nearest before each onset (within 180 d), else "noise"
    t = truth.assign(date=pd.to_datetime(truth["date"]))
    out = []
    for _, e in ms.iterrows():
        tt = t[(t["patient"] == e["patient"]) & (t["date"] <= pd.Timestamp(e["onset"]))
               & (t["date"] > pd.Timestamp(e["onset"]) - pd.Timedelta(days=180))]
        out.append(tt.sort_values("date")["cause"].iloc[-1] if len(tt) else "noise")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    a = ap.parse_args()
    rows = []
    for name, kw in PY.items():
        allm = []
        for seed in a.seeds:
            t = simulate(Model(n_per_arm=a.n, seed=seed))
            v, rel = t["visits"], t["relapses"]
            m = ms_attr(v, rel).merge(py_attr(v, rel, kw), on=["patient", "onset"], how="outer")
            m["truth"] = truth_cause(t["truth"], m)
            allm.append(m.assign(seed=seed))
        m = pd.concat(allm)
        both = m.dropna(subset=["ms_class", "py_class"])
        tab = pd.crosstab(both["ms_class"], both["py_class"])
        print(f"\n== {name}: {len(both)} shared events, {m['py_class'].isna().sum()} msclinical-only, "
              f"{m['ms_class'].isna().sum()} pymsprog-only\n{tab}")
        for tool in ("ms_class", "py_class"):
            s = both[both["truth"].isin(["PIRA", "RAW"])]
            acc = (s[tool] == s["truth"]).mean()
            print(f"   {tool} matches simulator truth on {acc:.1%} of {len(s)} events with a true cause")
        rows.append({"setting": name, "shared": len(both), "same_class": int((both["ms_class"] == both["py_class"]).sum()),
                     "agreement": round((both["ms_class"] == both["py_class"]).mean(), 4)})
        m[m["ms_class"] != m["py_class"]].to_csv(OUT / f"attr_disagreements_{name}.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(OUT / "attribution.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
