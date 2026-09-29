"""pirma visits.tsv relapses.tsv [--mri mri.tsv] [--confirm-days 84] [--reference roving|fixed]"""

import argparse
import json
import sys

import pandas as pd

from msclinical.edssplus.core import Rules, events

from .core import Windows, attribute, summary


def _read(p):
    return pd.read_csv(p, sep="," if p.endswith(".csv") else "\t")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="pirma", description=__doc__)
    ap.add_argument("visits")
    ap.add_argument("relapses", help="patient, date (relapse onset)")
    ap.add_argument("--mri", help="patient, date, new_t2, enlarging_t2, gd")
    ap.add_argument("--followup-start", help="patient, date: start of relapse ascertainment")
    ap.add_argument("--confirm-days", type=int, default=84)
    ap.add_argument("--reference", choices=("roving", "fixed"), default="roving")
    ap.add_argument("--component", default="edss", help="edss (default), t25fw, hpt9_dom, hpt9_nondom or all")
    ap.add_argument("--raw-days", type=int, default=90)
    ap.add_argument("--mri-window", choices=("reference", "onset"), default="reference")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    vis = _read(a.visits)
    ev = events(vis, Rules(confirm_days=a.confirm_days, reference=a.reference))
    first = pd.to_datetime(vis["date"]).groupby(vis["patient"]).min()
    if a.component != "all":
        ev = ev[ev["component"] == a.component]
    att = attribute(ev, _read(a.relapses), _read(a.mri) if a.mri else None,
                    _read(a.followup_start) if a.followup_start else None,
                    Windows(raw_days=a.raw_days, mri_window=a.mri_window), visits_first=first)
    att.to_csv(a.out or sys.stdout, sep="\t", index=False)
    print(json.dumps(summary(att)), file=sys.stderr)
    return 0
