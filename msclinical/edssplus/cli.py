"""edssplus events visits.tsv [--confirm-days 84|168] [--reference roving|fixed] [--composite]"""

import argparse
import sys

import pandas as pd

from .core import Rules, edss_plus, events


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="edssplus", description=__doc__)
    ap.add_argument("visits", help="TSV/CSV with patient, date and any of edss, t25fw, hpt9_dom, hpt9_nondom")
    ap.add_argument("--confirm-days", type=int, default=84)
    ap.add_argument("--reference", choices=("roving", "fixed"), default="roving")
    ap.add_argument("--no-rove-on-improvement", action="store_true")
    ap.add_argument("--composite", action="store_true", help="one row per patient: first confirmed event on any component")
    ap.add_argument("--out", help="write TSV here (default stdout)")
    a = ap.parse_args(argv)
    sep = "," if a.visits.endswith(".csv") else "\t"
    v = pd.read_csv(a.visits, sep=sep)
    r = Rules(confirm_days=a.confirm_days, reference=a.reference, rove_on_improvement=not a.no_rove_on_improvement)
    df = edss_plus(v, r) if a.composite else events(v, r)
    df.to_csv(a.out or sys.stdout, sep="\t", index=False)
    return 0
