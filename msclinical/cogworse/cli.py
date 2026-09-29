"""cogworse visits.tsv (patient, date, sdmt) [--points 4 --pct 10 --combine or|and --confirm-days 84]"""

import argparse
import sys

import pandas as pd

from .core import Rules, events


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="cogworse", description=__doc__)
    ap.add_argument("visits")
    ap.add_argument("--points", type=float, default=4.0)
    ap.add_argument("--pct", type=float, default=10.0)
    ap.add_argument("--combine", choices=("or", "and"), default="or")
    ap.add_argument("--confirm-days", type=int, default=84)
    ap.add_argument("--baseline-visits", type=int, default=2)
    ap.add_argument("--fixed", action="store_true")
    a = ap.parse_args(argv)
    v = pd.read_csv(a.visits, sep="," if a.visits.endswith(".csv") else "\t")
    events(v, Rules(a.points, a.pct, a.combine, a.confirm_days, a.baseline_visits, not a.fixed)).to_csv(sys.stdout, sep="\t", index=False)
    return 0
