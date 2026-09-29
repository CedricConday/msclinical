"""dmtepochs treatments.tsv followup.tsv [--relapses relapses.tsv] [--rebaseline-months 6 --carryover-days 90]"""

import argparse
import sys

import pandas as pd

from .core import Rules, arr_by_drug, epochs


def _read(p):
    return pd.read_csv(p, sep="," if p.endswith(".csv") else "\t")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dmtepochs", description=__doc__)
    ap.add_argument("treatments")
    ap.add_argument("followup", help="patient, start, end")
    ap.add_argument("--relapses")
    ap.add_argument("--rebaseline-months", type=int, default=6)
    ap.add_argument("--carryover-days", type=int, default=90)
    ap.add_argument("--gap-days", type=int, default=30)
    a = ap.parse_args(argv)
    ep = epochs(_read(a.treatments), _read(a.followup), Rules(a.gap_days, a.rebaseline_months, a.carryover_days))
    ep.to_csv(sys.stdout, sep="\t", index=False)
    if a.relapses:
        print(arr_by_drug(_read(a.relapses), ep).to_string(), file=sys.stderr)
    return 0
