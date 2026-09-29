"""nedacalc visits.tsv relapses.tsv mri.tsv [--starts starts.tsv] [--rebaseline-months 6] [--years 2]"""

import argparse
import json
import sys

import pandas as pd

from .core import Params, cohort, neda


def _read(p):
    return pd.read_csv(p, sep="," if p.endswith(".csv") else "\t")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="nedacalc", description=__doc__)
    ap.add_argument("visits", help="patient, date, edss")
    ap.add_argument("relapses", help="patient, date")
    ap.add_argument("mri", help="patient, date, new_t2, enlarging_t2, gd [, pbvc (percent change since the previous scan)]")
    ap.add_argument("--starts", help="patient, date: start of observation (e.g. treatment start)")
    ap.add_argument("--window-months", type=int, default=12)
    ap.add_argument("--rebaseline-months", type=int, default=0)
    ap.add_argument("--years", type=int, default=1, help="cohort summary over the first N windows")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    t = neda(_read(a.visits), _read(a.relapses), _read(a.mri), _read(a.starts) if a.starts else None,
             Params(window_months=a.window_months, rebaseline_months=a.rebaseline_months))
    t.to_csv(a.out or sys.stdout, sep="\t", index=False)
    print(json.dumps(cohort(t, a.years)), file=sys.stderr)
    return 0
