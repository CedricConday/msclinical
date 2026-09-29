"""relapsecheck relapses.tsv [--followup f.tsv] [--out clean.tsv] [--removed removed.tsv]"""

import argparse
import json
import sys

import pandas as pd

from .core import check


def _read(p):
    return pd.read_csv(p, sep="," if p.endswith(".csv") else "\t")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="relapsecheck", description=__doc__)
    ap.add_argument("relapses", help="patient, date [, fever, infection, duration_hours]")
    ap.add_argument("--followup", help="patient, start, end")
    ap.add_argument("--min-gap-days", type=int, default=30)
    ap.add_argument("--keep-pseudo", action="store_true")
    ap.add_argument("--dayfirst", action="store_true", help="read 05.01.2020 as 5 January")
    ap.add_argument("--date-format", help='e.g. "%%d.%%m.%%Y"')
    ap.add_argument("--out")
    ap.add_argument("--removed")
    a = ap.parse_args(argv)
    clean, removed, s = check(_read(a.relapses), _read(a.followup) if a.followup else None, a.min_gap_days, not a.keep_pseudo,
                              dayfirst=a.dayfirst, date_format=a.date_format)
    clean.to_csv(a.out or sys.stdout, sep="\t", index=False)
    if a.removed:
        removed.to_csv(a.removed, sep="\t", index=False)
    print(json.dumps(s), file=sys.stderr)
    return 0
