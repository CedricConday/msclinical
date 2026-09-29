"""edssnoise visits.tsv (patient, date, edss) [--boot 200]"""

import argparse
import json

import pandas as pd

from .core import estimate


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="edssnoise", description=__doc__)
    ap.add_argument("visits")
    ap.add_argument("--boot", type=int, default=200)
    a = ap.parse_args(argv)
    v = pd.read_csv(a.visits, sep="," if a.visits.endswith(".csv") else "\t")
    r = estimate(v, a.boot)
    r.pop("calibration")
    print(json.dumps(r, indent=2))
    return 0
