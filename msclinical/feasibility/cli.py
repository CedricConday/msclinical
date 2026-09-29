"""msfeasibility FOLDER --criteria '{"age":[18,55],"edss":[0,5.5],"activity":{"relapses":1,"gd":1}}' [--index-date 2025-01-01]"""

import argparse
import json

from msclinical.tables.core import load

from .core import evaluate


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="msfeasibility", description=__doc__)
    ap.add_argument("folder", help="patients/visits/relapses/mri tables (mstables schema)")
    ap.add_argument("--criteria", required=True, help="JSON string or path to a JSON file")
    ap.add_argument("--index-date")
    a = ap.parse_args(argv)
    c = a.criteria
    if c.strip().startswith("{"):
        crit = json.loads(c)
    else:
        with open(c) as fh:
            crit = json.load(fh)
    _, funnel = evaluate(load(a.folder), crit, a.index_date)
    print(funnel.to_string(index=False))
    return 0
