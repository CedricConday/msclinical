"""msreportparse reports.tsv (patient, date, text) > mri.tsv   |   msreportparse --text "..." """

import argparse
import json
import sys

import pandas as pd

from .core import parse


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="msreportparse", description=__doc__)
    ap.add_argument("reports", nargs="?")
    ap.add_argument("--text")
    ap.add_argument("--no-evidence", action="store_true")
    a = ap.parse_args(argv)
    if a.text:
        print(json.dumps(parse(a.text), indent=2, ensure_ascii=False))
        return 0
    df = pd.read_csv(a.reports, sep="\t")
    rows = [{"patient": r["patient"], "date": r["date"], **parse(str(r["text"]))} for _, r in df.iterrows()]
    out = pd.DataFrame(rows)
    if a.no_evidence:
        out = out[[c for c in out if not c.endswith(("_evidence", "_approx"))]]
    out.to_csv(sys.stdout, sep="\t", index=False)
    return 0
