"""mstables check FOLDER   (patients.tsv, visits.tsv, relapses.tsv, mri.tsv, treatments.tsv; .csv also read)"""

import argparse
import sys

from .core import load, validate


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mstables", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("folder")
    a = ap.parse_args(argv)
    t = load(a.folder)
    if not t:
        print("no known tables found", file=sys.stderr)
        return 2
    issues = validate(t)
    for i in issues:
        print(f"{i['level']:7s} {i['table']}: {i['message']}" + (f" (rows {i['rows'][:5]})" if i["rows"] else ""))
    errors = sum(i["level"] == "error" for i in issues)
    print(f"{len(t)} tables, {errors} errors, {len(issues) - errors} other", file=sys.stderr)
    return 1 if errors else 0
