"""msmultiverse visits.tsv relapses.tsv [--mri mri.tsv] [--followup-start s.tsv] [--out table.tsv]"""

import argparse

import pandas as pd

from .core import influence, run


def _read(p):
    return pd.read_csv(p, sep="," if p.endswith(".csv") else "\t")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="msmultiverse", description=__doc__)
    ap.add_argument("visits")
    ap.add_argument("relapses")
    ap.add_argument("--mri")
    ap.add_argument("--followup-start")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    t = run(_read(a.visits), _read(a.relapses), _read(a.mri) if a.mri else None, _read(a.followup_start) if a.followup_start else None)
    if a.out:
        t.to_csv(a.out, sep="\t", index=False)
    print(f"{len(t)} definitions; events {t['events'].min()}-{t['events'].max()}; "
          f"PIRA share {t['pira_share'].min():.2f}-{t['pira_share'].max():.2f}")
    print(influence(t).to_string(index=False))
    return 0
