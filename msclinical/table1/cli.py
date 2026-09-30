"""mstable1 FOLDER [--group arm] [--markdown]"""

import argparse

from msclinical.tables.core import load

from .core import table1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mstable1", description=__doc__)
    ap.add_argument("folder")
    ap.add_argument("--group", default="arm")
    ap.add_argument("--markdown", action="store_true")
    a = ap.parse_args(argv)
    t = table1(load(a.folder), a.group)
    print(t.to_markdown() if a.markdown else t.to_string())
    return 0
