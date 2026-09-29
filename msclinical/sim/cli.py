"""trialtruth simulate --out dir [--n 200 --years 2 --seed 0]"""

import argparse

from .core import Model, simulate, write


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="trialtruth", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--out", required=True)
    s.add_argument("--n", type=int, default=200, help="patients per arm")
    s.add_argument("--years", type=float, default=2.0)
    s.add_argument("--noise", type=float, default=0.15, help="probability a visit's EDSS is misrated by 0.5")
    s.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    t = simulate(Model(n_per_arm=a.n, years=a.years, rating_noise_p=a.noise, seed=a.seed))
    write(t, a.out)
    print({k: len(v) for k, v in t.items()})
    return 0
