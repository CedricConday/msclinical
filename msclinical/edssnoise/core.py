"""EDSS rating noise from visit sequences alone.

True disability rarely goes down and then straight back up within a few months; a half-point
misrating does that all the time. The reversal rate of a dataset (visits whose EDSS differs from both
neighbours in the same direction, a local peak or trough) therefore carries the misrating rate. The
mapping from reversal rate to misrating probability depends on visit spacing and true progression, so
it is calibrated by simulation (trialtruth) with the dataset's own number of patients, visits and
spacing, and the estimate is the misrating probability whose simulated reversal rate matches the
observed one, with a bootstrap interval over patients.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def reversal_rate(visits: pd.DataFrame) -> tuple[float, int]:
    n_rev, n_mid = 0, 0
    for _, g in visits.sort_values(["patient", "date"]).groupby("patient"):
        e = g["edss"].to_numpy(dtype=float)
        if len(e) < 3:
            continue
        a, b, c = e[:-2], e[1:-1], e[2:]
        rev = ((b > a) & (b > c)) | ((b < a) & (b < c))
        n_rev += int(rev.sum())
        n_mid += len(b)
    return (n_rev / n_mid if n_mid else float("nan")), n_mid


def _design(visits: pd.DataFrame) -> dict:
    v = visits.assign(date=pd.to_datetime(visits["date"])).sort_values(["patient", "date"])
    gaps = v.groupby("patient")["date"].diff().dt.days.dropna()
    span = v.groupby("patient")["date"].agg(lambda s: (s.max() - s.min()).days / 365.25)
    return {"patients": v["patient"].nunique(), "visit_months": max(1, round(gaps.median() / 30.44)),
            "years": max(0.5, float(span.median()))}


GRID = tuple(float(x) for x in np.linspace(0, 0.4, 17))


def calibrate(design: dict, grid=GRID, reps: int = 5, seed: int = 0) -> pd.DataFrame:
    from msclinical.sim.core import Model, simulate

    rows = []
    n = max(50, design["patients"] // 2)
    if design["years"] * 12 / design["visit_months"] < 2:
        raise ValueError(f"design has fewer than 3 visits per patient ({design}); reversals cannot be counted")
    for p in grid:
        rs = []
        for r in range(reps):
            t = simulate(Model(n_per_arm=n, years=design["years"], visit_months=design["visit_months"],
                               rating_noise_p=float(p), seed=seed + r))
            rs.append(reversal_rate(t["visits"])[0])
        rows.append({"misrating_p": float(p), "reversal_rate": float(np.mean(rs))})
    cal = pd.DataFrame(rows)
    # the same seeds at every grid point (trialtruth keeps the cohort fixed across noise levels), then forced
    # monotone so that inverting it is well defined
    cal["reversal_rate"] = np.maximum.accumulate(cal["reversal_rate"].to_numpy())
    return cal


def estimate(visits: pd.DataFrame, boot: int = 200, seed: int = 0) -> dict:
    obs, n_mid = reversal_rate(visits)
    if not n_mid:
        raise ValueError("no patient has three or more EDSS visits; the misrating rate cannot be estimated")
    d = _design(visits)
    cal = calibrate(d)
    inv = lambda r: float(np.interp(r, cal["reversal_rate"], cal["misrating_p"]))
    rng = np.random.default_rng(seed)
    pats = visits["patient"].unique()
    bs = []
    for _ in range(boot):
        pick = rng.choice(pats, len(pats), replace=True)
        bv = pd.concat([visits[visits["patient"] == p].assign(patient=f"{p}_{i}") for i, p in enumerate(pick)])
        bs.append(inv(reversal_rate(bv)[0]))
    return {"reversal_rate": obs, "interior_visits": n_mid, "design": d, "misrating_p": inv(obs),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "calibration": cal.to_dict("records"),
            "note": "misrating = a visit rated half a point off the true EDSS; the model and its progression rates are trialtruth defaults"}
