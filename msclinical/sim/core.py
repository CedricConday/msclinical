"""A generative model with every cause recorded.

Per patient (arm effects are rate ratios):
* relapses: Poisson process, ``relapse_rate`` per year; each relapse leaves a residual EDSS step of
  0.5 or 1.0 with probability ``raw_residual_p`` (relapse-associated worsening, RAW);
* relapse-independent progression: Poisson process, ``pira_rate`` per year, each a 0.5 or 1.0 step
  (PIRA); a PIRA episode has MRI activity in its window with probability ``pira_mri_active_p``;
* MRI: new T2 lesions per scan interval, negative binomial around ``lesion_rate`` per year, plus the
  lesions each relapse brings (``relapse_lesions``);
* EDSS rating noise: each visit's rating is the true EDSS plus or minus 0.5 with probability
  ``rating_noise_p`` (inter-rater variability is of this order), clipped to 0..10.
Every true step is written to ``truth`` with its date, size and cause, so an attribution tool's output
can be scored against it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Model:
    n_per_arm: int = 200
    years: float = 2.0
    visit_months: int = 3
    mri_months: int = 6
    baseline_edss: tuple = (1.0, 5.5)
    relapse_rate: float = 0.35
    raw_residual_p: float = 0.35
    pira_rate: float = 0.12
    pira_mri_active_p: float = 0.3
    lesion_rate: float = 2.0
    lesion_k: float = 1.0
    relapse_lesions: float = 2.0
    rating_noise_p: float = 0.15
    arms: dict = field(default_factory=lambda: {"placebo": 1.0, "active": 0.5})  # rate ratio on relapses and lesions
    pira_effect: dict = field(default_factory=lambda: {"placebo": 1.0, "active": 0.8})
    start: str = "2024-01-01"
    seed: int = 0


def simulate(m: Model | None = None) -> dict:
    m = m or Model()
    rng = np.random.default_rng(m.seed)
    noise_rng = np.random.default_rng([m.seed, 1])  # own stream: the cohort is identical at every noise level
    t0 = pd.Timestamp(m.start)
    pts, vis, rel, mri, truth = [], [], [], [], []
    pid = 0
    days = int(m.years * 365.25)
    for arm, rr in m.arms.items():
        for _ in range(m.n_per_arm):
            pid += 1
            p = f"s{pid:04d}"
            e0 = rng.integers(int(m.baseline_edss[0] * 2), int(m.baseline_edss[1] * 2) + 1) / 2
            steps = []
            for d in np.sort(rng.uniform(0, days, rng.poisson(m.relapse_rate * rr * m.years))):
                rel.append({"patient": p, "date": (t0 + pd.Timedelta(days=int(d))).date()})
                if rng.random() < m.raw_residual_p:
                    size = rng.choice([0.5, 1.0])
                    if int(d) + 30 <= days:  # a residual step after follow-up ends is not observable truth
                        steps.append((int(d) + 30, size, "RAW", True))
            for d in np.sort(rng.uniform(0, days, rng.poisson(m.pira_rate * m.pira_effect[arm] * m.years))):
                steps.append((int(d), rng.choice([0.5, 1.0]), "PIRA", bool(rng.random() < m.pira_mri_active_p)))
            steps.sort()
            for d, s, cause, act in steps:
                truth.append({"patient": p, "arm": arm, "date": (t0 + pd.Timedelta(days=d)).date(), "step": float(s),
                              "cause": cause, "mri_active": act})
            def true_edss(day, e0=e0, steps=steps):
                return min(10.0, e0 + sum(s for d, s, _, _ in steps if d <= day))
            for k in range(0, int(m.years * 12) + 1, m.visit_months):
                day = int(k * 30.44)
                r = true_edss(day)
                u, sign = noise_rng.random(), noise_rng.choice([-0.5, 0.5])  # always drawn: the stream never shifts
                if u < m.rating_noise_p:
                    r = r + sign
                    if r == 0.5:  # 0.5 is not an EDSS step: a half-point misrating of 1.0 downward reads as 0
                        r = 0.0
                vis.append({"patient": p, "date": (t0 + pd.Timedelta(days=day)).date(), "edss": float(np.clip(r, 0, 10))})
            prev = 0
            rel_days = [int((pd.Timestamp(x["date"]) - t0).days) for x in rel if x["patient"] == p]
            for k in range(0, int(m.years * 12) + 1, m.mri_months):
                day = int(k * 30.44)
                if k == 0:
                    mri.append({"patient": p, "date": t0.date(), "new_t2": 0, "enlarging_t2": 0, "gd": 0})
                    continue
                dt = (day - prev) / 365.25
                mu = m.lesion_rate * rr * dt + m.relapse_lesions * sum(prev < d <= day for d in rel_days)
                n = int(rng.negative_binomial(m.lesion_k, m.lesion_k / (m.lesion_k + mu))) if mu > 0 else 0
                in_int = [(c, a) for d, _, c, a in steps if prev < d <= day]
                relapse_here = any(prev < d <= day for d in rel_days)
                if any(c == "PIRA" and a for c, a in in_int):
                    n = max(n, 1)
                elif any(c == "PIRA" and not a for c, a in in_int) and not relapse_here:
                    n = 0  # "no MRI activity" for an inactive PIRA, unless a relapse in the same interval brings lesions
                mri.append({"patient": p, "date": (t0 + pd.Timedelta(days=day)).date(), "new_t2": n, "enlarging_t2": 0,
                            "gd": int(rng.random() < min(0.9, 0.3 * n))})
                prev = day
            pts.append({"patient": p, "arm": arm, "start": t0.date(), "end": (t0 + pd.Timedelta(days=days)).date(), "baseline_edss": e0})
    return {k: pd.DataFrame(v) for k, v in (("patients", pts), ("visits", vis), ("relapses", rel), ("mri", mri), ("truth", truth))}


def write(tables: dict, out) -> None:
    from pathlib import Path

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for k, v in tables.items():
        v.to_csv(out / f"{k}.tsv", sep="\t", index=False)


def score_attribution(att: pd.DataFrame, truth: pd.DataFrame, window_days: int = 120) -> dict:
    """Match each attributed event (pirma output) to the nearest true step within ``window_days`` of onset."""
    tr = truth.assign(date=pd.to_datetime(truth["date"]))
    rows = []
    for _, e in att.iterrows():
        o = pd.Timestamp(e["onset_date"])
        c = tr[(tr["patient"] == e["patient"]) & ((tr["date"] - o).dt.days.abs() <= window_days)]
        if not len(c):
            rows.append({"called": e["class"], "true": "none (rating noise)"})
            continue
        best = c.iloc[(c["date"] - o).dt.days.abs().argmin()]
        true = best["cause"] if best["cause"] == "RAW" else ("PIRMA" if not best["mri_active"] else "PIRA_MRI_active")
        rows.append({"called": e["class"], "true": true})
    df = pd.DataFrame(rows)
    if not len(df):
        return {"events": 0}
    ct = pd.crosstab(df["called"], df["true"])
    return {"events": len(df), "agreement": float((df["called"] == df["true"]).mean()),
            "false_events_from_noise": int((df["true"] == "none (rating noise)").sum()), "confusion": ct.to_dict()}
