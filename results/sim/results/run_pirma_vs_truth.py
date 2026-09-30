"""Regenerates results/pirma_vs_truth.txt: pirma attribution against trialtruth's recorded causes."""
import pandas as pd
from edssplus.core import Rules, events
from pirma.core import attribute
from trialtruth.core import Model, score_attribution, simulate

for noise in (0.0, 0.15, 0.3):
    for cd in (84, 168):
        t = simulate(Model(n_per_arm=300, rating_noise_p=noise, seed=3))
        ev = events(t["visits"], Rules(confirm_days=cd), components=("edss",))
        fs = t["patients"][["patient"]].assign(date=pd.Timestamp("2023-06-01"))
        att = attribute(ev, t["relapses"], t["mri"], fs)
        s = score_attribution(att, t["truth"])
        print(f"noise {noise:.2f} confirm {cd}d: true steps {len(t['truth'])}, events {s['events']}, "
              f"from noise {s['false_events_from_noise']}, agreement {s['agreement']:.2f}, called {att['class'].value_counts().to_dict()}")
