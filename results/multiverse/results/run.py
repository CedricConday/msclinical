"""Regenerates results/: 32 definitions on one trialtruth cohort (600 patients, 15 % misrating, seed 3)."""
from msmultiverse.core import influence, run
from trialtruth.core import Model, simulate

t = simulate(Model(n_per_arm=300, seed=3))
fs = t["patients"][["patient"]].assign(date="2023-06-01")
tab = run(t["visits"], t["relapses"], t["mri"], fs)
tab.to_csv("results/multiverse_trialtruth.tsv", sep="\t", index=False)
print(f"{len(tab)} definitions; events {tab['events'].min()}-{tab['events'].max()}; "
      f"PIRA share {tab['pira_share'].min():.2f}-{tab['pira_share'].max():.2f}; PIRMA {tab['pirma'].min()}-{tab['pirma'].max()}")
print("influence on PIRA share\n" + influence(tab).to_string(index=False))
print("influence on the PIRMA count\n" + influence(tab, "pirma").to_string(index=False))
