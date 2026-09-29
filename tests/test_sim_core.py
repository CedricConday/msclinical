from msclinical.sim.core import Model, simulate


def test_tables_and_truth():
    t = simulate(Model(n_per_arm=50, seed=1))
    assert set(t) == {"patients", "visits", "relapses", "mri", "truth"}
    assert len(t["patients"]) == 100 and t["visits"]["edss"].between(0, 10).all()
    assert set(t["truth"]["cause"]) <= {"RAW", "PIRA"}
    arr = t["relapses"].merge(t["patients"], on="patient").groupby("arm").size()
    assert arr["active"] < arr["placebo"]


def test_noise_free_edss_moves_only_with_truth():
    t = simulate(Model(n_per_arm=20, rating_noise_p=0.0, seed=2))
    v = t["visits"].sort_values(["patient", "date"])
    moved = v.groupby("patient")["edss"].apply(lambda s: (s.diff().fillna(0) != 0).any())
    has_step = t["patients"]["patient"].isin(t["truth"]["patient"]).values
    assert not (moved.values & ~has_step).any()


def test_edss_on_scale():
    t = simulate(Model(n_per_arm=100, rating_noise_p=0.5, seed=4))
    assert not (t["visits"]["edss"] == 0.5).any()


def test_noise_does_not_change_the_cohort():
    a = simulate(Model(n_per_arm=100, rating_noise_p=0.0, seed=7))
    b = simulate(Model(n_per_arm=100, rating_noise_p=0.3, seed=7))
    assert a["truth"].equals(b["truth"]) and a["relapses"].equals(b["relapses"]) and a["mri"].equals(b["mri"])
