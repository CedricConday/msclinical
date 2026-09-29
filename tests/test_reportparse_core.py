from msclinical.reportparse.core import parse

CASES = [
    (("Compared with the prior study there are two new T2 hyperintense lesions in the periventricular white matter. "
      "No gadolinium-enhancing lesions."), {"new_t2": 2, "gd": 0, "enlarging_t2": None}),
    ("Stable lesion load. No new or enlarging lesions. No abnormal enhancement.", {"new_t2": 0, "enlarging_t2": 0, "gd": 0}),
    ("One enhancing lesion in the left centrum semiovale; three new FLAIR lesions and one enlarging lesion.",
     {"gd": 1, "new_t2": 3, "enlarging_t2": 1}),
    ("Im Vergleich zur Voruntersuchung zwei neue T2-hyperintense Läsionen. Keine KM-aufnehmenden Herde.", {"new_t2": 2, "gd": 0}),
    ("Befund konstant. Keine neuen Läsionen, keine Kontrastmittel-Anreicherung.", {"new_t2": 0, "gd": 0}),
    ("Mehrere neue Läsionen, davon eine kontrastmittelaufnehmende Läsion.", {"new_t2": 2, "gd": 1}),
    ("Brain MRI with contrast. Normal ventricles.", {"new_t2": None, "gd": None, "enlarging_t2": None}),
]


def test_cases():
    for text, want in CASES:
        got = parse(text)
        for k, v in want.items():
            assert got[k] == v, (text, k, got[k], v)
    assert parse("Several new lesions.")["new_t2_approx"] is True



def test_review_regressions():
    assert parse("Compared with 2023 new lesions are seen.")["new_t2"] is None
    assert parse("Im Vergleich zu 2021 neue Läsionen.")["new_t2"] is None
    assert parse("Findings: two new lesions. Impression: two new lesions.")["new_t2"] == 2
    r = parse("No new or enlarging T2 lesions.")
    assert r["new_t2"] == 0 and r["enlarging_t2"] == 0
    assert parse("No new enhancing lesions.")["gd"] == 0
    r = parse("Stable disease. One new enhancing lesion.")
    assert r["gd"] == 1 and r["new_t2"] is None
    assert parse("Several enhancing lesions were seen on the prior study; none now.")["gd"] is None
