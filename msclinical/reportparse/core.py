"""Rule-based extraction of MS MRI activity from report text.

Per report: new_t2, enlarging_t2, gd, each an integer, 0 when the report states their absence, or None when
the report does not mention them. Every value keeps the sentence it came from, so a reviewer checks
evidence instead of trusting a number. Number words (one..ten, ein..zehn) and "multiple/several/mehrere"
(recorded as ">=2", count 2, flagged) are understood. Rules are regular expressions in RULES; a sentence
matched by a negation pattern sets the field to 0 unless a count in the same sentence says otherwise.
"""

from __future__ import annotations

import re

NUM = {"one": 1, "a": 1, "an": 1, "single": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7,
       "eight": 8, "nine": 9, "ten": 10, "ein": 1, "eine": 1, "einer": 1, "einen": 1, "zwei": 2, "drei": 3, "vier": 4,
       "fünf": 5, "sechs": 6, "sieben": 7, "acht": 8, "neun": 9, "zehn": 10}
SEVERAL = r"(?:multiple|several|numerous|mehrere|multiple|zahlreiche)"
COUNT = r"\b(\d{1,2}(?!\d)|" + "|".join(sorted(NUM, key=len, reverse=True)) + "|" + SEVERAL + r")"

LESION = r"(?:lesions?|foci|focus|herde?n?|läsionen|läsion|hyperintensit(?:y|ies|äten))"
RULES = {
    "gd": {
        "pos": [COUNT + r"\s+(?:new\s+)?(?:gadolinium[- ]|gd[- ]?|contrast[- ])?enhancing\s+" + LESION,
                COUNT + r"\s+(?:kontrastmittel|km|gadolinium)[- ]?(?:aufnehmende[nr]?|anreichernde[nr]?)\s+" + LESION],
        "neg": [r"\bno\s+(?:new\s+)?(?:gadolinium[- ]|gd[- ]?|contrast[- ])?enhanc", r"without\s+(?:any\s+)?enhancement",
                r"\bkeine?\s+(?:kontrastmittel|km|gadolinium)[- ]?(?:aufnehmende|anreichernde|aufnahme|anreicherung)",
                r"ohne\s+(?:kontrastmittel|km)[- ]?(?:aufnahme|anreicherung)", r"no\s+abnormal\s+enhancement"],
    },
    "new_t2": {
        "pos": [COUNT + r"\s+new\s+(?:t2[- ]?(?:/flair)?[- ]?(?:hyperintense\s+)?|flair\s+|white\s+matter\s+|demyelinating\s+)?" + LESION,
                COUNT + r"\s+neue[nr]?\s+(?:t2[- ]?(?:/flair)?[- ]?(?:hyperintense[nr]?\s+)?|flair[- ]?)?" + LESION],
        "neg": [r"\bno\s+new\s+(?:or\s+enlarging\s+)?(?:t2|flair|white|demyelinating|lesion|foci)", r"\bkeine?\s+neue[nr]?\s",
                r"no\s+(?:interval\s+)?(?:change|progression)", r"unchanged\s+lesion\s+(?:load|burden)", r"befund\s+(?:konstant|stabil|unverändert)",
                r"stable\s+(?:lesion|white\s+matter|disease)", r"unveränderte[rn]?\s+(?:läsion|herd)"],
    },
    "enlarging_t2": {
        "pos": [COUNT + r"\s+(?:enlarging|enlarged|growing)\s+" + LESION, COUNT + r"\s+(?:größenprogrediente[nr]?|vergrößerte[nr]?)\s+" + LESION],
        "neg": [r"no\s+(?:new\s+or\s+)?(?:enlarging|enlarged|growing)", r"keine?\s+(?:größenprogrediente|vergrößerte)", r"no\s+(?:interval\s+)?(?:change|progression)",
                r"befund\s+(?:konstant|stabil|unverändert)", r"stable\s+(?:lesion|white\s+matter|disease)"],
    },
}


def _num(tok: str) -> tuple[int, bool]:
    t = tok.lower()
    if t.isdigit():
        return int(t), False
    if t in NUM:
        return NUM[t], False
    return 2, True  # "several": at least two


def sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.;!?\n])\s+", text) if s.strip()]


PRIOR = re.compile(r"(?:seen|noted|present|visible|described|were|was)\s+(?:on|in)\s+(?:the\s+)?(?:prior|previous|earlier)\s+"
                   r"(?:study|scan|mri|examination)|in\s+der\s+voruntersuchung\s+(?:gesehen|beschrieben|abgrenzbar)|vorbekannt")
STABLE_ONLY = re.compile(r"stable|konstant|unverändert|no\s+(?:interval\s+)?(?:change|progression)")


def parse(text: str) -> dict:
    """Counts per sentence; the report value is the largest sentence count (Findings and Impression usually
    repeat the same lesions, so summing would double them). Sentences about a prior study are skipped."""
    out = {}
    sents = [x for x in sentences(text) if not PRIOR.search(x.lower())]
    positives = {f: any(re.search(p, x.lower()) for x in sents for p in r["pos"]) for f, r in RULES.items()}
    for field, r in RULES.items():
        val, ev, approx = None, [], False
        for s in sents:
            low = s.lower()
            hits = [m for p in r["pos"] for m in re.finditer(p, low)]
            if hits:
                n = 0
                for m in hits:
                    k, a = _num(m.group(1))
                    n += k
                    approx |= a
                val = max(val or 0, n)
                ev.append(s)
                continue
            neg = [p for p in r["neg"] if re.search(p, low)]
            if neg and val is None:
                # a generic "stable" statement does not zero a field when the same report counts another
                # kind of activity: "Stable disease. One new enhancing lesion." is not zero new lesions
                generic = all(STABLE_ONLY.search(re.search(p, low).group(0)) for p in neg)
                if generic and any(positives[f] for f in positives if f != field):
                    continue
                val = 0
                ev.append(s)
        out[field] = val
        out[f"{field}_evidence"] = " | ".join(ev)
        out[f"{field}_approx"] = approx
    return out
