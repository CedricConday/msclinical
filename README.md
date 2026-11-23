# msclinical

The clinical-rule layer of MS research as one package. Every rule is in one place, every tool reads and writes the same tables (the `tables` schema), and every rule has been checked against simulated data with known causes (`sim`).

```bash
pip install msclinical
msclinical <tool> --help
```

**No install:** [worsening, PIRA, NEDA, relapse cleaning and SDMT worsening in your browser](https://cedricconday.github.io/msclinical/browser/)
(Pyodide; nothing is uploaded). Overview of all six MS packages: https://cedricconday.github.io/msclinical/

| Tool | Command | What it does |
|---|---|---|
| [`edssplus`](#edssplus) | `msclinical edssplus` | confirmed worsening on EDSS, T25FW, 9HPT; EDSS-plus; roving or fixed reference |
| [`pirma`](#pirma) | `msclinical pirma` | RAW, PIRA and PIRMA attribution of every confirmed worsening event |
| [`cogworse`](#cogworse) | `msclinical cogworse` | confirmed SDMT worsening with a practice-effect-aware reference |
| [`neda`](#neda) | `msclinical neda` | NEDA-3 and NEDA-4 per patient and year, unknowns kept unknown |
| [`relapses`](#relapses) | `msclinical relapses` | relapse table hygiene: 30-day merge, pseudo-relapses, follow-up bounds |
| [`epochs`](#epochs) | `msclinical epochs` | treatment epochs with re-baselining and carry-over; ARR per drug |
| [`edssnoise`](#edssnoise) | `msclinical edssnoise` | EDSS misrating rate of a registry, from visit reversals |
| [`multiverse`](#multiverse) | `msclinical multiverse` | every common PIRA definition on one dataset, and which choice moves the answer |
| [`tables`](#tables) | `msclinical tables` | the shared table schema and a validator |
| [`feasibility`](#feasibility) | `msclinical feasibility` | registry funnel for a trial's inclusion criteria |
| [`table1`](#table1) | `msclinical table1` | Table 1 of a cohort paper, per arm, with missing counts |
| [`reportparse`](#reportparse) | `msclinical reportparse` | MRI activity from English and German report text, with evidence sentences |
| [`sim`](#sim) | `msclinical sim` | simulated cohorts with every cause recorded (formerly sim) |

Research software, not a medical device. MIT licence. Each section below keeps the tool's own results and the
corrections that three independent reviews made to them.

## edssplus


Confirmed disability worsening (CDW) in multiple sclerosis from a plain visit table: EDSS, Timed
25-Foot Walk and 9-Hole Peg Test, one component at a time or as the **EDSS-plus** composite, with a
fixed or **roving** reference and 12- or 24-week confirmation. Every rule is in one dataclass.

```bash
pip install msclinical
msclinical edssplus visits.tsv --confirm-days 84 --reference roving           # every event, every component
msclinical edssplus visits.tsv --composite                                     # first event per patient (time-to-event input)
```

`visits.tsv`: `patient, date` plus any of `edss, t25fw, hpt9_dom, hpt9_nondom` (seconds for the timed tests).

Rules (`msclinical/edssplus/core.py`, `Rules`): EDSS +1.5 from 0, +1.0 from 1.0 to 5.5, +0.5 above 5.5; timed
tests at least 20 % slower (Cadavid et al. 2017); the worsening must hold at every visit up to the first
visit at least `confirm_days` after onset. A worsening still open at the end of follow-up is
`unconfirmed`, never an event. Roving reference as in Kappos et al. 2018, including re-baselining on
confirmed improvement: a lower value becomes the reference when every visit up to the first one at least
`confirm_days` later stays at or below it, and the search resumes from that lower visit.

Two choices here differ from other open tools and change results. The +1.0 step runs to a reference of 5.5;
[pymsprog](https://github.com/noemimontobbio/pymsprog) stops it at 5.0 (Lorscheider et al. 2016), so patients
at 5.5 get an earlier event there. Set `Rules(edss_steps=((0.0, 1.5), (5.0, 1.0), (10.0, 0.5)))` to match it.
After a confirmed improvement, pymsprog moves the reference to the confirmation visit instead; on simulated
cohorts this changes about one event in 600.

Relapse attribution (RAW, PIRA, PIRMA) is the job of the companion package
[pirma](https://github.com/CedricConday/msclinical#pirma); NEDA-3/4 is [neda](https://github.com/CedricConday/msclinical#neda).

## pirma


Every confirmed disability worsening event in a multiple sclerosis cohort, attributed:
**RAW** (relapse-associated), **PIRA** (progression independent of relapse activity) and
**PIRMA** (independent of relapse *and* MRI activity). Built on
[edssplus](https://github.com/CedricConday/msclinical#edssplus) for the events.

```bash
pip install msclinical
msclinical pirma visits.tsv relapses.tsv --mri mri.tsv --followup-start start.tsv --confirm-days 84
```

The MRI table has `patient, date, new_t2, enlarging_t2, gd`. The rules are in `msclinical/pirma/core.py`:
a relapse within 90 days before onset or between onset and confirmation makes the event RAW; otherwise
it is PIRA when relapse ascertainment covers that window, else `undetermined`. PIRA becomes PIRMA
only when at least one scan lies between the reference visit and confirmation and none shows new,
enlarging or enhancing lesions. PIRA with no scan in that window stays `PIRA_MRI_unknown`: the
package never infers "no MRI activity" from missing MRI.

Existing open code (Bsteh, Marti, Hoepner; drstrupf/multiple-sclerosis-disability-progression)
separates PIRA from RAW on EDSS; pirma adds the MRI-activity split and the timed-test components.

### Checked against pymsprog (`results/pymsprog/`)

[pymsprog](https://github.com/noemimontobbio/pymsprog) 1.0.7 (Montobbio et al., Mult Scler 2024) is an
independent implementation of the same rules. On the same 3,000 simulated patients from
[sim](https://github.com/CedricConday/msclinical#sim), with its settings matched to these:

| | Same events | Differ |
|---|---:|---:|
| First CDW, fixed reference | 539 | 0 |
| First CDW, fixed reference, irregular visits (dates ±20 days, 10 % missed, interim visits) | 492 | 0 |
| All CDW, roving reference | 601 | 1 (re-baselining after improvement, see edssplus) |
| RAW/PIRA class of those 601 events | 596 | 5 |

The five are all a relapse 85 to 92 days after onset, just before confirmation: pirma calls the event RAW,
pymsprog's window counted from onset leaves it unclassified. With every pymsprog setting at its default,
94.7 % of shared events get the same class; the rest come from its different windows, listed in
`results/pymsprog/compare_attribution.py`. Both tools give the simulator's true cause for about 81 % of
events. Rerun with `pip install pymsprog==1.0.7` and the two scripts in that folder.

## cogworse


Disability worsening in MS is usually counted on EDSS, which barely sees cognition. cogworse detects
**confirmed cognitive worsening** from Symbol Digit Modalities Test (SDMT) visits, the same way
[edssplus](https://github.com/CedricConday/msclinical#edssplus) does for EDSS. A decline of at least 4 points or 10 %
(Benedict et al. 2017; either or both, configurable) counts only if it holds at every visit up to a
confirmation at least 12 weeks later.

```bash
pip install msclinical
msclinical cogworse visits.tsv                          # patient, date, sdmt
msclinical cogworse visits.tsv --combine and --confirm-days 168
```

Practice effects are handled explicitly. The second SDMT is usually higher than the first, so the
reference is the best of the first two administrations. With a first-visit reference, a patient who goes
46, 52, 47 shows no decline; with the practice-aware reference the 5-point drop is visible (see the tests).

## neda

*formerly nedacalc*

**NEDA-3 and NEDA-4** for multiple sclerosis cohorts, per patient and per year, from three plain
tables: visits (EDSS), relapses, and MRI (new, enlarging and enhancing lesion counts, optional
percent brain volume change). Confirmed worsening comes from
[edssplus](https://github.com/CedricConday/msclinical#edssplus).

```bash
pip install msclinical
msclinical neda visits.tsv relapses.tsv mri.tsv --starts treatment_start.tsv --rebaseline-months 6 --years 2
```

Every component is `True` (no activity), `False` (activity) or `None` (not assessable). A window is
NEDA only when all components are `True`; a window with no MRI is `None`, never NEDA, and the cohort
summary reports assessable patients separately. That is the difference between a NEDA proportion and
a NEDA proportion that can be defended.

Definitions: Havrdova et al. 2009 (NEDA-3), Kappos et al. 2016 (NEDA-4, brain volume loss above
0.4 %/yr). Re-baselining (`--rebaseline-months`) ignores activity in the first months after a start
date, as most treatment studies do.

## relapses

*formerly relapsecheck*

Relapse tables from registries and trials carry the same errors again and again: one attack reported
twice a week apart, a fever episode logged as a relapse, an event dated before follow-up began. Every
downstream count (annualised relapse rate, RAW against PIRA, NEDA) inherits them. relapses applies the
clinical definition before any of that: at least 24 hours, no fever or infection, at least 30 days after
the previous attack, inside follow-up. It returns the clean table and every removed row with its reason.

```bash
pip install msclinical
msclinical relapses relapses.tsv --followup followup.tsv --out clean.tsv --removed removed.tsv
```

Optional columns: `fever`, `infection`, `duration_hours`. Feed `clean.tsv` to
[pirma](https://github.com/CedricConday/msclinical#pirma) or [neda](https://github.com/CedricConday/msclinical#neda).

## epochs

*formerly dmtepochs*

MS registries follow patients across several therapies. Relapse rate, NEDA and worsening computed "per
patient" mix those therapies; computed "per drug" they need epochs with an onset-of-action delay and a
carry-over after stopping. epochs builds those epochs from a treatments table and slices any dated
table into them, recording the rule that ended each epoch.

```bash
pip install msclinical
msclinical epochs treatments.tsv followup.tsv --relapses relapses.tsv --rebaseline-months 6 --carryover-days 90
```

Courses of the same drug less than 30 days apart merge. An epoch's window starts after the re-baselining
months and ends at the next drug, the stop date plus carry-over, or the end of follow-up. Untreated time is
its own epoch. Feed the windows to [neda](https://github.com/CedricConday/msclinical#neda) (`--starts`) or
[endpoints](https://github.com/CedricConday/mstrial#endpoints) (`patients.tsv` start/end, `arm` = drug).

## edssnoise


How noisy is the EDSS in your registry or trial? edssnoise estimates the probability that a visit is
rated half a point off the true score, from the visit sequences alone. A real disability level rarely
drops and returns within months; a misrating does that constantly, so the rate of local peaks and
troughs carries the misrating rate. The mapping is calibrated by simulating the same design (patients,
visit spacing, follow-up) with [sim](https://github.com/CedricConday/msclinical#sim) at known noise,
and the answer comes with a bootstrap interval over patients.

```bash
pip install msclinical
msclinical edssnoise visits.tsv          # patient, date, edss
```

Why it matters: [sim](https://github.com/CedricConday/msclinical#sim) shows that with 30 % of visits
misrated, a quarter of 12-week confirmed worsening events are noise, and they are called PIRA. Knowing
your noise level tells you which confirmation window your PIRA numbers need.

Checked on simulation (tests): 20 % misrating is recovered within 6 points. On MSLesSeg's clinical file
(115 EDSS ratings, but only 15 visits with a visit on each side) the observed reversal rate is 1 in 15,
far too few to estimate anything; registries with quarterly or yearly EDSS over several years are the
intended input. Research software. MIT.

## multiverse

*formerly msmultiverse*

How much does a PIRA share depend on the definition? multiverse runs every common combination of
disability-worsening and attribution choices on one dataset, using
[edssplus](https://github.com/CedricConday/msclinical#edssplus) and [pirma](https://github.com/CedricConday/msclinical#pirma):
12 or 24-week confirmation, roving or fixed reference, 30 or 90-day relapse window, MRI window from
reference or onset, EDSS or EDSS-plus. It reports every result and, per choice, how far flipping only
that choice moves the answer.

```bash
pip install msclinical
msclinical multiverse visits.tsv relapses.tsv --mri mri.tsv --followup-start start.tsv --out multiverse.tsv
```

### On a simulated cohort (`results/multiverse/results/`)

600 patients from [sim](https://github.com/CedricConday/msclinical#sim), 32 definitions
(`results/multiverse/results/`, regenerate with `results/multiverse/results/run.py`): **PIRA share from 0.44 to 0.86**, events 68 to 99, PIRMA
events 6 to 43. The relapse window alone moves the PIRA share by 32 points on average (30 against 90 days),
confirmation by 8, the reference by 3. The MRI window leaves the PIRA share unchanged, because the share pools
PIRMA with PIRA with MRI activity, but it moves the PIRMA count by 12 events on average (up to 23), as much as
the confirmation window does. EDSS-plus changes nothing here only because the simulation has no timed tests.

Two published cohorts that differ only in the relapse window are not reporting the same quantity. This
tool makes that visible before the comparison is made. Research software. MIT.

## tables

*formerly mstables*

One schema for MS clinical tables, and a validator that catches what breaks disability and relapse
analyses before they run:

* EDSS values off the scale (0.5 does not exist) and jumps of 3 points or more between visits;
* timed tests outside 0-600 s (wrong unit), negative or fractional lesion counts;
* unparseable dates, duplicate patient-dates;
* visits within 30 days of a relapse onset, where EDSS is not a baseline;
* patients missing from the patient table.

```bash
pip install msclinical
msclinical tables check data/        # patients.tsv, visits.tsv, relapses.tsv, mri.tsv, treatments.tsv (.csv also read)
```

The schema (`msclinical/tables/schema.py`) is the one read by [edssplus](https://github.com/CedricConday/msclinical#edssplus),
[pirma](https://github.com/CedricConday/msclinical#pirma), [neda](https://github.com/CedricConday/msclinical#neda),
[relapses](https://github.com/CedricConday/msclinical#relapses), [endpoints](https://github.com/CedricConday/mstrial#endpoints)
and written by [sim](https://github.com/CedricConday/msclinical#sim). Its first catch was in that family:
sim's rating noise could write EDSS 0.5, which is fixed.

## feasibility

*formerly msfeasibility*

"How many of our patients would qualify?" is the first question of every MS trial feasibility request, and
registries usually answer it by hand. feasibility evaluates the usual inclusion criteria on registry tables
([tables](https://github.com/CedricConday/msclinical#tables) schema): age, EDSS range at the index date, disease
duration, course, relapses or enhancing lesions in a look-back window, and stability before screening. It
prints the funnel in the order given, so the binding criterion is visible, and counts unknowns (a missing
EDSS) separately instead of dropping them silently.

```bash
pip install msclinical
msclinical feasibility registry/ --criteria '{"age":[18,55],"course":["RRMS"],"edss":[0,5.5],"activity":{"relapses":1,"gd":1},"no_relapse_days":30}'
```

## table1

*formerly mstable1*

The descriptive "Table 1" of an MS cohort paper, from the shared clinical tables
([tables](https://github.com/CedricConday/msclinical#tables) schema), per arm or group:

* age, disease duration and EDSS as median (IQR), never mean (EDSS is ordinal);
* relapses as a cohort ARR (total relapses over total patient-years), not the mean of per-patient rates;
* every row states how many patients contributed, so missing data are visible.

```bash
pip install msclinical
msclinical table1 cohort/ --group arm --markdown
```

Example on a simulated two-arm trial ([sim](https://github.com/CedricConday/msclinical#sim)): `results/table1/example.md`.

## reportparse

*formerly msreportparse*

Most MS registries have MRI **reports**, not lesion counts. reportparse reads English and German report
text and fills the structured MRI table (`new_t2`, `enlarging_t2`, `gd`). Each value is a count, 0 when the
report states absence, or empty when the report says nothing. Every value carries the sentence it came
from, so a reviewer checks evidence instead of trusting a number.

```bash
pip install msclinical
msclinical reportparse reports.tsv > mri.tsv                 # patient, date, text
msclinical reportparse --text "Zwei neue T2-hyperintense Läsionen. Keine KM-aufnehmenden Herde."
```

The output feeds [neda](https://github.com/CedricConday/msclinical#neda) (NEDA-3) and
[pirma](https://github.com/CedricConday/msclinical#pirma) (PIRMA) directly. It is rule-based on purpose: auditable,
no model, no data leaves the machine. The rules are regular expressions in `msclinical/reportparse/core.py`. The test
set is a handful of written examples, not a validation on real reports; measure it on 50 of your own before
trusting a count.

## sim

*formerly trialtruth*

Simulated multiple sclerosis trial cohorts **with every cause recorded**: relapses, relapse-associated
disability steps (RAW), relapse-independent steps (PIRA) with or without MRI activity, new T2 lesions,
and EDSS rating noise. The tables have exactly the columns
[edssplus](https://github.com/CedricConday/msclinical#edssplus), [pirma](https://github.com/CedricConday/msclinical#pirma),
[neda](https://github.com/CedricConday/msclinical#neda) and [endpoints](https://github.com/CedricConday/mstrial#endpoints)
read, so any attribution or endpoint pipeline can be scored against the truth that produced its input.

```bash
pip install msclinical
msclinical sim simulate --out sim --n 300 --years 2 --noise 0.15
```

### First result: how often is PIRA really PIRA? (`results/sim/results/pirma_vs_truth.txt`)

600 simulated patients, two years, EDSS every three months, MRI every six. Events from edssplus, attributed
by pirma, matched to the true step nearest their onset:

| EDSS misrating | Confirmation | Events | From rating noise alone | Attribution agrees with truth |
|---|---|---:|---:|---:|
| none | 12 weeks | 99 | 0 | 71 % |
| none | 24 weeks | 79 | 0 | 68 % |
| 15 % of visits | 12 weeks | 99 | 9 | 65 % |
| 15 % of visits | 24 weeks | 73 | 6 | 60 % |
| 30 % of visits | 12 weeks | 122 | 33 | 52 % |
| 30 % of visits | 24 weeks | 72 | 8 | 57 % |

The same cohort (230 true steps) under every noise level: rating noise has its own random stream, so only the
ratings differ. With a rating error on 30 % of visits (inter-rater EDSS disagreement of half a point is common),
27 % of 12-week confirmed events are rating noise, and they are called PIRA, because noise has no relapse next
to it. 24-week confirmation cuts them from 33 to 8. Even without noise, about 30 % of events get a cause
different from the true one, mostly because the step and the relapse or scan that should explain it fall in
different windows. PIRA shares reported from registries carry both errors. Regenerate with
`results/sim/results/run_pirma_vs_truth.py`.

The model is simple on purpose and its parameters are in `msclinical/sim/core.py`. Research software. MIT.
