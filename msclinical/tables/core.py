"""Validate a set of MS tables against the schema and against the checks that matter downstream."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .schema import EDSS_STEPS, SCHEMA


def load(folder) -> dict:
    out = {}
    for name in SCHEMA:
        for ext, sep in ((".tsv", "\t"), (".csv", ",")):
            p = Path(folder) / f"{name}{ext}"
            if p.exists():
                out[name] = pd.read_csv(p, sep=sep)
    return out


def validate(tables: dict) -> list[dict]:
    issues = []
    def add(table, level, msg, rows=None):
        issues.append({"table": table, "level": level, "message": msg, "rows": rows or []})
    for name, df in tables.items():
        spec = SCHEMA.get(name)
        if spec is None:
            add(name, "warning", "unknown table name")
            continue
        miss = [c for c in spec["required"] if c not in df]
        if miss:
            add(name, "error", f"missing required columns {miss}")
            continue
        extra = [c for c in df if c not in spec["required"] and c not in spec["optional"]]
        if extra:
            add(name, "info", f"columns not in the schema (kept, ignored by the tools): {extra}")
        for col, typ in {**spec["required"], **spec["optional"]}.items():
            if col not in df:
                continue
            s = df[col]
            if typ == "date":
                bad = pd.to_datetime(s, errors="coerce").isna() & s.notna()
                if bad.any():
                    add(name, "error", f"{col}: unparseable dates", list(df.index[bad][:20]))
            elif typ == "edss":
                v = pd.to_numeric(s, errors="coerce")
                nonnum = v.isna() & s.notna()
                if nonnum.any():
                    add(name, "error", f"{col}: not a number (decimal comma?)", list(df.index[nonnum][:20]))
                bad = v.notna() & ~v.isin(EDSS_STEPS | {10.0})
                if bad.any():
                    add(name, "error", f"{col}: values not on the EDSS scale (0, 1.0, 1.5 ... 10; 0.5 does not exist)", list(df.index[bad][:20]))
            elif typ == "seconds":
                v = pd.to_numeric(s, errors="coerce")
                bad = v.notna() & ((v <= 0) | (v > 600))
                if bad.any():
                    add(name, "warning", f"{col}: outside 0-600 s (unit? 180 s is the usual test cap)", list(df.index[bad][:20]))
            elif typ == "count":
                v = pd.to_numeric(s, errors="coerce")
                bad = v.notna() & ((v < 0) | (v % 1 != 0))
                if bad.any():
                    add(name, "error", f"{col}: counts must be non-negative integers", list(df.index[bad][:20]))
        if "patient" in df and "date" in df:
            dup = df.duplicated(["patient", "date"], keep=False)
            if dup.any() and name in ("visits", "mri"):
                add(name, "warning", "several rows for one patient and date", list(df.index[dup][:20]))
    v = tables.get("visits")
    if v is not None and {"patient", "date", "edss"} <= set(v):
        w = v.assign(date=pd.to_datetime(v["date"], errors="coerce"), edss=pd.to_numeric(v["edss"], errors="coerce")).sort_values(["patient", "date"])
        jump = w.groupby("patient")["edss"].diff().abs() >= 3
        if jump.any():
            add("visits", "warning", "EDSS changes by 3 points or more between consecutive visits", list(w.index[jump][:20]))
    r, p = tables.get("relapses"), tables.get("patients")
    if v is not None and r is not None and {"patient", "date"} <= set(r):
        rr = r.assign(date=pd.to_datetime(r["date"], errors="coerce"))
        vv = v.assign(date=pd.to_datetime(v["date"], errors="coerce")).reset_index(names="row")
        m = vv.merge(rr, on="patient", suffixes=("", "_rel"))
        during = m[(m["date"] >= m["date_rel"]) & ((m["date"] - m["date_rel"]).dt.days < 30)]
        if len(during):
            add("visits", "warning", f"{len(during)} visits within 30 days after a relapse onset: EDSS there is not a baseline "
                "(edssplus reference choice and pirma RAW window should see them)", sorted(set(during["row"]))[:20])
    for name in ("visits", "relapses", "mri"):
        t = tables.get(name)
        if t is None or p is None or not {"patient", "start", "end"} <= set(p):
            continue
        unknown = set(t["patient"]) - set(p["patient"])
        if unknown:
            add(name, "error", f"patients not in patients table: {sorted(unknown)[:10]}")
    return issues
