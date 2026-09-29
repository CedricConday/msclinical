"""Column schema shared by edssplus, pirma, nedacalc, msendpoints, relapsecheck, trialtruth."""

SCHEMA = {
    "patients": {"required": {"patient": "str", "start": "date", "end": "date"},
                 "optional": {"arm": "str", "sex": "str", "birth_year": "int", "onset_date": "date", "course": "str"}},
    "visits": {"required": {"patient": "str", "date": "date"},
               "optional": {"edss": "edss", "t25fw": "seconds", "hpt9_dom": "seconds", "hpt9_nondom": "seconds", "sdmt": "int"}},
    "relapses": {"required": {"patient": "str", "date": "date"},
                 "optional": {"fever": "bool", "infection": "bool", "duration_hours": "float", "steroids": "bool"}},
    "mri": {"required": {"patient": "str", "date": "date"},
            "optional": {"new_t2": "count", "enlarging_t2": "count", "gd": "count", "pbvc": "float", "t2_volume_ml": "float"}},
    "treatments": {"required": {"patient": "str", "drug": "str", "start": "date"}, "optional": {"stop": "date"}},
}
EDSS_STEPS = {0.0, 1.0} | {x / 2 for x in range(3, 20)}  # 0, 1.0, 1.5 ... 9.5 (0.5 does not exist)
