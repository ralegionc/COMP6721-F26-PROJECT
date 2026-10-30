"""Q3. Is the positive class stated?

We have two yes/no tasks in this project:

1. Signature presence. Positive = the signature explains at least
   PRESENCE_THRESHOLD (5%) of the tumor's mutations.
2. Smoking check (E5, lung adenocarcinoma). Positive = current or reformed
   smoker. Negative = lifelong non-smoker. We exclude missing or unknown
   status (NaN) rather than guess it.
"""
import numpy as np
import pandas as pd

from .config import PRESENCE_THRESHOLD

# Plain-language definitions that we print in the run summary.
POSITIVE_CLASS = {
    "presence": f"signature weight >= {PRESENCE_THRESHOLD:.0%}",
    "smoking": "current or reformed smoker (negative: lifelong non-smoker; unknown excluded)",
}


def presence_labels(exposures, threshold=PRESENCE_THRESHOLD):
    """We return 1 where a signature is present (weight >= threshold), else 0."""
    return (np.asarray(exposures, dtype=float) >= threshold).astype(int)


# GDC field exposures.tobacco_smoking_status (legacy TCGA name:
# tobacco_smoking_history). The GDC currently returns text; older exports use
# codes 1-7:
#   1 lifelong non-smoker, 2 current smoker, 3 reformed > 15 years,
#   4 reformed <= 15 years, 5 reformed (duration unknown),
#   6 smoker at diagnosis, 7 smoking history not documented.
# clinical.unmapped_smoking_values() lists any value we do not recognise.
# Each entry maps to (label, quit more than 15 years ago).
_TEXT_MAP = {
    "lifelong non-smoker": (0, False),
    "current smoker": (1, False),
    "current reformed smoker for > 15 yrs": (1, True),
    "current reformed smoker for < or = 15 yrs": (1, False),
    "current reformed smoker, duration not specified": (1, False),
    "smoker at diagnosis": (1, False),
}
# Code 7 ("not documented") is deliberately absent, so it becomes NaN.
_CODE_MAP = {1: (0, False), 2: (1, False), 3: (1, True), 4: (1, False), 5: (1, False), 6: (1, False)}


def smoking_labels(values, exclude_long_reformed=False):
    """We map smoking history to 1 (smoker), 0 (never smoked) or NaN (excluded).

    With exclude_long_reformed=True we also drop people who quit more than
    15 years ago, for our sensitivity analysis.
    """
    out = []
    for v in pd.Series(values):
        label, long_reformed = None, False
        # First we try the text form; failing that, a number stored as text ("3").
        if isinstance(v, str):
            key = v.strip().lower()
            if key in _TEXT_MAP:
                label, long_reformed = _TEXT_MAP[key]
            else:
                try:
                    v = float(key)
                except ValueError:
                    v = None
        # Then the numeric code form; we only accept whole numbers we know.
        if label is None and isinstance(v, (int, float, np.integer, np.floating)) and not pd.isna(v):
            code = int(v)
            if code in _CODE_MAP and code == v:
                label, long_reformed = _CODE_MAP[code]
        if label is None or (exclude_long_reformed and long_reformed):
            out.append(np.nan)
        else:
            out.append(float(label))
    return pd.Series(out, index=pd.Series(values).index, name="smoker")
