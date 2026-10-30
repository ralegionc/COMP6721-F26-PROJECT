"""Mutational Signature Detective: our data hygiene and evaluation code.

We organised the package so that each module answers one question of the
course checklist:

    Q1 target defined          -> target.py
    Q2 no identifier leakage   -> splits.py, simulate.py
    Q3 positive class stated   -> labels.py
    Q4 majority-class check    -> metrics.py, methods.py (B0 baseline)
    Q5 same split for all      -> splits.py (save/load), storage.py (sim sets)
    Q6 where results are kept  -> storage.py
"""
