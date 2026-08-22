#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validation test for condition_retention() — runs WITHOUT any downloaded data.

Why this file exists: the README claims the diagnostic was validated on synthetic data with
a known, injected condition effect. A claim like that is worthless unless the test is in the
repo and anyone can re-run it. This is that test.

    python mvp/test_condition_retention.py

It builds a synthetic dataset with the SAME nesting as CPJUMP1 (plate nested inside
cell_line x timepoint), injects a known condition effect, and asserts that the diagnostic:
  1. reports enrichment >> 1 when the condition structure is present (raw);
  2. falls towards 1.0 once that structure is removed by correction;
  3. sits at ~1.0 when there is no condition structure at all (negative control).
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jump_mvp import condition_retention                      # noqa: E402

TOL_PRESENT, TOL_ERASED = 1.30, 1.10                          # enrichment thresholds


def synth(n_per_plate=250, effect=3.0, dims=50, seed=0):
    """CPJUMP1-shaped: 2 cell lines x 2 timepoints x 4 plates each, plate nested in condition."""
    rng = np.random.RandomState(seed)
    X, cell, time, plate = [], [], [], []
    for ci, c in enumerate(["A549", "U2OS"]):
        for ti, tp in enumerate(["24h", "48h"]):
            centre = np.zeros(dims)
            centre[0] = effect * (1 if ci else -1)             # cell-line axis
            centre[1] = effect * (1 if ti else -1)             # timepoint axis
            for pl in range(4):
                X.append(rng.randn(n_per_plate, dims) + centre)
                cell += [c] * n_per_plate
                time += [tp] * n_per_plate
                plate += [f"{c}-{tp}-p{pl}"] * n_per_plate
    return np.vstack(X), np.array(cell), np.array(time), np.array(plate)


def strip_condition(X, cell, time):
    """Oracle 'perfect correction': centre every condition group at the origin."""
    Z = X.copy()
    for c in np.unique(cell):
        for t in np.unique(time):
            m = (cell == c) & (time == t)
            Z[m] -= Z[m].mean(0)
    return Z


def main():
    ok = True
    X, cell, time, plate = synth()

    raw_c = condition_retention(X, cell)
    raw_t = condition_retention(X, time)
    print(f"raw            cell_line={raw_c:6.2f}   timepoint={raw_t:6.2f}   (expect >> 1)")
    ok &= raw_c > TOL_PRESENT and raw_t > TOL_PRESENT

    Z = strip_condition(X, cell, time)
    cor_c = condition_retention(Z, cell)
    cor_t = condition_retention(Z, time)
    print(f"corrected      cell_line={cor_c:6.2f}   timepoint={cor_t:6.2f}   (expect ~ 1)")
    ok &= cor_c < TOL_ERASED and cor_t < TOL_ERASED

    N, ncell, ntime, _ = synth(effect=0.0, seed=1)             # negative control: no condition effect
    neg_c = condition_retention(N, ncell)
    neg_t = condition_retention(N, ntime)
    print(f"no-effect ctrl cell_line={neg_c:6.2f}   timepoint={neg_t:6.2f}   (expect ~ 1)")
    ok &= neg_c < TOL_ERASED and neg_t < TOL_ERASED

    # a label the embedding cannot possibly encode must also sit at 1.0
    rnd = np.random.RandomState(7).permutation(cell)
    shuf = condition_retention(X, rnd)
    print(f"shuffled label                       {shuf:6.2f}   (expect ~ 1)")
    ok &= shuf < TOL_ERASED

    print("\nPASS — the diagnostic detects injected condition structure and its removal."
          if ok else "\nFAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
