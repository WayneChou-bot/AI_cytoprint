#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the page's neighbour tables and headline metrics come from the same computation.
Recomputes Top-1/Top-5/random/fold from the STORED per-drug neighbour lists (third-party style)
and asserts they equal the stored `metrics` block. Run:  python test_consistency.py"""
import json, numpy as np
from pathlib import Path
from collections import Counter

# locate webdata.json whether run from package root (web/) or a flattened folder
here = Path(__file__).resolve().parent
cand = [here/"web"/"webdata.json", here/"webdata.json", Path("web/webdata.json"), Path("webdata.json")]
path = next((p for p in cand if p.exists()), cand[0])
d = json.load(open(path, encoding="utf-8"))
drugs, met = d["drugs"], d["metrics"]
ann = [x for x in drugs if x["moa"] != "unannotated"]
P = len(ann); cnt = Counter(x["moa"] for x in ann)
elig = [x for x in ann if cnt[x["moa"]] >= 2]

t1 = sum(1 for x in elig if x["nb"] and x["nb"][0]["moa"] == x["moa"]) / len(elig)
t5 = sum(1 for x in elig if any(n["moa"] == x["moa"] for n in x["nb"][:5])) / len(elig)
rnd = float(np.mean([(cnt[x["moa"]] - 1) / (P - 1) for x in elig]))
fold = round(t1 / rnd, 1)

print(f"recomputed : n_eval={len(elig)} top1={t1:.3f} top5={t5:.3f} rand={rnd:.4f} fold={fold}x")
print(f"stored     : n_eval={met['n_eval']} top1={met['acc1']} top5={met['top5_recall']} rand={met['random']} fold={met['fold']}x")

assert len(elig) == met["n_eval"], "n_eval mismatch"
assert abs(t1 - met["acc1"]) < 0.002, "top1 mismatch"
assert abs(t5 - met["top5_recall"]) < 0.002, "top5 mismatch"
assert abs(rnd - met["random"]) < 0.001, "random-baseline mismatch"
print("OK — neighbour tables and headline metrics are consistent.")
