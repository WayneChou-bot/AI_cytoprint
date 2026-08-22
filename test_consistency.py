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
# ---- per-MoA breakdown: recompute Top-1/Top-5 per class from the same stored tables ----
pm = d.get("per_moa")
if pm:
    per = {}
    for x in elig:
        r = per.setdefault(x["moa"], {"n": cnt[x["moa"]], "t1": 0, "t5": 0, "q": 0})
        r["q"] += 1
        r["t1"] += int(bool(x["nb"]) and x["nb"][0]["moa"] == x["moa"])
        r["t5"] += int(any(n["moa"] == x["moa"] for n in x["nb"][:5]))
    stored = {r["moa"]: r for r in pm["rows"]}
    assert set(stored) == set(per), "per-MoA class set mismatch"
    for m, r in per.items():
        s = stored[m]
        assert s["n"] == r["n"] and s["queries"] == r["q"], f"{m}: class size mismatch"
        assert abs(s["top1"] - r["t1"]/r["q"]) < 0.002, f"{m}: per-MoA top1 mismatch"
        assert abs(s["top5"] - r["t5"]/r["q"]) < 0.002, f"{m}: per-MoA top5 mismatch"
    # the headline must be the query-weighted mean of the per-class rates it is built from
    w = sum(s["top1"]*s["queries"] for s in pm["rows"]) / sum(s["queries"] for s in pm["rows"])
    assert abs(w - met["acc1"]) < 0.002, "per-MoA rows do not average to the headline Top-1"
    s = pm["summary"]
    assert s["n_classes"] == len(pm["rows"])
    assert s["n_classes_zero_top1"] == sum(1 for r in pm["rows"] if r["top1"] == 0)
    assert s["n_classes_perfect_top1"] == sum(1 for r in pm["rows"] if r["top1"] == 1.0)
    # hits/queries must be integers that reproduce the rates, not just a percentage
    for r in pm["rows"]:
        assert abs(r["hits1"]/r["queries"] - r["top1"]) < 0.002, f"{r['moa']}: hits1/queries != top1"
        assert abs(r["hits5"]/r["queries"] - r["top5"]) < 0.002, f"{r['moa']}: hits5/queries != top5"
        assert r["top1_ci"][0] <= r["top1"] <= r["top1_ci"][1], f"{r['moa']}: Wilson CI excludes the estimate"
    # BOTH averages must reproduce — the section's whole claim is that they differ
    mac1 = sum(r["top1"] for r in pm["rows"]) / len(pm["rows"])
    mac5 = sum(r["top5"] for r in pm["rows"]) / len(pm["rows"])
    assert abs(mac1 - s["top1_macro"]) < 0.002, "macro top1 mismatch"
    assert abs(mac5 - s["top5_macro"]) < 0.002, "macro top5 mismatch"
    assert abs(s["top1_query_weighted"] - met["acc1"]) < 0.002, "query-weighted top1 != headline"
    assert s["top1_macro"] < s["top1_query_weighted"], "macro should be below query-weighted here"
    assert s["top1_hits_total"] == sum(r["hits1"] for r in pm["rows"])
    print(f"per-MoA    : {s['n_classes']} classes recomputed and consistent; "
          f"rows average to headline top1={w:.3f}")
    print(f"             top1 query-weighted={s['top1_query_weighted']} vs macro={s['top1_macro']}; "
          f"top5 {s['top5_query_weighted']} vs {s['top5_macro']}")
    # ---- activity link: recompute the headline evidence independently from the stored tables ----
    al = pm.get("activity_link")
    if al:
        c, m = al["correct_retrieved"], al["mis_retrieved"]
        hit_ok = [x for x in elig if x["nb"] and x["nb"][0]["moa"] == x["moa"]]
        hit_no = [x for x in elig if not (x["nb"] and x["nb"][0]["moa"] == x["moa"])]
        r_ok = sum(1 for x in hit_ok if x["hit"]) / len(hit_ok)
        r_no = sum(1 for x in hit_no if x["hit"]) / len(hit_no)
        assert c["n"] + m["n"] == met["n_eval"], "activity-link split does not cover the eligible set"
        assert (c["n"], m["n"]) == (len(hit_ok), len(hit_no)), "activity-link split sizes mismatch"
        assert abs(r_ok - c["candidate_active_rate"]) < 0.002, "retrieved candidate-active rate mismatch"
        assert abs(r_no - m["candidate_active_rate"]) < 0.002, "missed candidate-active rate mismatch"
        assert r_ok > r_no, "correctly-retrieved compounds should not be LESS often candidate-active"
        print(f"activity   : candidate-active {r_ok:.3f} (retrieved, n={len(hit_ok)}) "
              f"vs {r_no:.3f} (missed, n={len(hit_no)}) — recomputed, matches stored")

        # Spearman rho over the reliable classes, recomputed from the per-class rows
        rel_n = pm["summary"]["reliable_min_n"]
        sub = [r for r in pm["rows"] if r["n"] >= rel_n and "hit_rate" in r]
        stored = al.get("reliable_classes", {}).get("vs_candidate_active_rate")
        if stored and len(sub) >= 4:
            assert stored["n_classes"] == len(sub), "reliable-class count mismatch"
            def rank(v):                                    # average ranks, ties shared
                order = sorted(range(len(v)), key=lambda i: v[i])
                rk = [0.0]*len(v); i = 0
                while i < len(order):
                    j = i
                    while j+1 < len(order) and v[order[j+1]] == v[order[i]]: j += 1
                    avg = (i+j)/2 + 1
                    for k2 in range(i, j+1): rk[order[k2]] = avg
                    i = j+1
                return rk
            rx, ry = rank([r["top1"] for r in sub]), rank([r["hit_rate"] for r in sub])
            mx, my = sum(rx)/len(rx), sum(ry)/len(ry)
            num = sum((a2-mx)*(b2-my) for a2, b2 in zip(rx, ry))
            den = (sum((a2-mx)**2 for a2 in rx) * sum((b2-my)**2 for b2 in ry)) ** 0.5
            rho = num/den
            assert abs(rho - stored["spearman_rho"]) < 0.005, \
                f"Spearman rho mismatch: recomputed {rho:.3f} vs stored {stored['spearman_rho']}"
            print(f"             Spearman rho={rho:.3f} over {len(sub)} classes with n>={rel_n} "
                  f"— recomputed, matches stored")

print("OK — neighbour tables, headline metrics and the per-MoA breakdown are consistent.")
