#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MorphoProfile — authoritative, reproducible pipeline (single source of truth).

Everything the web page shows is produced HERE, offline, from the SAME
standardized representation, so the neighbour tables and the headline metrics
are guaranteed to agree. Run:

    pip install -r requirements.txt
    python pipeline.py            # -> web/webdata.json
    python build_page.py          # -> web/MorphoProfile.html

Data: LINCS Cell Painting pilot (Cell Painting Gallery cpg0004), A549,
batch 2016_04_01_a549_48hr_batch1. Per-plate CellProfiler profiles are hosted
on GitHub (no AWS account needed). Pinned to COMMIT below for reproducibility.
The same pipeline points at JUMP cpg0016 by swapping the loader (see run_on_JUMP.py).

Metric definitions (all computed on standardized per-compound consensus profiles):
  * neighbours  : cosine similarity, top-5, on standardized features.
  * MoA top-1   : nearest neighbour shares the query's (first) MoA.
  * MoA top-5   : any of the 5 nearest neighbours shares it.
  * random top-1: expected value under random ranking, EXCLUDING self, matched
                  to the eligible-query MoA-class sizes:  mean_i (size(MoA_i)-1)/(P-1),
                  where P = number of MoA-annotated compounds (candidate pool).
  * candidate-active ("hit"): consensus L2 distance from the DMSO centre exceeds
                  the 95th percentile of a DMSO-only null (random 5-well medians).
                  This is a heuristic screen, NOT a significance test (no plate-
                  matched null, no whitening, no per-compound FDR). Labelled
                  "candidate active" in the UI accordingly.
  * de-biasing mAP: information-retrieval mAP for same-MoA retrieval, computed
                  identically for raw / sphered / sphered+Harmony, with a
                  compound-level bootstrap 95% CI on each and on the delta.
"""
import json, requests, warnings, io
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
warnings.filterwarnings("ignore")
RNG = np.random.RandomState(42)

COMMIT = "da8ae6a3bc103346095d61b4ee02f08fc85a5d98"          # lincs-cell-painting pin
RAW    = f"https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/{COMMIT}"
MEDIA  = f"https://media.githubusercontent.com/media/broadinstitute/lincs-cell-painting/{COMMIT}"
BATCH  = "2016_04_01_a549_48hr_batch1"
N_MAPS = 10
DATA = Path("data"); WEB = Path("web"); DATA.mkdir(exist_ok=True); WEB.mkdir(exist_ok=True)

# ─────────────────────────── data loading (cached) ───────────────────────────
def load_wells():
    cache = DATA/"lincs_norm_10maps.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    import concurrent.futures as cf
    bc = pd.read_csv(f"{RAW}/metadata/platemaps/{BATCH}/barcode_platemap.csv")
    grp = bc.groupby("Plate_Map_Name")["Assay_Plate_Barcode"].apply(list)
    plates = [p for m in list(grp.index)[:N_MAPS] for p in grp[m]]
    def get(p):
        d = pd.read_csv(f"{MEDIA}/profiles/{BATCH}/{p}/{p}_normalized.csv.gz")
        d["Metadata_Assay_Plate_Barcode"] = p; return d
    with cf.ThreadPoolExecutor(max_workers=10) as ex:
        frames = list(ex.map(get, plates))
    meta = sorted({c for f in frames for c in f.columns if c.startswith("Metadata")})
    feat = [c for c in frames[0].columns if not c.startswith("Metadata")]
    for f in frames:
        for m in meta:
            if m not in f.columns: f[m] = np.nan
    df = pd.concat([f[meta+feat] for f in frames], ignore_index=True)
    df[feat] = df[feat].astype("float32"); df.to_parquet(cache)
    return df

MAXABS_CAP = 100.0   # features whose max |value| exceeds this are MAD-blowup artifacts → DROP
WINSOR = 15.0        # winsorize legit heavy tails to ±WINSOR
VAR_MIN, CORR_MAX = 1e-3, 0.9

def sanitize_select(df, feat):
    """SANITIZE FIRST, then select. Returns (sel, clean_frame) where clean_frame holds the
    CLEANED feature values (finite, |x|<=WINSOR) that ALL downstream steps must use.
    Fixes the P0 bug where clipping was applied only to pick column names, not the analysis matrix."""
    X = df[feat].replace([np.inf, -np.inf], np.nan)
    maxabs = X.abs().max()
    broken = list(maxabs[maxabs > MAXABS_CAP].index)          # e.g. near-zero-MAD blowups (values up to 1e19)
    X = X.drop(columns=broken).dropna(axis=1)                 # drop broken + all-NaN features
    X = X.clip(-WINSOR, WINSOR)                               # winsorize remaining tails, WRITTEN BACK
    X = X[X.var()[lambda s: s > VAR_MIN].index]               # variance filter on CLEAN data
    corr = np.corrcoef(X.values.T); cols = list(X.columns); drop = set()
    for i in range(len(cols)):
        if i in drop: continue
        for j in np.where(np.abs(corr[i]) > CORR_MAX)[0]:
            if j > i: drop.add(j)
    sel = [cols[i] for i in range(len(cols)) if i not in drop]
    Xc = X[sel]
    assert np.isfinite(Xc.values).all(), "sanitize failed: non-finite remain"
    assert np.abs(Xc.values).max() <= WINSOR + 1e-6, "sanitize failed: values exceed winsor cap"
    meta = [c for c in df.columns if c.startswith("Metadata")]
    clean = df[meta].join(Xc)
    print(f"sanitize: dropped {len(broken)} blowup features (max|x|>{MAXABS_CAP}); "
          f"selected {len(sel)}; clean max|x|={np.abs(Xc.values).max():.1f}")
    return sel, clean

def moa_name_map():
    nm = pd.read_csv(io.StringIO(requests.get(
        f"{RAW}/metadata/moa/repurposing_info_external_moa_map_resolved.tsv").text), sep="\t")
    return dict(zip(nm.broad_sample, nm.pert_iname))

# ─────────────────────────── metrics ───────────────────────────
def average_precision(rels):
    rels = np.asarray(rels, float); npos = rels.sum()
    if npos == 0: return np.nan
    cum = np.cumsum(rels); prec = cum/np.arange(1, len(rels)+1)
    return float((prec*rels).sum()/npos)

def moa_map(Z, moa):
    """information-retrieval mAP for same-MoA retrieval over the annotated pool."""
    sim = cosine_similarity(Z); np.fill_diagonal(sim, -np.inf)
    pool = [i for i in range(len(moa)) if isinstance(moa[i], str)]
    cnt = Counter(moa[i] for i in pool)
    elig = [i for i in pool if cnt[moa[i]] >= 2]
    aps = []
    for i in elig:
        order = [j for j in np.argsort(-sim[i]) if j in set(pool) and j != i]
        rels = [1 if moa[j] == moa[i] else 0 for j in order]
        aps.append(average_precision(rels))
    return float(np.nanmean(aps)), aps, elig

def retrieval_metrics(Z, moa):
    sim = cosine_similarity(Z); np.fill_diagonal(sim, -np.inf)
    pool = [i for i in range(len(moa)) if isinstance(moa[i], str)]
    P = len(pool); cnt = Counter(moa[i] for i in pool)
    elig = [i for i in pool if cnt[moa[i]] >= 2]
    t1 = t5 = 0
    for i in elig:
        order = [j for j in np.argsort(-sim[i]) if j in set(pool) and j != i]
        if moa[order[0]] == moa[i]: t1 += 1
        if any(moa[order[k]] == moa[i] for k in range(min(5, len(order)))): t5 += 1
    # random top-1 EXCLUDING self, matched to eligible-query class sizes
    rnd = float(np.mean([(cnt[moa[i]]-1)/(P-1) for i in elig]))
    return dict(acc1=round(t1/len(elig), 3), top5_recall=round(t5/len(elig), 3),
                random=round(rnd, 4), fold=round((t1/len(elig))/rnd, 1), n_eval=len(elig))

def wilson(k, n, z=1.96):
    """Wilson score interval for a binomial rate — behaves sensibly at 0/n and n/n,
    which a normal-approximation interval does not, and most classes here are tiny."""
    if n == 0: return [None, None]
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return [round(max(0.0, c-h), 3), round(min(1.0, c+h), 3)]

def per_moa_breakdown(Z, moa, hit=None, act=None):
    """Per-MoA retrieval breakdown.

    A single headline Top-1 is a QUERY-WEIGHTED average over MoA classes with wildly
    different difficulty and size, so it hides the shape of the result: a few classes with
    strong, distinctive morphology are retrieved almost perfectly, most sit at zero. The
    class-balanced (macro) average of the same data is a different number, and both are
    reported here — quoting one alone is choosing whichever definition flatters the result.

    If `hit` (candidate-active flag) and `act` (activity score) are supplied, each class also
    carries its candidate-active rate and median activity, which is what lets the page ask
    whether retrieval failure tracks *phenotype strength* rather than asserting it.

    Returns (rows, summary). Each row carries hits/queries, not just the percentage, plus a
    Wilson 95% interval — 57 of these classes hold only 2 compounds, where a rate is 0 or 1.
    """
    sim = cosine_similarity(Z); np.fill_diagonal(sim, -np.inf)
    pool = [i for i in range(len(moa)) if isinstance(moa[i], str)]
    poolset = set(pool); P = len(pool); cnt = Counter(moa[i] for i in pool)
    elig = [i for i in pool if cnt[moa[i]] >= 2]
    per = {}
    for i in elig:
        order = [j for j in np.argsort(-sim[i]) if j in poolset and j != i]
        rels = [1 if moa[j] == moa[i] else 0 for j in order]
        d = per.setdefault(moa[i], {"n": cnt[moa[i]], "t1": 0, "t5": 0, "ap": [], "q": 0,
                                    "hit": 0, "act": []})
        d["q"] += 1
        d["t1"] += int(rels[0] == 1)
        d["t5"] += int(any(rels[:5]))
        d["ap"].append(average_precision(rels))
        if hit is not None: d["hit"] += int(bool(hit[i]))
        if act is not None: d["act"].append(float(act[i]))
    rows = []
    for m, d in per.items():
        r = {"moa": m, "n": d["n"], "queries": d["q"], "hits1": d["t1"], "hits5": d["t5"],
             "top1": round(d["t1"]/d["q"], 3), "top5": round(d["t5"]/d["q"], 3),
             "top1_ci": wilson(d["t1"], d["q"]),
             "mAP": round(float(np.nanmean(d["ap"])), 3),
             "rand": round((d["n"]-1)/(P-1), 4)}
        if hit is not None: r["hit_rate"] = round(d["hit"]/d["q"], 3)
        if act is not None: r["median_act"] = round(float(np.median(d["act"])), 1)
        rows.append(r)
    rows.sort(key=lambda r: (-r["top1"], -r["mAP"], -r["n"]))

    RELIABLE = 5                                     # a rate from <5 compounds is nearly uninformative
    rel = [r for r in rows if r["n"] >= RELIABLE]
    tot1 = sum(r["hits1"] for r in rows)
    summary = {
        "n_classes": len(rows), "n_queries": len(elig),
        "n_classes_zero_top1": sum(1 for r in rows if r["top1"] == 0),
        "n_classes_perfect_top1": sum(1 for r in rows if r["top1"] == 1.0),
        "n_classes_size2": sum(1 for r in rows if r["n"] == 2),
        "reliable_min_n": RELIABLE, "n_classes_reliable": len(rel),
        # the two averages of the SAME data — the whole point of this section
        "top1_query_weighted": round(tot1/len(elig), 3),
        "top1_macro": round(float(np.mean([r["top1"] for r in rows])), 3),
        "top5_query_weighted": round(sum(r["hits5"] for r in rows)/len(elig), 3),
        "top5_macro": round(float(np.mean([r["top5"] for r in rows])), 3),
        "mAP_query_weighted": round(float(np.average([r["mAP"] for r in rows],
                                                    weights=[r["queries"] for r in rows])), 3),
        "mAP_macro_per_class": round(float(np.mean([r["mAP"] for r in rows])), 3),
        "top1_hits_total": tot1,
        "top1_share_from_top5_classes": round(sum(r["hits1"] for r in rows[:5])/max(1, tot1), 3),
        "top1_share_from_top5_classes_reliable": round(
            sum(r["hits1"] for r in rel[:5])/max(1, sum(r["hits1"] for r in rel)), 3) if rel else None,
        "hits_from_small_classes": sum(r["hits1"] for r in rows[:5] if r["n"] < RELIABLE),
    }
    return rows, summary

def activity_link(rows, drugs_hit, drugs_t1, min_n=5):
    """Does retrieval failure track PHENOTYPE STRENGTH, or is it unexplained?

    This is the layer that separates "these classes were not retrieved" (which the breakdown
    shows) from "these classes produced no detectable morphological signal" (which it does
    NOT show on its own). It is descriptive support, not independent biological validation:
    the activity score is derived from the SAME morphology features as the retrieval, so a
    correlation between them cannot prove the compounds are biologically inert.
    """
    from scipy import stats
    out = {"caveat": ("activity is derived from the same morphology features as retrieval, "
                      "so this is descriptive support, not independent biological validation")}
    n1 = sum(1 for h, t in zip(drugs_hit, drugs_t1) if t)
    out["correct_retrieved"] = {"n": n1,
        "candidate_active_rate": round(sum(h for h, t in zip(drugs_hit, drugs_t1) if t)/max(1, n1), 3)}
    n0 = sum(1 for t in drugs_t1 if not t)
    out["mis_retrieved"] = {"n": n0,
        "candidate_active_rate": round(sum(h for h, t in zip(drugs_hit, drugs_t1) if not t)/max(1, n0), 3)}
    for tag, sub in (("all_classes", rows), ("reliable_classes", [r for r in rows if r["n"] >= min_n])):
        if len(sub) < 4 or "hit_rate" not in sub[0]: continue
        x = [r["top1"] for r in sub]
        for key, y in (("vs_candidate_active_rate", [r["hit_rate"] for r in sub]),
                       ("vs_median_activity", [r["median_act"] for r in sub])):
            s = stats.spearmanr(x, y)
            out.setdefault(tag, {})[key] = {"n_classes": len(sub),
                                            "spearman_rho": round(float(s.statistic), 3),
                                            "p": float(f"{s.pvalue:.3g}")}
    return out

# ─────────────────────────── de-biasing (raw / sphered / +harmony) ───────────
def zca(Xt, dmso):
    mu = dmso.mean(0)
    cov = np.cov((dmso-mu).T) + 1e-3*np.eye(dmso.shape[1])
    U, S, _ = np.linalg.svd(cov)
    return (Xt-mu) @ (U @ np.diag(1/np.sqrt(S)) @ U.T)

def consensus(Zwell, bs):
    fc = [f"f{i}" for i in range(Zwell.shape[1])]
    d = pd.DataFrame(np.nan_to_num(Zwell), columns=fc); d["bs"] = bs
    g = d.groupby("bs")[fc].median()
    return g.index.tolist(), g.values

def bootstrap_map(Z, moa, n=150):
    _, aps, _ = moa_map(Z, moa)
    aps = np.array(aps);
    boots = [np.nanmean(RNG.choice(aps, len(aps), replace=True)) for _ in range(n)]
    return float(np.nanmean(aps)), (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))), aps

def debias(df, sel):
    trt = df[df.Metadata_broad_sample != "DMSO"]
    bs = trt["Metadata_broad_sample"].values
    moa_of = trt.groupby("Metadata_broad_sample")["moa1"].first()
    Xt = trt[sel].values; dmso = df[df.Metadata_broad_sample == "DMSO"][sel].values
    reps = {"raw": StandardScaler().fit_transform(Xt), "sphered": zca(Xt, dmso)}
    out = {"mAP": {}, "ci": {}, "aps": {}}
    for name in ["raw", "sphered"]:
        ids, V = consensus(reps[name], bs)
        moa = moa_of.reindex(ids).values
        m, ci, aps = bootstrap_map(V, moa)
        out["mAP"][name] = round(m, 3); out["ci"][name] = [round(ci[0], 3), round(ci[1], 3)]
        out["aps"][name] = aps
    # +Harmony at (compound,plate) consensus, batch = plate (light on memory)
    try:
        import harmonypy
        sph = reps["sphered"]; fc = [f"f{i}" for i in range(sph.shape[1])]
        d = pd.DataFrame(sph, columns=fc); d["bs"] = bs; d["plate"] = trt["plate"].values
        cp = d.groupby(["bs", "plate"])[fc].median().reset_index()
        ho = harmonypy.run_harmony(cp[fc].values.astype("float32"), cp[["plate"]], ["plate"], max_iter_harmony=10)
        Zc = np.asarray(ho.Z_corr); Zc = Zc.T if Zc.shape[0] != len(cp) else Zc
        corr = pd.DataFrame(Zc, columns=fc); corr["bs"] = cp["bs"].values
        g = corr.groupby("bs")[fc].median()
        m, ci, _ = bootstrap_map(g.values, moa_of.reindex(g.index).values)
        out["mAP"]["sphered+harmony"] = round(m, 3); out["ci"]["sphered+harmony"] = [round(ci[0], 3), round(ci[1], 3)]
    except Exception as e:
        print("harmony skipped:", repr(e)[:100])
    # paired bootstrap delta (sphered - raw) over the shared eligible queries
    ar, asf = out["aps"]["raw"], out["aps"]["sphered"]
    L = min(len(ar), len(asf)); ar, asf = np.array(ar[:L]), np.array(asf[:L])
    deltas = []
    for _ in range(2000):
        idx = RNG.randint(0, L, L); deltas.append(np.nanmean(asf[idx]) - np.nanmean(ar[idx]))
    out["delta"] = {"value": round(float(np.nanmean(asf) - np.nanmean(ar)), 3),
                    "ci": [round(float(np.percentile(deltas, 2.5)), 3), round(float(np.percentile(deltas, 97.5)), 3)],
                    "p_gt0": round(float(np.mean(np.array(deltas) > 0)), 3)}
    out.pop("aps")
    return out

# ─────────────────────────── main ───────────────────────────
DERIVED_VERSION = "v2-sanitized"   # bump when sanitize rule changes → invalidates derived cache

def main():
    raw = load_wells()
    feat_raw = [c for c in raw.columns if not c.startswith("Metadata")]
    # P0 FIX: sanitize BEFORE selection; `clean` carries the CLEANED values used everywhere.
    sel, clean = sanitize_select(raw, feat_raw)
    name = moa_name_map()
    clean["moa1"] = clean["Metadata_moa"].map(lambda m: m.split("|")[0].strip() if isinstance(m, str) else None)
    clean["plate"] = clean["Metadata_Assay_Plate_Barcode"]
    clean["drug"] = clean["Metadata_broad_sample"].map(lambda b: name.get(b) if isinstance(name.get(b), str) and name.get(b) else b)
    # versioned cache + manifest so downstream/JUMP never read a stale or contaminated matrix
    clean.to_parquet(DATA/f"lincs_selected_{DERIVED_VERSION}.parquet")
    json.dump({"derived_version": DERIVED_VERSION, "lincs_commit": COMMIT, "batch": BATCH,
               "plate_maps": N_MAPS, "n_plates": int(clean["Metadata_Assay_Plate_Barcode"].nunique()),
               "sanitize": {"drop_maxabs_gt": MAXABS_CAP, "winsorize": WINSOR, "var_gt": VAR_MIN, "corr_lt": CORR_MAX},
               "features_raw": len(feat_raw), "features_selected": len(sel),
               "clean_maxabs": float(np.abs(clean[sel].values).max())},
              open(DATA/"manifest.json", "w", encoding="utf-8"), indent=2)

    df = clean   # everything below reads the cleaned matrix
    dmso = df[df.Metadata_broad_sample == "DMSO"][sel].values
    trt = df[df.Metadata_broad_sample != "DMSO"]
    cons = trt.groupby("drug")[sel].median()
    moa = trt.groupby("drug")["moa1"].first().reindex(cons.index).values
    ids = list(cons.index)

    # candidate-active screen (heuristic; NOT a significance test)
    mu = np.median(dmso, axis=0); act = np.linalg.norm(cons.values - mu, axis=1)
    null = [np.linalg.norm(np.median(dmso[RNG.choice(len(dmso), 5, replace=False)], axis=0) - mu) for _ in range(5000)]
    thr = float(np.percentile(null, 95)); hit = act > thr

    # ONE standardized representation drives neighbours + metrics + UMAP
    Z = StandardScaler().fit_transform(cons.values)
    sim = cosine_similarity(Z); np.fill_diagonal(sim, -np.inf)
    met = retrieval_metrics(Z, moa)

    moa_rows, moa_summary = per_moa_breakdown(Z, moa, hit=hit, act=act)
    _pool = [i for i in range(len(moa)) if isinstance(moa[i], str)]
    _cnt = Counter(moa[i] for i in _pool)
    _elig = [i for i in _pool if _cnt[moa[i]] >= 2]
    _sim = cosine_similarity(Z); np.fill_diagonal(_sim, -np.inf)
    _t1, _hit = [], []
    for i in _elig:
        _o = [j for j in np.argsort(-_sim[i]) if j in set(_pool) and j != i]
        _t1.append(int(moa[_o[0]] == moa[i])); _hit.append(int(bool(hit[i])))
    moa_activity = activity_link(moa_rows, _hit, _t1, min_n=moa_summary["reliable_min_n"])
    import umap
    emb = umap.UMAP(n_neighbors=15, min_dist=0.4, metric="cosine", random_state=42).fit_transform(Z)
    top_moa = [m for m, _ in Counter([x for x in moa if isinstance(x, str)]).most_common(10)]
    grp = lambda m: m if (isinstance(m, str) and m in top_moa) else "other"

    def ftype(c):
        for t in ["AreaShape", "Intensity", "Texture", "Granularity", "RadialDistribution", "Correlation", "Neighbors"]:
            if t in c: return t
        return "Other"
    fp_types = ["AreaShape", "Intensity", "Texture", "Granularity", "RadialDistribution", "Correlation", "Neighbors", "Other"]
    fts = np.array([ftype(c) for c in sel])

    drugs = []
    for i, d in enumerate(ids):
        order = [j for j in np.argsort(-sim[i]) if isinstance(moa[j], str)][:5]
        nb = [{"name": ids[j], "sim": round(float(sim[i][j]), 2), "moa": str(moa[j]), "match": bool(moa[j] == moa[i])} for j in order]
        votes = Counter(x["moa"] for x in nb[:3]).most_common(1)[0] if nb else ("?", 0)
        row = cons.iloc[i].values
        fp = [round(float(row[fts == t].mean()), 2) if (fts == t).any() else 0.0 for t in fp_types]
        drugs.append({"name": d, "moa": (moa[i] if isinstance(moa[i], str) else "unannotated"),
                      "x": round(float(emb[i, 0]), 2), "y": round(float(emb[i, 1]), 2),
                      "grp": grp(moa[i]), "hit": bool(hit[i]), "act": round(float(act[i]), 1),
                      "nb": nb, "pred": votes[0], "votes": f"{votes[1]}/3", "fp": fp})

    out = {
      "provenance": {"demo_dataset": f"LINCS Cell Painting pilot (cpg0004) · {BATCH}",
        "commit": COMMIT, "cell_line": "A549", "plates": N_MAPS*5, "plate_maps": N_MAPS,
        "wells": int(df.shape[0]), "dmso_wells": int((df.Metadata_broad_sample == 'DMSO').sum()),
        "drugs": int(len(ids)), "features_raw": len(feat_raw), "features_selected": len(sel),
        "sanitized": True, "clean_maxabs": round(float(np.abs(cons.values).max()), 1),
        "target_dataset": "JUMP-CP cpg0016 · ~116k compounds + gene perturbations (~15k unique genes) · U2OS"},
      "metrics": {"hit_frac": round(float(hit.mean()), 3), "n_hits": int(hit.sum()),
        "acc1": met["acc1"], "top5_recall": met["top5_recall"], "random": met["random"],
        "fold": met["fold"], "n_eval": met["n_eval"], "thr": round(thr, 1),
        "hit_metric": "candidate-active: consensus L2 distance from DMSO centre > DMSO-null P95 (heuristic screen, not a significance test)"},
      "per_moa": {"rows": moa_rows, "summary": moa_summary, "activity_link": moa_activity},
      "classes": top_moa,
      "moa_sizes": dict(Counter([x for x in moa if isinstance(x, str)]).most_common(12)),
      "drugs": drugs, "fp_types": fp_types,
      "activity_null": [round(float(x), 1) for x in RNG.choice(null, 150)],
      "debias": debias(df, sel),
    }
    json.dump(out, open(WEB/"webdata.json", "w", encoding="utf-8"))
    dd = out["debias"]
    ms = out["per_moa"]["summary"]
    print(f"per-MoA: {ms['n_classes']} classes  perfect-top1={ms['n_classes_perfect_top1']}  "
          f"zero-top1={ms['n_classes_zero_top1']}")
    print(f"  top1 query-weighted={ms['top1_query_weighted']} macro={ms['top1_macro']}  |  "
          f"top5 qw={ms['top5_query_weighted']} macro={ms['top5_macro']}  |  "
          f"mAP qw={ms['mAP_query_weighted']} macro={ms['mAP_macro_per_class']}")
    al = out["per_moa"]["activity_link"]
    print(f"  candidate-active rate: correctly retrieved {al['correct_retrieved']['candidate_active_rate']} "
          f"(n={al['correct_retrieved']['n']}) vs mis-retrieved {al['mis_retrieved']['candidate_active_rate']} "
          f"(n={al['mis_retrieved']['n']})")
    rc = al.get("reliable_classes", {}).get("vs_candidate_active_rate", {})
    if rc: print(f"  Spearman top1~candidate-active over n>={ms['reliable_min_n']} classes: "
                 f"rho={rc['spearman_rho']} p={rc['p']} (k={rc['n_classes']})")
    print(f"drugs={len(ids)}  hit={100*out['metrics']['hit_frac']:.0f}%  "
          f"top1={100*met['acc1']:.1f}% ({met['fold']}x rand={100*met['random']:.2f}%)  top5={100*met['top5_recall']:.1f}%")
    print(f"debias mAP raw={dd['mAP']['raw']} {dd['ci']['raw']}  sphered={dd['mAP']['sphered']} {dd['ci']['sphered']}  "
          f"+harmony={dd['mAP'].get('sphered+harmony','NA')}")
    print(f"delta(sphered-raw)={dd['delta']['value']} CI{dd['delta']['ci']} P(>0)={dd['delta']['p_gt0']}")
    print("wrote web/webdata.json")

if __name__ == "__main__":
    main()
