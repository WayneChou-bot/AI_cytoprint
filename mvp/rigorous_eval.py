#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rigorous batch-correction evaluation engine (the "next rung").

Implements the reviewer's asks:
  (1) SAME-SAMPLE comparison of representations in a common reduced space
      (raw / sphered / sphered+Harmony), so differences are the method, not the inputs.
  (2) DUAL metrics — you must measure BOTH:
        · biology retention  = MoA-retrieval mAP  (higher = better; keeps signal)
        · batch removal      = plate-mixing enrichment (→1 = well mixed; >1 = batch structure remains)
      Reporting only one is how people fool themselves: Harmony optimises mixing BY
      CONSTRUCTION, so the real question is what it costs biology.
  (3) PLATE-aware BLOCK bootstrap for 95% CI — resample whole plates (the technical
      unit / block), NOT individual compounds, so the CI reflects batch-level variability.

POINT ESTIMATES vs BOOTSTRAP (important, and a bug that was fixed here):
  mAP is a non-linear function of the *pool of compounds present*. A block bootstrap changes
  that pool in every draw, so the mean of the bootstrap draws is a BIASED estimate of the
  full-data value — on LINCS the mean-of-draws delta for sphered+Harmony was -0.001 while the
  observed difference was +0.008, i.e. the reported sign was wrong. Therefore:
      · every `*_mAP` / `batch_mixing` value is the OBSERVED statistic on the full data;
      · the bootstrap is used ONLY to estimate block-level DISPERSION, and intervals are
        recentred on the observed value (`observed + (quantile - mean)` of the draws);
      · `bootstrap_mean` is reported alongside so the bias is visible, not hidden.

BLOCK RESOLUTION (a hard limit, disclosed rather than papered over):
  a block bootstrap over k plates has only C(2k-1, k) distinct resamples — 35 for k=4.
  Raising B past that refines the weighting, not the resolution. Every result JSON carries
  `n_distinct_block_draws`; with k=4 (the CPJUMP1 strata) treat intervals as EXPLORATORY.

Runs on REAL data. Here the block is the PLATE: both datasets used in this repo — LINCS
batch1 and the CPJUMP1 pilot — are single-source, so plate is the available technical unit.
Genuinely multi-source correction (block = data-generating SOURCE) is future work; nothing
in this repo demonstrates it.

Caveat (stated, not hidden): correction transforms are fit ONCE on the full data; the
bootstrap resamples the EVALUATION over plates (given fixed embeddings). Re-fitting the
correction inside every bootstrap draw is the gold standard but far slower.
"""
import json, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from collections import Counter
warnings.filterwarnings("ignore")
RNG = np.random.RandomState(42)
NPC, KMIX, SUBS, B = 50, 20, 3000, 1000  # PCA dims, kNN for mixing, subsample size, bootstrap draws

def average_precision(rels):
    rels = np.asarray(rels, float); npos = rels.sum()
    if npos == 0: return np.nan
    cum = np.cumsum(rels); return float(((cum/np.arange(1, len(rels)+1))*rels).sum()/npos)

def moa_map_consensus(emb, comp, moa):
    """biology: per-compound consensus in `emb`, mAP for same-MoA retrieval.
    Vectorised; verified to reproduce the original per-query loop exactly (diff = 0)."""
    df = pd.DataFrame(emb); df["c"] = comp
    cons = df.groupby("c").median(); ids = list(cons.index)
    m = pd.Series(moa).groupby(comp).first().reindex(ids).values
    V = cons.values
    pool = np.array([i for i in range(len(ids)) if isinstance(m[i], str)])
    if len(pool) < 2: return np.nan
    Vp = V[pool]; lab = np.asarray([m[i] for i in pool])
    norm = np.linalg.norm(Vp, axis=1); norm[norm == 0] = 1.0
    Vn = Vp / norm[:, None]
    sim = Vn @ Vn.T; np.fill_diagonal(sim, -np.inf)
    cnt = Counter(lab); elig = np.array([i for i in range(len(pool)) if cnt[lab[i]] >= 2])
    if not len(elig): return np.nan
    order = np.argsort(-sim[elig], axis=1, kind="stable")
    rel = (lab[order] == lab[elig][:, None]).astype(float)
    rel[:, -1] = 0.0                                     # self ranks last (-inf), never relevant
    prec = np.cumsum(rel, axis=1) / np.arange(1, rel.shape[1] + 1)
    npos = rel.sum(axis=1)
    ap = (prec * rel).sum(axis=1) / np.where(npos == 0, np.nan, npos)
    return float(np.nanmean(ap))

def plate_mixing(emb, plate, comp, k=KMIX, subs=SUBS, rng=RNG):
    """batch: kNN same-plate enrichment among DIFFERENT-compound neighbours only
    (so preserved biological replicates are NOT counted as batch structure).
    1.0 = perfectly mixed; >1 = plate structure remains."""
    emb = np.ascontiguousarray(emb)
    # NOTE: duplicated wells from a block-bootstrap draw are harmless here — a duplicate of
    # well w carries w's compound, and same-compound neighbours are excluded below.
    n = len(emb); idx = rng.choice(n, min(subs, n), replace=False)
    E = emb[idx]; P = np.asarray(plate)[idx]; C = np.asarray(comp)[idx]
    nn = NearestNeighbors(n_neighbors=k+1).fit(E)
    _, ind = nn.kneighbors(E)
    fr = []
    for i in range(len(E)):
        neigh = [j for j in ind[i, 1:] if C[j] != C[i]]   # exclude same-compound replicates
        if neigh: fr.append(np.mean([P[j] == P[i] for j in neigh]))
    if not fr: return np.nan
    same = float(np.mean(fr))
    cnt = Counter(P); S = len(P); exp = sum(c*(c-1) for c in cnt.values())/(S*(S-1))
    return float(same/exp) if exp > 0 else np.nan

def build_reps(Xw, dmso, plate):
    """common 50-dim space for fair comparison."""
    Zt = StandardScaler().fit_transform(Xw)
    raw = PCA(NPC, random_state=0).fit_transform(Zt)
    # ZCA sphering on DMSO controls
    mu = dmso.mean(0); cov = np.cov((dmso-mu).T) + 1e-3*np.eye(dmso.shape[1])
    U, S, _ = np.linalg.svd(cov); W = U @ np.diag(1/np.sqrt(S)) @ U.T
    sph = PCA(NPC, random_state=0).fit_transform((Xw-mu) @ W)
    reps = {"raw": raw, "sphered": sph}
    try:
        import harmonypy
        ho = harmonypy.run_harmony(sph.astype("float32"), pd.DataFrame({"plate": plate}), ["plate"], max_iter_harmony=10)
        Zc = np.asarray(ho.Z_corr); reps["sphered+harmony"] = (Zc.T if Zc.shape[0] != len(sph) else Zc)
    except Exception as e:
        print("harmony skipped:", repr(e)[:90])
    return reps

def evaluate(reps, comp, moa, plate):
    """COMMON, PAIRED plate-block bootstrap.

    Every draw resamples the SAME set of plates and evaluates BOTH metrics on ALL
    representations, so intervals are block-level and method deltas are paired.

    Point estimates are always the OBSERVED statistics on the full data. The bootstrap
    supplies dispersion only; intervals are recentred on the observed value, because a
    block bootstrap changes the compound pool and therefore biases mean-of-draws (see the
    module docstring). `bootstrap_mean` is reported so that bias stays visible.
    """
    comp = np.asarray(comp); moa = np.asarray(moa); plate = np.asarray(plate)
    plates = np.array(sorted(set(plate))); names = list(reps)
    k = len(plates)
    point = {n: {"bio": moa_map_consensus(reps[n], comp, moa),
                 "bat": plate_mixing(reps[n], plate, comp, rng=np.random.RandomState(0))} for n in names}
    draws = {n: {"bio": [], "bat": []} for n in names}
    dref = "raw" if "raw" in names else names[0]
    delta = {n: {"bio": [], "bat": []} for n in names if n != dref}
    for b in range(B):
        samp = RNG.choice(plates, k, replace=True)                                 # one shared block draw
        mask = np.concatenate([np.where(plate == p)[0] for p in samp])             # multiplicity PRESERVED
        per = {}
        for n in names:                                                            # both metrics see the same draw
            bio = moa_map_consensus(reps[n][mask], comp[mask], moa[mask])
            bat = plate_mixing(reps[n][mask], plate[mask], comp[mask], rng=np.random.RandomState(b))
            draws[n]["bio"].append(bio); draws[n]["bat"].append(bat); per[n] = (bio, bat)
        for n in delta:
            delta[n]["bio"].append(per[n][0] - per[dref][0])
            delta[n]["bat"].append(per[n][1] - per[dref][1])

    def centred_ci(a, observed, nd=3):
        """bootstrap gives the SPREAD; the interval is anchored at the observed statistic."""
        a = np.asarray(a, float); a = a[np.isfinite(a)]
        if a.size == 0 or not np.isfinite(observed): return [None, None]
        mu = float(np.mean(a))
        lo, hi = np.percentile(a, 2.5) - mu, np.percentile(a, 97.5) - mu
        return [round(float(observed + lo), nd), round(float(observed + hi), nd)]

    out = {}
    for n in names:
        pb, pt = float(point[n]["bio"]), float(point[n]["bat"])
        out[n] = {"biology_mAP": round(pb, 3), "biology_ci": centred_ci(draws[n]["bio"], pb),
                  "batch_mixing": round(pt, 2), "batch_ci": centred_ci(draws[n]["bat"], pt, nd=2)}
        if n in delta:
            dv = np.array(delta[n]["bio"], float)
            obs = pb - float(point[dref]["bio"])                                   # OBSERVED paired difference
            fin = dv[np.isfinite(dv)]
            out[n]["delta_bio_vs_" + dref] = {
                "value": round(float(obs), 3),
                "ci": centred_ci(dv, obs),
                "bootstrap_mean": round(float(np.nanmean(dv)), 3),                 # bias kept visible
                "p_gt0": round(float(np.mean((fin - fin.mean() + obs) > 0)), 3) if fin.size else None}
        print(f"{n:16s} biology mAP={out[n]['biology_mAP']} CI{out[n]['biology_ci']}  "
              f"batch-mixing={out[n]['batch_mixing']} CI{out[n]['batch_ci']}"
              + (f"  \u0394bio_vs_{dref}={out[n]['delta_bio_vs_'+dref]['value']} CI{out[n]['delta_bio_vs_'+dref]['ci']}"
                 if n in delta else ""), flush=True)
    from math import comb
    nd = int(min(B, comb(2*k - 1, k))) if k <= 20 else B
    out["_bootstrap"] = {"draws": B, "n_blocks": int(k), "n_distinct_block_draws": nd,
                         "point_estimate": "observed on full data",
                         "interval": "bootstrap dispersion recentred on the observed value",
                         "note": ("with only %d blocks the block bootstrap has %d distinct resamples — "
                                  "intervals are EXPLORATORY, not a strict frequentist guarantee" % (k, nd))
                                 if nd < 200 else "block resolution adequate"}
    if nd < 200:
        print(f"  ! only {k} blocks -> {nd} distinct bootstrap resamples; intervals are exploratory", flush=True)
    return out

def main():
    # read the SANITIZED matrix produced by pipeline.py (P0-clean; never the old contaminated file)
    data = Path(__file__).resolve().parent.parent/"data"
    clean = sorted(data.glob("lincs_selected_v*-sanitized.parquet"))
    src = clean[-1] if clean else data/"lincs_selected.parquet"
    print("reading:", src.name)
    df = pd.read_parquet(src)
    HELPER = {"moa1", "plate", "drug"}
    feat = [c for c in df.columns if not c.startswith("Metadata") and c not in HELPER]
    assert np.isfinite(df[feat].values).all() and np.abs(df[feat].values).max() <= 15.0001, \
        "input matrix is not sanitized — run pipeline.py first"
    df["moa1"] = df["Metadata_moa"].map(lambda m: m.split("|")[0].strip() if isinstance(m, str) else None)
    trt = df[df.Metadata_broad_sample != "DMSO"].reset_index(drop=True)
    dmso = df[df.Metadata_broad_sample == "DMSO"][feat].values
    reps = build_reps(trt[feat].values, dmso, trt["Metadata_Assay_Plate_Barcode"].values)
    res = evaluate(reps, trt["Metadata_broad_sample"].values, trt["moa1"].values,
                   trt["Metadata_Assay_Plate_Barcode"].values)
    payload = {"dataset": "LINCS Cell Painting batch1 (cpg0004) · single-centre → block = plate",
               "n_wells": int(len(trt)), "n_plates": int(trt["Metadata_Assay_Plate_Barcode"].nunique()),
               "n_features": len(feat), "pca_dims": NPC, "bootstrap_draws": B,
               "metrics": {"biology": "MoA-retrieval mAP (higher=better)",
                           "batch": "plate-mixing enrichment (1.0=mixed, >1=batch structure)"},
               "results": res}
    outp = Path(__file__).resolve().parent/"mvp_results.json"
    json.dump(payload, open(outp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("wrote", outp)

if __name__ == "__main__":
    main()
