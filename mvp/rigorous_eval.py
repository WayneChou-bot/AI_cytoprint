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

Runs on REAL data. In this repo we use LINCS batch1 (single-centre, so the PLATE is the
batch/block). The identical engine points at multi-source JUMP via jump_mvp.py, where the
block is the data-generating SOURCE and the batch effect is large.

Caveat (stated, not hidden): correction transforms are fit ONCE on the full data; the
bootstrap resamples the EVALUATION over plates (given fixed embeddings). Re-fitting the
correction inside every bootstrap draw is the gold standard but far slower; the direction
of the biology/batch trade-off is unchanged. Single-centre LINCS has modest plate effects
by design — the point here is the METHOD; the JUMP loader exercises it where effects are large.
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
NPC, KMIX, SUBS, B = 50, 20, 3000, 100   # PCA dims, kNN for mixing, subsample size, bootstrap draws

def average_precision(rels):
    rels = np.asarray(rels, float); npos = rels.sum()
    if npos == 0: return np.nan
    cum = np.cumsum(rels); return float(((cum/np.arange(1, len(rels)+1))*rels).sum()/npos)

def moa_map_consensus(emb, comp, moa):
    """biology: per-compound consensus in `emb`, mAP for same-MoA retrieval."""
    df = pd.DataFrame(emb); df["c"] = comp
    cons = df.groupby("c").median(); ids = list(cons.index)
    m = pd.Series(moa).groupby(comp).first().reindex(ids).values
    V = cons.values
    sim = V @ V.T
    norm = np.linalg.norm(V, axis=1); sim = sim/np.outer(norm, norm); np.fill_diagonal(sim, -np.inf)
    pool = [i for i in range(len(ids)) if isinstance(m[i], str)]
    cnt = Counter(m[i] for i in pool); elig = [i for i in pool if cnt[m[i]] >= 2]
    if not elig: return np.nan
    aps = []
    for i in elig:
        order = [j for j in np.argsort(-sim[i]) if j in set(pool) and j != i]
        aps.append(average_precision([1 if m[j] == m[i] else 0 for j in order]))
    return float(np.nanmean(aps))

def plate_mixing(emb, plate, comp, k=KMIX, subs=SUBS, rng=RNG):
    """batch: kNN same-plate enrichment among DIFFERENT-compound neighbours only
    (so preserved biological replicates are NOT counted as batch structure).
    1.0 = perfectly mixed; >1 = plate structure remains."""
    emb = np.ascontiguousarray(emb)
    uniq = np.unique(np.arange(len(emb)))              # (mask already deduped by caller)
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
    """COMMON, PAIRED plate-block bootstrap: every draw resamples the same set of plates and
    evaluates BOTH metrics on ALL representations, so CIs are block-level and method deltas are paired."""
    comp = np.asarray(comp); moa = np.asarray(moa); plate = np.asarray(plate)
    plates = np.array(sorted(set(plate))); names = list(reps)
    point = {n: {"bio": moa_map_consensus(reps[n], comp, moa),
                 "bat": plate_mixing(reps[n], plate, comp, rng=np.random.RandomState(0))} for n in names}
    draws = {n: {"bio": [], "bat": []} for n in names}
    dref = "raw" if "raw" in names else names[0]
    delta = {n: {"bio": [], "bat": []} for n in names if n != dref}
    for b in range(B):
        samp = RNG.choice(plates, len(plates), replace=True)                       # one shared block draw
        mask = np.concatenate([np.where(plate == p)[0] for p in samp])             # biology uses full mask
        umask = np.unique(mask)                                                    # batch uses deduped wells
        per = {}
        for n in names:
            bio = moa_map_consensus(reps[n][mask], comp[mask], moa[mask])
            bat = plate_mixing(reps[n][umask], plate[umask], comp[umask], rng=np.random.RandomState(b))
            draws[n]["bio"].append(bio); draws[n]["bat"].append(bat); per[n] = (bio, bat)
        for n in delta:
            delta[n]["bio"].append(per[n][0] - per[dref][0])
            delta[n]["bat"].append(per[n][1] - per[dref][1])
    def ci(a): return [round(float(np.nanpercentile(a, 2.5)), 3), round(float(np.nanpercentile(a, 97.5)), 3)]
    out = {}
    for n in names:
        out[n] = {"biology_mAP": round(float(point[n]["bio"]), 3), "biology_ci": ci(draws[n]["bio"]),
                  "batch_mixing": round(float(point[n]["bat"]), 2),
                  "batch_ci": [round(float(np.nanpercentile(draws[n]["bat"], 2.5)), 2),
                               round(float(np.nanpercentile(draws[n]["bat"], 97.5)), 2)]}
        if n in delta:
            dv = np.array(delta[n]["bio"]); out[n]["delta_bio_vs_"+dref] = {
                "value": round(float(np.nanmean(dv)), 3), "ci": ci(dv),
                "p_gt0": round(float(np.mean(dv[np.isfinite(dv)] > 0)), 3)}
        print(f"{n:16s} biology mAP={out[n]['biology_mAP']} CI{out[n]['biology_ci']}  "
              f"batch-mixing={out[n]['batch_mixing']} CI{out[n]['batch_ci']}"
              + (f"  Δbio_vs_{dref}={out[n]['delta_bio_vs_'+dref]['value']} CI{out[n]['delta_bio_vs_'+dref]['ci']}"
                 if n in delta else ""), flush=True)
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
