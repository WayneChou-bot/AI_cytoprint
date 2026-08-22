#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JUMP rigorous MVP — runs the SAME engine (rigorous_eval.py) on real JUMP data.

Two ready targets (pick with MODE):

  MODE = "cpjump1"   JUMP-Target set from the JUMP pilot (cpg0000 / CPJUMP1).
      260 compounds, 2 cell lines (U2OS, A549), 2 timepoints, 16 compound plates.
      Real technical batch structure. Profiles live as git-LFS in
      github.com/jump-cellpainting/2021_Chandrasekaran_submitted.
      Fetch only the compound plates you need (each ~1 MB), e.g.:
        git clone https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted
        cd 2021_Chandrasekaran_submitted
        git lfs pull --include="profiles/2020_11_04_CPJUMP1/BR001169*/*_normalized.csv.gz"
      (This sandbox cannot fetch that repo's LFS; run it on your machine / Colab.)

  MODE = "jump_moa"  JUMP-MOA 90-compound set (many replicates per compound → high power).
      Metadata is in github.com/jump-cellpainting/JUMP-MOA (90 compounds + MoA).
      Profiles are on S3 (Cell Painting Gallery); read anonymously (no AWS account):
        s3://cellpainting-gallery/cpg0016-jump/... (see jump-cellpainting/datasets)
      Use s3fs anon=True + polars scan_parquet; join MoA by broad_sample.

What the engine reports (see rigorous_eval.py): for raw / sphered / sphered+Harmony,
in a common 50-D space, with a SOURCE/PLATE-aware BLOCK bootstrap 95% CI:
  · biology retention = MoA-retrieval mAP        (higher is better)
  · batch removal     = batch-mixing enrichment  (→1 is better; >1 = batch remains)
Set BATCH_COL to the block you want to correct/measure: 'plate', 'cell_line', or 'batch'.
"""
import sys, io, requests, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rigorous_eval import build_reps, evaluate   # reuse the identical engine

MODE      = "cpjump1"                 # "cpjump1" | "jump_moa"
BATCH_COL = "plate"                   # block/batch: 'plate'(技術批,預設) | 'batch'.
                                      # 注意:'cell_line' 是「生物」不是純技術 batch,對它校正會移除真實細胞株差異,
                                      #      只在明確要做「跨細胞株調和」時才用,且須把細胞株差異視為 domain 而非 noise。
CPJUMP1_REPO = "path/to/2021_Chandrasekaran_submitted"   # fallback; better: pass as arg → python jump_mvp.py <repo>

def resolve_repo():
    """Repo path from argv[1], else CPJUMP1_REPO, else auto-detect a sibling clone."""
    cands = ([sys.argv[1]] if len(sys.argv) > 1 else [CPJUMP1_REPO])
    here = Path(__file__).resolve()
    for base in [here.parent, here.parent.parent, here.parent.parent.parent, Path.cwd()]:
        cands.append(str(base/"2021_Chandrasekaran_submitted"))
    for c in cands:
        if (Path(c)/"datasplits"/"cpjump1_metadata.csv").exists():
            return Path(c)
    raise SystemExit(
        "找不到 CPJUMP1 repo。請先在 code 資料夾:\n"
        "  git clone https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted\n"
        '  cd 2021_Chandrasekaran_submitted && git lfs pull --include="profiles/2020_11_04_CPJUMP1/*/*_normalized.csv.gz" && cd ..\n'
        "再執行:python mvp\\jump_mvp.py 2021_Chandrasekaran_submitted")
REPURPOSING = ("https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/"
               "da8ae6a3bc103346095d61b4ee02f08fc85a5d98/metadata/moa/"
               "repurposing_info_external_moa_map_resolved.tsv")

def moa_by_inchikey():
    nm = pd.read_csv(io.StringIO(requests.get(REPURPOSING).text), sep="\t")
    return dict(zip(nm["InChIKey14"], nm["moa"]))

def load_cpjump1(repo):
    """Load CPJUMP1 compound plates (normalized full features) + metadata; join MoA."""
    repo = Path(repo)
    meta = pd.read_csv(repo/"datasplits"/"cpjump1_metadata.csv")
    comp = meta[meta["Metadata_experiment_type"] == "Compound"].copy()
    ik2moa = moa_by_inchikey()
    comp["moa"] = comp["Metadata_InChIKey"].str.slice(0, 14).map(ik2moa)
    frames = []
    for bc in comp["Metadata_Assay_Plate_Barcode"].unique():
        # profiles live under profiles/<batch>/<barcode>/<barcode>_normalized.csv.gz
        hits = list(repo.glob(f"profiles/*/{bc}/{bc}_normalized.csv.gz"))
        if not hits:
            continue
        if hits[0].stat().st_size < 2000:     # still a git-LFS pointer file, not real data
            raise SystemExit(
                f"{hits[0].name} 只有 {hits[0].stat().st_size} bytes = LFS 指標檔(尚未下載內容)。\n"
                '請在 repo 內執行:  git lfs pull --include="profiles/2020_11_04_CPJUMP1/*/*_normalized.csv.gz"')
        d = pd.read_csv(hits[0]); d["Metadata_Assay_Plate_Barcode"] = bc
        frames.append(d)
    if not frames:
        raise SystemExit("沒讀到任何 compound profile。確認已 `git lfs pull` 抓下 2020_11_04_CPJUMP1 的 *_normalized.csv.gz。")
    prof = pd.concat(frames, ignore_index=True)
    print(f"  loaded {len(frames)} compound plates, {len(prof)} wells")
    # CPJUMP1 profiles already carry some Metadata_* cols → drop overlaps before merge so no _x/_y suffixes
    addcols = ["Metadata_cell_line", "Metadata_timepoint", "moa", "Metadata_target", "Metadata_pert_type"]
    key = ["Metadata_Assay_Plate_Barcode", "Metadata_Well"]
    if "Metadata_Well" not in prof.columns:
        raise SystemExit(f"profile 缺 Metadata_Well 欄。實際欄位範例:{[c for c in prof.columns if c.startswith('Metadata')][:12]}")
    prof = prof.drop(columns=[c for c in addcols if c in prof.columns], errors="ignore")
    prof = prof.merge(comp[key + addcols], on=key, how="left")
    n_trt = int((prof["Metadata_pert_type"] == "trt").sum())
    print(f"  merged metadata; trt wells={n_trt}, wells with MoA={int(prof['moa'].notna().sum())}")
    if n_trt == 0:
        raise SystemExit("merge 後沒有 trt wells——多半是 profile 與 metadata 的 Well 格式對不上(如 'A01' vs 'A1')。"
                         f"\nprofile Well 範例:{sorted(prof['Metadata_Well'].dropna().unique())[:5]}"
                         f"\nmetadata Well 範例:{sorted(comp['Metadata_Well'].dropna().unique())[:5]}")
    return prof

def sanitize(df, feat):
    """Same P0-safe recipe as pipeline.sanitize_select: drop blowup(>100)/NaN features,
    winsorize to ±15, WRITE BACK, then variance+correlation selection. Returns (sel, df_with_clean_values)."""
    X = df[feat].replace([np.inf, -np.inf], np.nan)
    maxabs = X.abs().max()
    X = X.drop(columns=list(maxabs[maxabs > 100].index)).dropna(axis=1).clip(-15, 15)
    X = X[X.var()[lambda s: s > 1e-3].index]
    corr = np.corrcoef(X.values.T); cols = list(X.columns); drop = set()
    for i in range(len(cols)):
        if i in drop: continue
        for j in np.where(np.abs(corr[i]) > 0.9)[0]:
            if j > i: drop.add(j)
    sel = [cols[i] for i in range(len(cols)) if i not in drop]
    df = df.drop(columns=[c for c in feat if c not in sel]).copy()
    df[sel] = X[sel]                                   # <-- write CLEANED values back (the P0 fix)
    assert np.isfinite(df[sel].values).all() and np.abs(df[sel].values).max() <= 15.0001
    print(f"  sanitize: selected {len(sel)} clean features (max|x|={np.abs(df[sel].values).max():.1f})")
    return sel, df

def condition_retention(emb, labels, k=25, subs=3000, seed=0):
    """校正後,這個條件還分得出來嗎? kNN 同標籤富集度(1.0=已完全抹平, 越高=結構保留越多)。

    量的是 **condition-associated structure**,不是「純生物訊號」。
    因為 plate 完全巢狀於 cell_line x timepoint,這份結構同時包含
    (a) 真實的細胞株/時間點生物差異 與 (b) 該條件所屬 plate 的技術效應,
    兩者在本設計中 **無法識別(unidentifiable)**。
    所以此指標下降只能說「條件關聯結構被移除了」,
    不能反推「已證明刪掉的是真生物」——只能說生物成分有受損風險。
    """
    from sklearn.neighbors import NearestNeighbors
    from collections import Counter
    rng = np.random.RandomState(seed)
    emb = np.ascontiguousarray(emb); labels = np.asarray(labels)
    idx = rng.choice(len(emb), min(subs, len(emb)), replace=False)
    E, Lb = emb[idx], labels[idx]
    nn = NearestNeighbors(n_neighbors=min(k+1, len(E))).fit(E)
    _, ind = nn.kneighbors(E)
    same = np.mean([(Lb[ind[i, 1:]] == Lb[i]).mean() for i in range(len(E))])
    cnt = Counter(Lb); S = len(Lb)
    exp = sum(c*(c-1) for c in cnt.values())/(S*(S-1)) if S > 1 else np.nan
    return float(same/exp) if exp and exp > 0 else np.nan

def run_block(trt, feat, neg, block_col, label=""):
    """在一組樣本上跑三種表徵 + 雙指標評估。"""
    reps = build_reps(trt[feat].values,
                      neg[feat].values if len(neg) else trt[feat].values,
                      trt[block_col].values)
    res = evaluate(reps,
                   trt["Metadata_broad_sample"].values if "Metadata_broad_sample" in trt else trt.index.values,
                   trt["moa1"].values, trt[block_col].values)
    return reps, res

def main():
    if MODE == "cpjump1":
        df = load_cpjump1(resolve_repo())
    else:
        raise SystemExit("jump_moa: 尚未實作 S3 讀取。請用 MODE='cpjump1'(JUMP-Target),或自行以 s3fs(anon=True)+"
                         "polars scan_parquet 讀 cpg0016 profiles、併 jump-cellpainting/JUMP-MOA 的 90 化合物 metadata。")
    feat0 = [c for c in df.columns if not c.startswith("Metadata") and c not in ("moa",)]
    feat, df = sanitize(df, feat0)
    df["moa1"] = df["moa"].map(lambda m: m.split("|")[0].strip() if isinstance(m, str) else None)
    PLATE = "Metadata_Assay_Plate_Barcode"
    trt = df[df["Metadata_pert_type"] == "trt"].reset_index(drop=True)
    neg = df[df["Metadata_pert_type"] != "trt"]

    # ── 確認並揭露混雜結構 ───────────────────────────────────────────────
    nest = df.groupby(PLATE).agg(cl=("Metadata_cell_line", "nunique"), tp=("Metadata_timepoint", "nunique"))
    confounded = bool((nest["cl"] == 1).all() and (nest["tp"] == 1).all())
    print(f"\n  plate 巢狀於 (cell_line × timepoint): {confounded}"
          f"  ({int((nest['cl']==1).sum())}/{len(nest)} plates 單一細胞株, "
          f"{int((nest['tp']==1).sum())}/{len(nest)} 單一時間點)")
    if confounded:
        print("  → 對 plate 校正無法與『移除細胞株/時間點的生物差異』區分;以下改採分層評估。")

    out = {"mode": MODE, "n_wells": int(len(trt)), "n_plates": int(trt[PLATE].nunique()),
           "plate_confounded_with_cell_line_and_timepoint": confounded,
           "status": "exploratory",
           "caveat": ("In CPJUMP1 every compound plate contains a single cell line and a single timepoint, "
                      "so plate is perfectly nested within cell_line x timepoint. A pooled plate-block analysis "
                      "cannot separate removal of technical variation from removal of condition-specific biology. "
                      "The stratified results below evaluate within each condition, where plate is a genuine "
                      "technical replicate.")}

    # ── 1. pooled(僅供對照,已知混雜)────────────────────────────────────
    print("\n=== POOLED (plate as block — CONFOUNDED, reported for reference only) ===")
    reps_pool, res_pool = run_block(trt, feat, neg, PLATE)
    out["pooled_confounded"] = res_pool

    # ── 2. 條件訊號保留:校正前後還分得出 cell line / timepoint 嗎 ──────
    print("\n=== CONDITION-SIGNAL RETENTION (pooled embeddings; 1.0 = signal erased) ===")
    cond = {}
    for name, emb in reps_pool.items():
        cond[name] = {
            "cell_line": round(condition_retention(emb, trt["Metadata_cell_line"].values), 2),
            "timepoint": round(condition_retention(emb, trt["Metadata_timepoint"].values), 2)}
        print(f"  {name:16s} cell_line={cond[name]['cell_line']:6.2f}   timepoint={cond[name]['timepoint']:6.2f}")
    out["condition_retention"] = cond
    print("  (>1 = 該條件的生物差異仍在;接近 1 = 已被校正抹平)")

    # ── 3. 分層評估:在每個 (cell_line, timepoint) 內,plate 才是純技術批次 ──
    print("\n=== STRATIFIED (within each cell_line x timepoint; plate = pure technical block) ===")
    strata = {}
    for (cl, tp), sub in trt.groupby(["Metadata_cell_line", "Metadata_timepoint"]):
        key = f"{cl}-{tp}h"
        sub = sub.reset_index(drop=True)
        sub_neg = neg[(neg["Metadata_cell_line"] == cl) & (neg["Metadata_timepoint"] == tp)]
        cmp_moa = sub.dropna(subset=["moa1"]).drop_duplicates("Metadata_broad_sample")
        vc = cmp_moa["moa1"].value_counts()
        n_eval = int(cmp_moa["moa1"].map(vc).ge(2).sum())   # MoA 類別內至少 2 個化合物才可評估
        print(f"\n  --- {key}: wells={len(sub)}, plates={sub[PLATE].nunique()}, "
              f"compounds={sub['Metadata_broad_sample'].nunique()}, evaluable queries={n_eval} ---")
        if sub[PLATE].nunique() < 2:
            print("     (少於 2 個 plate,無法做 plate-block bootstrap — 跳過)"); continue
        try:
            _, res = run_block(sub, feat, sub_neg, PLATE)
            strata[key] = {"n_wells": int(len(sub)), "n_plates": int(sub[PLATE].nunique()),
                           "n_compounds": int(sub["Metadata_broad_sample"].nunique()),
                           "n_evaluable_queries": n_eval, "results": res}
        except Exception as e:
            print("     FAILED:", repr(e)[:160])
    out["stratified"] = strata

    import json
    dst = Path(__file__).resolve().parent/"jump_mvp_results.json"
    json.dump(out, open(dst, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\nwrote {dst.name}  (status=exploratory; 分層結果為主,pooled 僅供對照)")

if __name__ == "__main__":
    main()
