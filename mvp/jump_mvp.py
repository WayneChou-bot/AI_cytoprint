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

def main():
    if MODE == "cpjump1":
        df = load_cpjump1(resolve_repo())
    else:
        raise SystemExit("jump_moa: 尚未實作 S3 讀取。請用 MODE='cpjump1'(JUMP-Target),或自行以 s3fs(anon=True)+"
                         "polars scan_parquet 讀 cpg0016 profiles、併 jump-cellpainting/JUMP-MOA 的 90 化合物 metadata。")
    feat0 = [c for c in df.columns if not c.startswith("Metadata") and c not in ("moa",)]
    feat, df = sanitize(df, feat0)
    df["moa1"] = df["moa"].map(lambda m: m.split("|")[0].strip() if isinstance(m, str) else None)
    block = {"plate": "Metadata_Assay_Plate_Barcode", "cell_line": "Metadata_cell_line",
             "batch": "Metadata_Assay_Plate_Barcode"}[BATCH_COL]
    trt = df[df["Metadata_pert_type"] == "trt"].reset_index(drop=True)
    neg = df[df["Metadata_pert_type"] != "trt"]
    reps = build_reps(trt[feat].values, neg[feat].values if len(neg) else trt[feat].values, trt[block].values)
    res = evaluate(reps, trt["Metadata_broad_sample"].values if "Metadata_broad_sample" in trt else trt.index.values,
                   trt["moa1"].values, trt[block].values)
    import json
    json.dump({"mode": MODE, "block": BATCH_COL, "n_wells": int(len(trt)), "results": res},
              open(Path(__file__).resolve().parent/"jump_mvp_results.json", "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("wrote jump_mvp_results.json")

if __name__ == "__main__":
    main()
