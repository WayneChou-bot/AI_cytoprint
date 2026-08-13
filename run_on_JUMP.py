#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_on_JUMP.py — 把 MorphoProfile 的管線指向「真正的 JUMP cpg0016」。

為什麼要另外一支：本沙盒環境的網路擋掉了 AWS S3，所以雲端這邊改用 LINCS pilot 跑。
JUMP 的 profiles 在 S3 開放桶上（免帳號），在**你自己的電腦或 Colab/Kaggle** 上可直接讀。

── 資料量與記憶體（已修正先前過度樂觀的說法）─────────────────────────────
• 你**不要**碰原始影像：cpg0016 原圖約 127 TB（2,380 盤、~110 萬 wells；2026-08 由
  datasets/stats/cpg0016_source_images_size.csv 加總）。
• 你**要用**的是各 subset 的分析就緒 profiles parquet（各 subset 前處理不同）：
    - CRISPR：sphering+Harmony+PCA correction   - ORF：sphering+Harmony
    - COMPOUND：featureselect+Harmony（檔名沒有 sphering）
  用本檔 `print_sizes()`（s3fs 只讀 metadata、不下載內容）或下列指令**先確認實際大小**：
    aws s3 ls --no-sign-request --recursive --summarize --human-readable \
      s3://cellpainting-gallery/cpg0016-jump-assembled/source_all/workspace/profiles_assembled/COMPOUND/
• **記憶體誠實話**：
    - 本檔 load_profiles 用 polars collect() 一次讀入，N_SAMPLE 只在讀完後抽樣，
      **省不了下載量與峰值記憶體**；要真的省請改 row-group 過濾/分塊讀取。
    - **不要對全量 ~11.6 萬化合物建完整 cosine 矩陣**：11.6 萬² ≈ 134 億個值，
      光 float64 就 >100 GB。全量檢索**必須**用 FAISS 或近似最近鄰／分塊 top-k。
    - 因此筆電適合「讀共識/抽樣 profile 做小規模分析」，或先跑 JUMP-MOA 90 化合物 /
      JUMP-Target control set 的嚴謹 MVP；不要在筆電上跑全量 all-pairs 檢索。
    - 只有要做影像級深度學習（自監督 embedding）才需要 Colab/Kaggle 的 GPU。
──────────────────────────────────────────────────────────────────────

需求套件：pip install pandas pyarrow s3fs polars scikit-learn umap-learn requests
（Colab：第一格 pip 安裝即可，S3 匿名讀取無須任何金鑰）
"""
import json, requests, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

# ── 設定 ─────────────────────────────────────────────
SUBSET   = "compound"     # 'compound'（做 MoA 檢索用這個）| 'crispr' | 'orf'（基因，較小，快測）
N_SAMPLE = 30000          # 安全預設：先抽樣跑通。設 None=全量前，請先看檔頭「記憶體誠實話」——
                          # 全量會 collect() 整份、且不可對 ~11.6 萬化合物建 dense cosine（>100GB）。
MAX_DENSE = 4000          # 化合物數超過此值就改用 top-k 最近鄰，不建整個相似度矩陣
INDEX = "https://raw.githubusercontent.com/jump-cellpainting/datasets/v0.11.0/manifests/profile_index.json"
REPURPOSING = ("https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/"
               "da8ae6a3bc103346095d61b4ee02f08fc85a5d98/metadata/moa/"
               "repurposing_info_external_moa_map_resolved.tsv")  # InChIKey14 -> MoA

def profile_url(subset):
    idx = requests.get(INDEX, timeout=60).json()
    return next(e["url"] for e in idx if e["subset"] == subset)

def print_sizes():
    """列出 JUMP 各 profiles 檔案大小（用 s3fs 讀 metadata，不下載內容）。"""
    import s3fs, pyarrow.parquet as pq
    fs = s3fs.S3FileSystem(anon=True)
    idx = requests.get(INDEX, timeout=60).json()
    for e in idx:
        key = e["url"].split(".amazonaws.com/")[1]
        try:
            info = fs.info(key); md = pq.ParquetFile(key, filesystem=fs).metadata
            print(f"{e['subset']:20s} {info['size']/1e9:6.2f} GB  rows={md.num_rows:>8}  cols={md.num_columns}")
        except Exception as ex:
            print(e["subset"], "ERR", repr(ex)[:80])

def load_profiles(subset, n_sample=None):
    """讀 JUMP well-level profiles。用 polars 掃描，避免一次載入全部。"""
    import polars as pl
    url = profile_url(subset)
    print("loading:", url)
    lf = pl.scan_parquet(url)
    df = lf.collect().to_pandas()
    if n_sample and len(df) > n_sample:
        df = df.sample(n=n_sample, random_state=0).reset_index(drop=True)
    meta = [c for c in df.columns if c.startswith("Metadata")]
    feat = [c for c in df.columns if not c.startswith("Metadata")]
    print(f"  {subset}: {df.shape[0]} wells × {len(feat)} features")
    return df, meta, feat

def attach_moa(df):
    """為化合物 well 併上 MoA：JUMP metadata(InChIKey) → Drug Repurposing Hub(InChIKey14→MoA)。"""
    comp = pd.read_csv("https://raw.githubusercontent.com/jump-cellpainting/datasets/"
                       "main/metadata/compound.csv.gz")  # Metadata_JCP2022, Metadata_InChIKey
    comp["ik14"] = comp["Metadata_InChIKey"].str.slice(0, 14)
    rep = pd.read_csv(REPURPOSING, sep="\t")
    ik2moa = dict(zip(rep["InChIKey14"], rep["moa"]))
    ik2name = dict(zip(rep["InChIKey14"], rep["pert_iname"]))
    comp["moa"] = comp["ik14"].map(ik2moa)
    comp["name"] = comp["ik14"].map(ik2name)
    j = df.merge(comp[["Metadata_JCP2022", "moa", "name"]], on="Metadata_JCP2022", how="left")
    print("  compounds with MoA annotation:", j["moa"].notna().sum(), "/", len(j))
    return j

def analyze(df, feat):
    """consensus per compound → 標準化 → 最近鄰 MoA 檢索。回報 Top-1/Top-5 與**排除自身**的隨機基準。
    註：此函式只做最近鄰檢索(不算 copairs mAP、也不寫 webdata.json)；要產完整網頁 JSON 請沿用 pipeline.py。
    為安全,化合物數 > MAX_DENSE 時改用 top-k 最近鄰,不建 dense 相似度矩陣。"""
    from sklearn.preprocessing import StandardScaler
    from sklearn.neighbors import NearestNeighbors
    from collections import Counter
    key = "name" if "name" in df.columns and df["name"].notna().any() else "Metadata_JCP2022"
    d = df.dropna(subset=["moa"]).copy()
    d["moa1"] = d["moa"].str.split("|").str[0].str.strip()
    cons = d.groupby(key)[feat].median()
    moa = d.groupby(key)["moa1"].first().reindex(cons.index).values
    ids = list(cons.index); P = len(ids)
    Z = StandardScaler().fit_transform(cons.values)
    Zn = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-9)   # cosine via normalised dot / kNN
    cnt = Counter(m for m in moa if isinstance(m, str))
    ev = [i for i in range(P) if isinstance(moa[i], str) and cnt[moa[i]] >= 2]
    K = 6
    if P > MAX_DENSE:
        print(f"  {P} compounds > MAX_DENSE({MAX_DENSE}) → top-k NearestNeighbors (no dense matrix)")
        nn = NearestNeighbors(n_neighbors=min(K+1, P), metric="cosine").fit(Z)
        _, ind = nn.kneighbors(Z)
        order_of = lambda i: [j for j in ind[i, 1:] if isinstance(moa[j], str)]
    else:
        sim = Zn @ Zn.T; np.fill_diagonal(sim, -2)
        order_of = lambda i: [j for j in np.argsort(-sim[i]) if isinstance(moa[j], str)][:5]
    t1 = np.mean([moa[order_of(i)[0]] == moa[i] for i in ev if order_of(i)])
    t5 = np.mean([any(moa[j] == moa[i] for j in order_of(i)[:5]) for i in ev])
    # random baseline EXCLUDING self, matched to eligible-query class sizes (not Σp²)
    rnd = float(np.mean([(cnt[moa[i]] - 1) / (P - 1) for i in ev]))
    print(f"  JUMP MoA retrieval: top1={t1:.3f} ({t1/rnd:.1f}x self-excluded), top5={t5:.3f}, n_eval={len(ev)}")
    return dict(top1=round(float(t1), 3), top5=round(float(t5), 3), fold=round(t1/rnd, 1),
                random=round(rnd, 4), n_eval=len(ev))

if __name__ == "__main__":
    # 1) 先看檔案大小（可選）
    # print_sizes()
    # 2) 載入 → 併 MoA（compound）→ 分析
    df, meta, feat = load_profiles(SUBSET, N_SAMPLE)
    if SUBSET == "compound":
        df = attach_moa(df)
        analyze(df, feat)
    print("done. 要產出 UMAP/檢索/命中 的完整 webdata.json，請沿用 pipeline.py 的輸出區塊(相同 schema)。")
