#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
image_pipeline.py — 真正的「影像 → 形態指紋」管線(這是作品的 CV 部分)。

輸入:CPJUMP1 repo 內的真實 Cell Painting 影像(example_images/,一般 git clone 就有,不需 LFS)
輸出:
  images/out/webimages.json  ── 每個化合物的通道圖(base64)、分割疊圖、per-cell 特徵摘要
  images/out/*.png           ── 中繼圖(可單獨查看)

流程(對應 CellProfiler 的簡化重製版):
  1. 讀 5 個螢光通道 TIFF(1080×1080, uint16)
  2. 以 DNA 通道分割細胞核(高斯 → Otsu → 去雜訊 → 分水嶺切開相黏核)
  3. 由核外擴得到細胞質環(近似細胞區域)
  4. 每顆細胞量化:形狀(面積/周長/離心率/緊緻度) + 各通道亮度(核內/質內) + 紋理(Laplacian/局部變異)
  5. 與 DMSO 對照比較,輸出 z-score 指紋

通道對應(重要,見 README 說明):
  ch5 = DNA        ← 視覺與分割確認(圓形核、背景乾淨)
  ch1 = Mito       ← 核內/質比最低、核仁對比最低(MitoTracker 不染核仁)
  ch2/3/4 = AGP / RNA / ER  ← 由上述兩端錨定,依「波長由長到短」推定;非官方文件確認
  ch6-8 = Brightfield(3 個 z 平面)
  權威對應在 repo 的 load_data csv(git-lfs),可自行驗證:
    git lfs pull --include="load_data_csv/2020_11_04_CPJUMP1/BR00117010/load_data.csv.gz"
    python -c "import pandas as pd;d=pd.read_csv('...load_data.csv.gz');print([c for c in d.columns if 'FileName' in c]);print(d.iloc[0])"
"""
import json, base64, io, sys
from pathlib import Path
import numpy as np, tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, morphology, measure, segmentation, exposure, feature

CH = {1: "Mito", 2: "AGP", 3: "RNA", 4: "ER", 5: "DNA"}          # 見檔頭說明
COLOR = {"DNA": (90, 134, 255), "RNA": (90, 224, 214), "ER": (90, 208, 122),
         "AGP": (255, 154, 70), "Mito": (255, 107, 208)}
ORDER = ["DNA", "RNA", "ER", "AGP", "Mito"]
OUT = Path(__file__).resolve().parent/"out"; OUT.mkdir(exist_ok=True)
THUMB, CHTHUMB = 460, 300

def find_repo():
    cands = ([sys.argv[1]] if len(sys.argv) > 1 else [])
    here = Path(__file__).resolve()
    for b in [here.parent, here.parent.parent, here.parent.parent.parent, Path.cwd()]:
        cands.append(str(b/"2021_Chandrasekaran_submitted"))
    for c in cands:
        if (Path(c)/"example_images").is_dir(): return Path(c)
    raise SystemExit("找不到 2021_Chandrasekaran_submitted/example_images。\n"
                     "請先 git clone https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted\n"
                     "（TIFF 影像一般 clone 即有，不需 git lfs pull）")

def load_field(folder):
    """回傳 {stain: 2D float array}，只取 5 個螢光通道。"""
    out = {}
    for f in sorted(Path(folder).glob("*.tiff")):
        try: ch = int(f.name.split("-ch")[1].split("sk")[0])
        except Exception: continue
        if ch in CH: out.setdefault(CH[ch], tifffile.imread(f).astype(np.float32))
    return out

def stretch(a, lo=1.0, hi=99.7):
    p1, p2 = np.percentile(a, lo), np.percentile(a, hi)
    return np.clip((a - p1)/max(p2 - p1, 1e-6), 0, 1)

def png_b64(arr_uint8, size, mode="L", fmt="JPEG", q=82):
    im = Image.fromarray(arr_uint8, mode=mode).resize((size, size), Image.LANCZOS)
    if fmt == "JPEG" and mode == "L": im = im.convert("L")
    buf = io.BytesIO(); im.save(buf, format=fmt, quality=q, optimize=True)
    return f"data:image/{'jpeg' if fmt=='JPEG' else 'png'};base64," + base64.b64encode(buf.getvalue()).decode()

GAIN = {"DNA": 0.95, "RNA": 0.42, "ER": 0.40, "AGP": 0.50, "Mito": 0.55}

def composite(chs):
    """五通道疊成彩色合成圖。用較低增益 + gamma,避免五通道相加後過曝成白色。"""
    h, w = next(iter(chs.values())).shape
    rgb = np.zeros((h, w, 3), np.float32)
    for st in ORDER:
        if st not in chs: continue
        g = stretch(chs[st], 3, 99.2) ** 0.85          # gamma:拉起暗部又不吃掉亮部
        rgb += g[..., None]*(np.array(COLOR[st], np.float32)/255.0)[None, None, :]*GAIN[st]
    rgb = rgb/max(np.percentile(rgb, 99.5), 1e-6)      # 依 99.5 百分位歸一,保留少量高光
    return (np.clip(rgb, 0, 1)*255).astype(np.uint8)

def seg_overlay(dna, lab):
    """分割檢視:DNA 通道(灰藍) + 黃色分割邊界,對比清楚。"""
    g = stretch(dna, 2, 99.5) ** 0.9
    rgb = np.stack([g*0.55, g*0.72, g*1.0], -1)        # 冷色調 DNA
    b = segmentation.find_boundaries(lab, mode="outer")
    b = morphology.dilation(b, morphology.square(2))   # 加粗,縮圖後才看得見
    rgb[b] = [1.0, 0.85, 0.15]                          # 黃色邊界
    return (np.clip(rgb, 0, 1)*255).astype(np.uint8)

MIN_NUC_AREA = 250

def drop_small(labels, min_area):
    """自寫版本,避免 skimage 各版本 remove_small_objects 參數語意變動。"""
    labels = np.asarray(labels)
    if labels.dtype == bool: labels = measure.label(labels)
    counts = np.bincount(labels.ravel())
    small = np.where(counts < min_area)[0]
    out = labels.copy(); out[np.isin(labels, small)] = 0
    return out

def segment_nuclei(dna):
    """高斯平滑 → Otsu → 去小物件/補洞 → 距離變換分水嶺切開相黏核。"""
    sm = filters.gaussian(dna, 2)
    mask = ndi.binary_fill_holes(sm > filters.threshold_otsu(sm))
    mask = drop_small(mask, MIN_NUC_AREA) > 0
    dist = ndi.distance_transform_edt(mask)
    peaks = feature.peak_local_max(dist, min_distance=14, labels=mask)
    mk = np.zeros(dist.shape, int); mk[tuple(peaks.T)] = np.arange(1, len(peaks)+1)
    lab = segmentation.watershed(-dist, mk, mask=mask)
    lab = drop_small(lab, MIN_NUC_AREA)
    return measure.label(segmentation.clear_border(lab))   # 去掉貼邊被切斷的核

def cyto_territory(lab, radius=18):
    """每顆細胞專屬的細胞質環。

    先前版本用 `dilation(所有核) & ~所有核` 取一個**全域**環帶,在密集視野中一顆細胞的
    「細胞質」會混進鄰近細胞的訊號。改為:把環帶內每個像素指派給**最近的細胞核**,
    讓每顆細胞擁有互斥的細胞質territory(等同對距離變換做最近鄰分割)。
    """
    dist, inds = ndi.distance_transform_edt(lab == 0, return_indices=True)
    nearest = lab[inds[0], inds[1]]
    ring = (dist > 0) & (dist <= radius)
    return np.where(ring, nearest, 0)

def features_per_cell(chs, lab):
    """每顆細胞一列特徵:形狀 + 各通道核內/質內亮度 + 紋理。"""
    dna = chs["DNA"]
    cyto = cyto_territory(lab)
    lap = {st: ndi.laplace(filters.gaussian(stretch(a), 1)) for st, a in chs.items()}
    rows, cols = [], None
    for p in measure.regionprops(lab, intensity_image=dna):
        if p.area < 250: continue
        cy, cx = p.centroid
        sl = (slice(max(0, int(cy)-45), int(cy)+45), slice(max(0, int(cx)-45), int(cx)+45))
        nm = lab[sl] == p.label
        cm = cyto[sl] == p.label            # 互斥領域:不含鄰近細胞
        if nm.sum() < 50 or cm.sum() < 50: continue
        f = {"AreaShape_Area": p.area, "AreaShape_Perimeter": p.perimeter,
             "AreaShape_Eccentricity": p.eccentricity, "AreaShape_Solidity": p.solidity,
             "AreaShape_Extent": p.extent,
             "AreaShape_FormFactor": 4*np.pi*p.area/max(p.perimeter**2, 1e-6),
             "AreaShape_MajorAxis": p.axis_major_length, "AreaShape_MinorAxis": p.axis_minor_length,
             "AreaShape_AspectRatio": p.axis_major_length/max(p.axis_minor_length, 1e-6)}
        for st in ORDER:
            if st not in chs: continue
            a = chs[st][sl]; bg = np.percentile(chs[st], 3)
            nv, cv = a[nm]-bg, a[cm]-bg
            f[f"Intensity_Nuc_Mean_{st}"] = nv.mean()
            f[f"Intensity_Nuc_Std_{st}"] = nv.std()
            f[f"Intensity_Nuc_P95_{st}"] = np.percentile(nv, 95)
            f[f"Intensity_Cyto_Mean_{st}"] = cv.mean()
            f[f"Intensity_Cyto_Std_{st}"] = cv.std()
            f[f"Intensity_NucCytoRatio_{st}"] = nv.mean()/max(cv.mean(), 1e-6)
            f[f"Intensity_IntegratedNuc_{st}"] = nv.sum()
            l = lap[st][sl]
            f[f"Texture_LapRMS_Nuc_{st}"] = float(np.sqrt((l[nm]**2).mean()))
            f[f"Texture_LapRMS_Cyto_{st}"] = float(np.sqrt((l[cm]**2).mean()))
            f[f"Texture_CV_Cyto_{st}"] = cv.std()/max(abs(cv.mean()), 1e-6)
        rows.append(f); cols = list(f)
    return rows, (cols or [])

def main():
    repo = find_repo(); root = repo/"example_images"
    comps = sorted([d.name for d in root.iterdir() if d.is_dir()])
    print(f"repo: {repo}\ncompounds: {comps}")
    per_cell, web = {}, {}
    for name in comps:
        chs = load_field(root/name)
        if len(chs) < 5: print(f"  skip {name}: only {len(chs)} channels"); continue
        lab = segment_nuclei(chs["DNA"])
        rows, cols = features_per_cell(chs, lab)
        per_cell[name] = rows
        # web assets
        imgs = {"composite": png_b64(composite(chs), THUMB, "RGB", "JPEG", 84),
                "segmented": png_b64(seg_overlay(chs["DNA"], lab), THUMB, "RGB", "JPEG", 86)}
        for st in ORDER:
            imgs[st] = png_b64((stretch(chs[st])*255).astype(np.uint8), CHTHUMB, "L", "JPEG", 78)
        web[name] = {"images": imgs, "n_cells": int(lab.max()), "n_measured": len(rows)}
        print(f"  {name:14s} nuclei={lab.max():4d}  cells_measured={len(rows):4d}")
    if "DMSO" not in per_cell: raise SystemExit("找不到 DMSO 對照,無法計算 z-score 指紋")

    # DMSO 為基準的 z-score 指紋(以 DMSO 的 per-cell 分布做標準化)
    cols = list(per_cell["DMSO"][0])
    D = np.array([[r[c] for c in cols] for r in per_cell["DMSO"]], float)
    mu, sd = D.mean(0), D.std(0) + 1e-9
    Zall = {}
    for name, rows in per_cell.items():
        X = np.array([[r[c] for c in cols] for r in rows], float)
        Z = (X - mu)/sd; Zall[name] = Z
        z = Z.mean(0)
        web[name]["fingerprint"] = {c: round(float(v), 2) for c, v in zip(cols, z)}
        web[name]["dist_from_dmso"] = round(float(np.linalg.norm(z)/np.sqrt(len(cols))), 3)
        ai = cols.index("AreaShape_Area")
        web[name]["nuc_area"] = [int(v) for v in X[:, ai][:250]]           # 原始尺度,畫分布圖
        web[name]["nuc_area_median"] = int(np.median(X[:, ai]))
    # 每個化合物的 top 驅動特徵(對 DMSO 的 Cohen's d)
    for name in web:
        if name == "DMSO": web[name]["top_features"] = []; continue
        A, B = Zall[name], Zall["DMSO"]
        s = np.sqrt(((len(A)-1)*A.var(0) + (len(B)-1)*B.var(0))/max(len(A)+len(B)-2, 1)) + 1e-9
        d = (A.mean(0) - B.mean(0))/s
        idx = np.argsort(-np.abs(d))[:6]
        web[name]["top_features"] = [{"f": cols[i], "d": round(float(d[i]), 2)} for i in idx]
    # 以「我的影像指紋」互相檢索:最相似的化合物
    # DMSO 的指紋依定義是零向量(自己對自己標準化),cosine 對零向量無意義 → 排除在排序之外
    names = [n for n in web if n != "DMSO"]
    V = np.array([[web[n]["fingerprint"][c] for c in cols] for n in names])
    Vn = V/(np.linalg.norm(V, axis=1, keepdims=True) + 1e-9)
    S = Vn @ Vn.T; np.fill_diagonal(S, -2)
    for i, n in enumerate(names):
        order = np.argsort(-S[i])[:3]
        web[n]["neighbors"] = [{"name": names[j], "sim": round(float(S[i][j]), 2)} for j in order]
    if "DMSO" in web:
        web["DMSO"]["neighbors"] = []       # 對照組本身:無最近鄰(零向量)
    payload = {"channel_map": {f"ch{k}": v for k, v in CH.items()},
               "channel_note": "ch5=DNA 與 ch1=Mito 為實證確認;ch2/3/4 依波長順序推定(見程式檔頭)",
               "channel_note_en": "ch5=DNA and ch1=Mito determined empirically; ch2/3/4 inferred from acquisition wavelength order (see module docstring)",
               "source": "JUMP pilot CPJUMP1 (cpg0000) example_images · plates BR00117010-13 · U2OS · 1 field/compound",
               "feature_names": cols, "n_features": len(cols), "order": ORDER, "compounds": web}
    (OUT/"webimages.json").write_text(json.dumps(payload), encoding="utf-8")
    # 另存完整 per-cell 特徵表(可重現/可驗證)
    import csv
    with open(OUT/"per_cell_features.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["compound"]+cols)
        for name, rows in per_cell.items():
            for r in rows: w.writerow([name]+[r[c] for c in cols])
    print(f"\nwrote {OUT/'webimages.json'} ({(OUT/'webimages.json').stat().st_size/1e6:.1f} MB, {len(cols)} features)")
    print(f"wrote {OUT/'per_cell_features.csv'} ({sum(len(v) for v in per_cell.values())} cells)")
    rank = sorted(((v["dist_from_dmso"], k) for k, v in web.items()), reverse=True)
    print("\n離 DMSO 最遠(表型最強):")
    for d, k in rank: print(f"  {k:14s} dist={d:5.3f}  核面積中位數={web[k]['nuc_area_median']:5d}  n={web[k]['n_measured']}")

if __name__ == "__main__":
    main()
