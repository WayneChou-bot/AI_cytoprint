# 嚴謹 MVP：批次去偏的雙指標評估（source/plate-aware block bootstrap）

回應朋友建議的「小而嚴謹」下一步。這個 MVP 不是再加功能，而是把**評估方法**補到能扛技術追問。

## 為什麼不是直接跑 JUMP-MOA？

JUMP / CPJUMP1 的 profiles 在 AWS S3 與 GitHub-LFS 上，**這個雲端沙盒的網路把兩者都擋掉了**
（S3 403、該 repo 的 LFS 未授權、media CDN 對它 404）。所以：

- **評估引擎在這裡用真正抓得到的資料（LINCS batch1）跑出真數字＋CI**；
- **同一支引擎**（`jump_mvp.py`）已寫好、metadata 已驗證，指向 **JUMP-Target（CPJUMP1）** 或 **JUMP-MOA 90**，
  在**你的電腦 / Colab**（S3 或 `git lfs pull` 可達）上即可跑，換的只是 loader 與 block 欄位。

已驗證可用的 JUMP metadata（不需 profiles 也能確認）：JUMP-MOA 90 化合物＋MoA；CPJUMP1 260 化合物、
2 細胞株（U2OS/A549）、2 時間點、16 compound plates，附 target/cell_line/plate。

## 三個嚴謹要點（`rigorous_eval.py`）

1. **同一樣本、同一 50 維空間**比較 raw / sphered(ZCA) / sphered+Harmony —— 差異來自方法，不是輸入。
2. **雙指標**（必須同時看，否則會自欺）：
   - 生物保存 = **MoA 檢索 mAP**（越高越好）
   - 批次移除 = **plate-mixing enrichment**（→1 越好，>1 表仍有批次結構）
   - Harmony 會「依定義」優化批次混合，所以真正該問的是它讓生物訊號付出多少代價。
3. **plate-aware BLOCK bootstrap** 95% CI —— 生物軸重抽「整個 plate（block）」而非單一化合物，
   CI 才反映批次層級的變異。（批次軸的 CI 用 evaluation 子抽樣，避免「整盤重複」對 enrichment 比值的偏差。）

## 真實結果（LINCS batch1，block = plate）

> 本表為**修掉極端值污染（P0）後**、用 clean matrix + paired plate-block bootstrap 重算的結果。

| 表徵 | 生物訊號 MoA mAP（95% CI） | 批次混合（95% CI，→1 越好） | Δbio vs raw（paired block CI） |
|---|---|---|---|
| raw（標準化） | 0.087 [0.081, 0.10] | 2.42 [2.21, 2.65] | — |
| sphered（ZCA） | 0.087 [0.075, 0.089] | 3.0 [2.83, 3.28] | −0.008 [−0.016, +0.001] |
| sphered+Harmony（依 plate） | **0.095** [0.083, 0.096] | **1.63** [1.55, 1.77] | −0.00 [−0.009, +0.010] |

**誠實解讀**：在這份**單一中心** pilot 上，三種方法的**生物訊號幾乎相同**（兩兩 paired block Δ 的 95% CI **都跨 0**，
無顯著差異）。批次軸上 **Harmony 明顯改善混合**（2.42× → 1.63×，往 1 靠），且**沒有付出可測的生物代價**；
單獨球化則略微變差。

**兩個關鍵教訓**：
1. **資料品質 gate 一切**：先前這張表曾出現「Harmony 兩軸皆輸（mixing 43×、biology 0.027）」——那是**極端值污染的假象**；
   清理後結論反轉。沒有乾淨矩陣，任何評估都不可信。
2. **雙指標＋block CI 缺一不可**：只看生物軸會漏掉「Harmony 改善了 batch 混合」；只看批次軸會漏掉「有沒有犧牲生物」。
   本結果為單一中心、批次效應小,差異多不顯著,**不能一般化**；到多來源、批次效應大的 **JUMP**（以 **source** 為 block）,
   同一引擎才是校正真正該發揮的場景。

## 怎麼在 JUMP 上跑（你的電腦 / Colab）

```bash
# 方式 A：JUMP-Target（CPJUMP1，git-lfs）
git clone https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted
cd 2021_Chandrasekaran_submitted
git lfs pull --include="profiles/2020_11_04_CPJUMP1/BR001169*/*_normalized.csv.gz"   # 只抓 compound plates
# 編輯 jump_mvp.py：MODE="cpjump1"、CPJUMP1_REPO=該 clone 路徑、BATCH_COL="cell_line" 或 "plate"
python jump_mvp.py          # 產 jump_mvp_results.json
python make_tradeoff.py     # 產 tradeoff 圖（改讀 jump_mvp_results.json 即可）

# 方式 B：JUMP-MOA 90（S3 匿名，需補上 cpg0016 profiles 讀取；metadata 見 jump-cellpainting/JUMP-MOA）
```

## 檔案

- `rigorous_eval.py` — 評估引擎（雙指標＋block bootstrap），跑 LINCS 真數字
- `jump_mvp.py` — 同引擎的 JUMP-Target/JUMP-MOA loader（你的電腦/Colab 執行）
- `make_tradeoff.py` — 由結果 JSON 畫雙指標權衡圖
- `mvp_results.json` / `tradeoff.html` / `tradeoff.png` — 產物
