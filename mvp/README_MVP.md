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

---

## ⚠ JUMP-Target 結果的重要限制（exploratory）

在 CPJUMP1 中,**16 個 compound plate 每一盤都只含單一 cell line 與單一 timepoint**
(已由 `datasplits/cpjump1_metadata.csv` 驗證:16/16)。因此 **plate 完全巢狀於 `cell_line × timepoint`**。

這代表:

> **Harmony improved plate mixing, while pooled MoA retrieval showed no measurable loss. Because plate
> is confounded with cell line and timepoint in CPJUMP1, this result cannot distinguish removal of
> technical variation from removal of condition-specific biology.**

pooled 的 MoA 檢索之所以「沒測到損失」,正是因為它本來就跨條件彙總——它偵測不到細胞株/時間點訊號被抹掉。

### 已加入的修正(`jump_mvp.py`)

1. **自動偵測並揭露巢狀結構**,若偵測到就明講 pooled 結果不可用於因果解讀。
2. **條件訊號保留診斷** `condition_retention()`:校正前後,還分不分得出 cell line 與 timepoint。
   數值 →1 表示該條件的差異已被抹平。在注入已知條件效應的合成資料上,此指標正確地顯示
   raw 保留訊號(~1.9)而 sphering/Harmony 抹平至 ~0.95–0.99,證明診斷有效。
3. **分層評估**:在四個 `cell_line × timepoint` 分層內各自評估,此時 plate 才是純技術重複。
4. 輸出 JSON 標記 `"status": "exploratory"` 並附 `caveat` 欄位。

### 真實結果(CPJUMP1,16 盤、4,102 個處理孔)

**合併分析(混雜,僅供對照)**

| 表徵 | 生物 MoA mAP(95% CI) | 批次混合(→1 越好) | Δbio vs raw |
|---|---|---|---|
| raw | 0.154 [0.132, 0.213] | 3.17 [2.30, 3.19] | — |
| sphered(ZCA) | 0.180 [0.114, 0.240] | 1.88 [1.66, 1.87] | +0.015 [−0.057, +0.096] |
| sphered + Harmony | 0.131 [0.107, 0.219] | 1.04 [1.01, 1.11] | −0.009 [−0.066, +0.074] |

**條件訊號保留度(1.0 = 該生物差異已被抹平)**

| 表徵 | 細胞株 | 時間點 |
|---|---|---|
| raw | 1.62 | 1.28 |
| sphered(ZCA) | 1.21 | 1.10 |
| sphered + Harmony | 1.11 | **1.02** |

這一張表就是整個 MVP 最重要的產出:Harmony 把批次混合從 3.17 壓到 1.04(近乎完美),
但同一組 embedding 的時間點訊號同時從 1.28 掉到 1.02——**幾乎完全被抹平**。
「批次移除」的漂亮分數,有相當部分是刪掉真生物換來的。

**分層評估(plate 此時才是純技術重複;每層 4 盤、約 1,000 wells、260 化合物、26 個合格查詢)**

| 分層 | wells | raw mAP | Δ sphered | Δ sphered+Harmony | raw 批次混合 |
|---|---|---|---|---|---|
| A549-24h | 1,039 | 0.159 | −0.027 [−0.123, +0.062] | −0.051 [−0.131, +0.068] | 1.46 |
| A549-48h | 1,040 | 0.187 | +0.006 [−0.046, +0.035] | 0.000 [−0.047, +0.047] | 1.11 |
| U2OS-24h | 983 | 0.154 | −0.011 [−0.045, +0.043] | −0.005 [−0.072, +0.075] | 1.19 |
| U2OS-48h | 1,040 | 0.166 | +0.001 [−0.042, +0.093] | +0.002 [−0.041, +0.082] | 1.59 |

**誠實解讀**:

1. **條件內的批次效應其實不大**——raw 批次混合 1.11–1.59,遠低於合併時的 3.17。
   合併分析裡看到的「嚴重批次效應」,主要是被混雜進去的生物條件差異。
   這符合 CPJUMP1 是單一 source、單一 site pilot 的預期;真正棘手的跨 source 批次效應,此資料集看不到。
2. **八個 Δ 的 95% CI 全部跨過 0**——在此設計與樣本量下,兩種校正都測不到生物效益,也測不到明確損害。
   點估計最負的是 A549-24h 的 sphering+Harmony(−0.051)。
3. **統計功效必須誠實揭露**——每層僅 26 個合格查詢(MoA 類別需 ≥2 個化合物),合計 104 個,
   CI 因此偏寬,不足以排除中等幅度的效應。結論是「在此測不到效益」,不是「去偏無效」。

### 結論等級

| 區塊 | 定位 |
|---|---|
| LINCS（主 pipeline） | 主要量化結果,可重現、第三方可驗證 |
| LINCS（`mvp/` 雙指標） | 輔助性量化結果 |
| JUMP-Target | **探索性**去偏比較,承認條件混雜;應以分層輸出為準 |
| 影像層 | 案例式端到端展示,非統計驗證 |
