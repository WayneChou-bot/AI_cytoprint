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
| raw（標準化） | 0.087 [0.080, 0.096] | 2.42 [2.10, 2.73] | — |
| sphered（ZCA） | 0.087 [0.079, 0.096] | 2.87 [2.60, 3.16] | −0.000 [−0.009, +0.008] |
| sphered+Harmony（依 plate） | **0.095** [0.089, 0.102] | **1.61** [1.46, 1.77] | **+0.008** [−0.001, +0.017] |

> 上表為 **B=1000、點估計採觀測值、區間以觀測值重新錨定**的版本（見下方「估計量的修正」）。
> 舊版把 bootstrap 抽樣平均當成 Δ 回報，Harmony 那格是 −0.00——與觀測值 **+0.008 正負號相反**。

**誠實解讀**：在這份**單一中心** pilot 上，三種方法的**生物訊號差異都不顯著**（paired block Δ 的 95% CI 都跨 0）。
批次軸上 **Harmony 明顯改善混合**（2.42× → 1.63×，往 1 靠），生物軸則是 **+0.008 [−0.001, +0.017]**
——方向為正、且 CI 幾乎不含 0（P(Δ>0)=0.96），但仍**未達顯著**。正確說法是「在此樣本量下測不到生物代價，
且有微弱的正向趨勢」，不是「證明沒有代價」，更不是「證明有幫助」。單獨球化則在批次軸上明顯變差（2.42× → 2.87×）。

**估計量的修正（重要）**：mAP 對「池子裡有哪些化合物」是非線性的，而 block bootstrap 每一次抽樣都在換池子，
因此**抽樣平均是有偏估計**。舊版 `value` 取的是抽樣平均，導致 Harmony 的 Δ 被報成 −0.00，
而觀測差值其實是 +0.008（B=1000 時抽樣平均為 +0.001，仍低估八倍）。現在的做法是：`value` 一律為**全資料觀測值**，bootstrap 只用來估計離散程度，
區間再以觀測值重新錨定，並額外輸出 `bootstrap_mean` 讓落差看得見。
另外每份 JSON 都會寫 `n_distinct_block_draws`——block bootstrap 對 k 個 plate 只有 `C(2k-1, k)` 種相異重抽，
k=4 時僅 35 種，低於 200 會直接印警告。

**三個關鍵教訓**：
1. **資料品質 gate 一切**：先前這張表曾出現「Harmony 兩軸皆輸（mixing 43×、biology 0.027）」——那是**極端值污染的假象**；
   清理後結論反轉。沒有乾淨矩陣，任何評估都不可信。
2. **估計量本身也會說謊**：同一份資料、同一個 bootstrap，只是把「抽樣平均」換成「觀測值」，
   Harmony 的 Δ 就從 −0.00 變成 +0.008（CI 也從跨 0 很深變成幾乎不含 0）。報 CI 之前，先確認點估計是什麼東西。
3. **雙指標＋block CI 缺一不可**：只看生物軸會漏掉「Harmony 改善了 batch 混合」；只看批次軸會漏掉「有沒有犧牲生物」。
   本結果為單一中心、批次效應小，差異多不顯著，**不能一般化**；同一引擎跑的 CPJUMP1 pilot **同樣是單一 source**，
   且 plate 與生物條件混雜（見下）。真正的多來源校正（以 **source** 為 block）仍是未來工作，本 repo 尚未示範。

## 怎麼在 JUMP 上跑（你的電腦 / Colab）

```bash
# 方式 A：JUMP-Target（CPJUMP1，git-lfs）
git clone https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted
cd 2021_Chandrasekaran_submitted
git lfs pull --include="profiles/2020_11_04_CPJUMP1/BR001169*/*_normalized.csv.gz"   # 只抓 compound plates
# 不需要編輯任何常數:直接把 clone 路徑當參數傳進去
python jump_mvp.py 2021_Chandrasekaran_submitted   # 產 jump_mvp_results.json
#   block 一律是 plate(唯一的技術單位);cell_line / timepoint 是「生物條件」,
#   會被用來分層,絕不可當成 batch 來校正。
python make_tradeoff.py     # 產 tradeoff 圖（改讀 jump_mvp_results.json 即可）

# 方式 B：JUMP-MOA 90（S3 匿名，需補上 cpg0016 profiles 讀取；metadata 見 jump-cellpainting/JUMP-MOA）
```

## 檔案

- `rigorous_eval.py` — 評估引擎（雙指標＋block bootstrap），跑 LINCS 真數字
- `jump_mvp.py` — 同引擎的 JUMP-Target/JUMP-MOA loader（你的電腦/Colab 執行）
- `test_condition_retention.py` — 條件診斷的合成資料驗證測試(`python mvp/test_condition_retention.py`,不需任何下載資料)
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
2. **條件關聯結構診斷** `condition_retention()`:校正前後,還分不分得出 cell line 與 timepoint。
   數值 →1 表示該條件已無法分辨。它量的是 **condition-associated structure**,不是「純生物訊號」
   ——因為 plate 與條件完全混雜,這份結構同時含真實生物差異與該條件所屬 plate 的技術效應,兩者**無法識別**。
   驗證方式已寫成可重跑的測試:`python mvp/test_condition_retention.py`(不需下載任何資料)。
   它用與 CPJUMP1 相同的巢狀結構造合成資料、注入已知條件效應,並斷言此指標在
   raw 上 >1.3、在「完美校正」後 <1.1,且在無效應對照與打亂標籤下皆 ~1.0。實測輸出:
   `raw 1.96 / corrected 1.00 / no-effect 1.00 / shuffled 1.00 → PASS`。
3. **分層評估**:在四個 `cell_line × timepoint` 分層內各自評估,此時 plate 才是純技術重複。
4. 輸出 JSON 標記 `"status": "exploratory"` 並附 `caveat` 欄位。

### 真實結果(CPJUMP1,16 盤、4,102 個處理孔)

**合併分析(混雜,僅供對照)**

| 表徵 | 生物 MoA mAP(95% CI) | 批次混合(→1 越好) | Δbio vs raw |
|---|---|---|---|
| raw | 0.154 [0.130, 0.201] | 3.08 [2.57, 3.57] | — |
| sphered(ZCA) | 0.180 [0.134, 0.246] | 1.85 [1.52, 2.14] | +0.025 [−0.044, +0.098] |
| sphered + Harmony | 0.131 [0.087, 0.193] | 1.04 [0.85, 1.24] | −0.023 [−0.084, +0.045] |

**條件關聯結構保留度(1.0 = 該條件已無法分辨)**

| 表徵 | 細胞株 | 時間點 |
|---|---|---|
| raw | 1.62 | 1.28 |
| sphered(ZCA) | 1.21 | 1.10 |
| sphered + Harmony | 1.11 | **1.02** |

這一張表就是整個 MVP 最重要的產出:Harmony 把批次混合從 3.08 壓到 1.04(近乎完美),
但同一組 embedding 的時間點條件結構同時從 1.28 掉到 1.02——**幾乎完全被抹平**。

措辭必須精確:可以說 Harmony **大幅移除了條件關聯結構**,而這份結構**可能包含真實生物差異**;
由於混雜,**本設計無法量化生物與技術各佔多少**,因此**不能宣稱已證明刪掉真生物**,
只能說批次分數是靠移除條件關聯結構換來的、其中生物成分有受損風險。

**分層評估(plate 此時才是純技術重複;每層 4 盤、約 1,000 wells、260 化合物、26 個合格查詢)**

| 分層 | wells | raw mAP | Δ sphered | Δ sphered+Harmony | raw 批次混合 |
|---|---|---|---|---|---|
| A549-24h | 1,039 | 0.159 | **−0.069** [−0.160, +0.028] | **−0.098** [−0.166, +0.028] | 1.45 |
| A549-48h | 1,040 | 0.187 | +0.020 [−0.030, +0.080] | +0.011 [−0.031, +0.056] | 1.10 |
| U2OS-24h | 983 | 0.154 | −0.011 [−0.046, +0.042] | −0.012 [−0.079, +0.068] | 1.18 |
| U2OS-48h | 1,040 | 0.166 | +0.014 [−0.032, +0.103] | +0.010 [−0.036, +0.088] | 1.57 |

> 舊估計量的影響在這裡看得最清楚:A549-24h 的 sphering Δ,舊版(抽樣平均)報 −0.027,
> 觀測值其實是 **−0.069**;+Harmony 從 −0.051 變成 **−0.098**。
> 每一格的 `bootstrap_mean` 都仍寫在 JSON 裡,大致都是觀測值的一半——舊估計量把效應系統性地往 0 縮。

**誠實解讀**:

1. **條件內的批次效應其實不大**——raw 批次混合 1.10–1.57,遠低於合併時的 3.08。
   合併分析裡看到的「嚴重批次效應」,主要由**條件之間的結構(between-condition structure)**主導;
   由於混雜,這份結構無法被拆解成技術與生物兩部分,因此不可直接稱之為「生物條件差異」。
   這符合 CPJUMP1 是單一 source、單一 site pilot 的預期;真正棘手的跨 source 批次效應,此資料集看不到。
2. **在目前這個探索性 bootstrap 下,八個 Δ 的 95% CI 全部跨過 0**——兩種校正都測不到生物效益,也測不到明確損害。
   但四層**並不齊一**:三層的 |Δ| ≤ 0.02,只有 A549-24h 是 −0.069 / −0.098,區間也只是勉強碰到 0。
   以每層 26 個查詢的功效,無法分辨這是該層真實的損害還是雜訊——但這是最值得追下去的訊號。
3. **統計功效必須誠實揭露**——每層僅 26 個合格查詢(MoA 類別需 ≥2 個化合物),合計 104 個,
   CI 因此偏寬,不足以排除中等幅度的效應。結論是「在此測不到效益」,不是「去偏無效」。
4. **bootstrap 本身的兩個限制**(已修正並揭露):
   - **點估計改用觀測值**。mAP 對「池子裡有哪些化合物」是非線性的,block bootstrap 每次都在換池子,
     所以抽樣平均是有偏估計。舊版把抽樣平均當成 Δ 回報——在 LINCS 上 sphered+Harmony 的抽樣平均是 −0.001,
     但觀測差值其實是 **+0.008,正負號相反**。現在 `value` 一律是觀測值,區間以觀測值重新錨定,
     並額外寫出 `bootstrap_mean` 讓落差看得見。
   - **block 解析度有硬上限**。對 k 個 plate 做 block bootstrap 只有 `C(2k-1, k)` 種相異重抽:
     k=4 時只有 **35 種**。提高抽樣次數(現為 B=1000)只細化權重,不提高解析度。
     每份 JSON 都寫入 `n_distinct_block_draws`,低於 200 時程式會直接印出警告。
   - batch 軸的 bootstrap 也已改為**保留 block 重複次數**(舊版用 `np.unique` 去重,
     等於退化成隨機挑 plate 子集)。重複的 well 不會污染指標,因為 plate-mixing 本來就排除同化合物鄰居。
   - **plate-mixing 的分子與分母現在條件一致**。舊版先取 k+1 個鄰居再剔除同化合物,
     導致部分查詢實際只用到不足 k 個鄰居、彼此權重不等;而期望值卻是對**所有**配對計算的,
     等於分子排除同化合物、分母沒排除。現在改成:先多取一些鄰居、再取**恰好 k 個不同化合物**的鄰居,
     期望值也只在**不同化合物配對**上計算。影響很小(LINCS 上 sphered 3.00 → 2.87、Harmony 1.63 → 1.61,
     raw 不變),結論不變——但這正是本專案在講的事:估計量的定義要對得起它宣稱在量的東西。

### 結論等級

| 區塊 | 定位 |
|---|---|
| LINCS（主 pipeline） | 主要量化結果,可重現、第三方可驗證 |
| LINCS（`mvp/` 雙指標） | 輔助性量化結果 |
| JUMP-Target | **探索性**去偏比較,承認條件混雜;應以分層輸出為準 |
| 影像層 | 案例式端到端展示,非統計驗證 |
