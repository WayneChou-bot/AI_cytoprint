# MorphoProfile — Image-Based Morphological Profiling & MoA Retrieval

[English](README.md) · [繁體中文](README.zh-TW.md)

**Predict a compound's mechanism of action from what cells *look like* — no molecular structure required.**

An end-to-end, reproducible pipeline on open [Cell Painting](https://doi.org/10.1038/nprot.2016.105)
data: raw 5-channel microscopy images → segmentation → morphological fingerprints → nearest-neighbour
MoA retrieval, plus a rigorously evaluated batch-correction comparison.

Everything is computed offline from public data by the scripts in this repo, then embedded in a
single self-contained `index.html`. **No fabricated or placeholder numbers.**

![trade-off](mvp/tradeoff.png)

---

## What this is

Most drug-discovery ML competes on *molecule → property* prediction — a crowded benchmark space.
This project takes the other input: **the image**. It builds morphological fingerprints from cell
images and uses **representation learning + retrieval** to infer mechanism, which is a far less
crowded niche for an individual developer and uses computer-vision skills directly.

The work is split into three layers, and the README is honest about what each one proves.

| Layer | What it does | Scale | Evidence level |
|---|---|---|---|
| **1. Image → fingerprint** (`images/`) | Real 5-channel images → my own nucleus segmentation → 59 morphological features → z-score vs DMSO | 9 compounds, 701 single cells | **Case study** — demonstrates the pipeline works and is biologically plausible. No statistical validation. |
| **2. Profile analysis** (`pipeline.py`) | Pre-computed CellProfiler profiles → cleaning → consensus → UMAP, MoA retrieval, candidate-activity screen | 19,200 wells, 561 compounds, 587 features | **Main quantitative result** — reproducible, third-party verifiable. |
| **3. Batch-correction evaluation** (`mvp/`) | Dual metrics + plate-block bootstrap + paired ΔCI | LINCS (50 plates); JUMP-Target (16 plates, 4,102 wells) | LINCS: supporting quantitative. JUMP: **exploratory** — plate is confounded with condition (see below). |

Being explicit about this is deliberate: not every section is a formal validation, and saying so is what
separates a trustworthy analysis from a demo dressed up as one.

Layer 2 uses profiles rather than pixels because that is the field's standard division of labour —
the Broad Institute already ran CellProfiler across all of JUMP and published the outputs. Downloading
the ~127 TB of raw JUMP images is neither necessary nor advisable. Layer 1 exists so the project can
still demonstrate the image-processing half end to end.

---

## Highlights

### Image pipeline recovers known biology

Using **my own 59 features** (not the official CellProfiler profiles), on real CPJUMP1 images:

| Compound | Mechanism | Median nucleus area | vs DMSO | Cell count |
|---|---|---|---|---|
| DMSO | control | 934 | — | 95 |
| **AMG900** | Aurora kinase inhibitor | **1389** | **+49 %** | 38 (40 %) |
| **FK-866** | NAMPT inhibitor | **451** | **−52 %** | 36 (38 %) |
| quinidine | ion-channel blocker | 928 | ≈ 0 % | 94 (99 %) |

The two strongest phenotypes move nuclear size in **opposite directions**, each consistent with its
known mechanism (Aurora inhibition → failed cytokinesis → polyploid nuclei; NAD⁺ depletion →
apoptotic chromatin condensation), while the weak-phenotype control barely moves. Suggestively, the
fingerprint is not purely a cell-death readout either: dexamethasone loses few cells yet ranks third,
its signature dominated by a drop in mitochondrial signal — consistent with glucocorticoid-receptor
action on mitochondria.

> **Scope of this layer:** an *illustrative end-to-end case study* using **one field per compound**.
> It demonstrates pipeline functionality and biological plausibility — **not statistical validation or
> generalization**. The individual cells are observations within a single field, not independent
> experimental replicates, so no significance testing is performed or implied.

### MoA retrieval on 561 compounds

Nearest-neighbour retrieval over standardized consensus profiles:
**Top-1 = 9.0 %**, **Top-5 recall = 14.9 %**, versus a self-excluded random baseline of **0.71 %**
(**≈ 12.7× chance**, n = 377 evaluable compounds).

Modest but clearly above chance — and honest about it. That gap is exactly the motivation for better
representations and validated batch correction.

### Batch correction: measured, not assumed

Correction is evaluated on **two axes simultaneously** (you can trivially win one by wrecking the
other), in a common 50-D space, with a **plate-block bootstrap** and **paired** Δ confidence intervals.

**Exploratory** JUMP-Target result (`mvp/jump_mvp_results.json`, 16 plates, 4,102 treated wells).

In CPJUMP1 **all 16 compound plates each contain a single cell line and a single timepoint**, so plate is
perfectly nested within `cell_line × timepoint` (verified in the metadata). A plate-block correction
therefore cannot be separated from erasing genuine U2OS-vs-A549 and 24 h-vs-48 h differences — and a
*pooled* MoA retrieval cannot detect that loss, because it pools across those very conditions.

Pooled (confounded — reported for reference only):

| Representation | Biology (MoA mAP, 95 % CI) | Batch mixing (→1 is better) |
|---|---|---|
| raw | 0.154 [0.132, 0.213] | 3.17 [2.30, 3.19] |
| sphered (ZCA) | 0.180 [0.114, 0.240] | 1.88 [1.66, 1.87] |
| sphered + Harmony | 0.131 [0.107, 0.219] | 1.04 [1.01, 1.11] |

> **Harmony improved plate mixing, while pooled MoA retrieval showed no measurable loss. Because plate
> is confounded with cell line and timepoint in CPJUMP1, this result cannot distinguish removal of
> technical variation from removal of condition-specific biology.**

**Condition-signal retention** measures exactly that (kNN same-label enrichment on the same embeddings;
1.0 means the difference is gone):

| Representation | Cell line | Timepoint |
|---|---|---|
| raw | 1.62 | 1.28 |
| sphered (ZCA) | 1.21 | 1.10 |
| sphered + Harmony | 1.11 | **1.02** |

Harmony's near-perfect batch mixing (3.17 → 1.04) coincides with the timepoint signal being **essentially
erased** (1.28 → 1.02) and the cell-line signal being more than halved in excess-over-chance terms
(1.62 → 1.11). A large share of that "batch removal" was paid for with real biology.

**Stratified** — evaluated within each `cell_line × timepoint` stratum, where plate *is* a genuine technical
replicate (4 plates, ~1,000 wells, 260 compounds, 26 eligible MoA queries each):

| Stratum | raw mAP | Δ sphered (95 % CI) | Δ sphered+Harmony (95 % CI) | raw batch mixing |
|---|---|---|---|---|
| A549-24h | 0.159 | −0.027 [−0.123, +0.062] | −0.051 [−0.131, +0.068] | 1.46 |
| A549-48h | 0.187 | +0.006 [−0.046, +0.035] | 0.000 [−0.047, +0.047] | 1.11 |
| U2OS-24h | 0.154 | −0.011 [−0.045, +0.043] | −0.005 [−0.072, +0.075] | 1.19 |
| U2OS-48h | 0.166 | +0.001 [−0.042, +0.093] | +0.002 [−0.041, +0.082] | 1.59 |

Two things follow. First, **the within-condition batch effect is modest**: raw batch mixing is 1.11–1.59,
not the 3.17 seen when pooling — so most of the apparent batch effect *was* the confounded biological
condition. That is what one expects from CPJUMP1, a single-source, single-site pilot. Second, **all eight
Δ intervals cross zero**: at this design and sample size (26 eligible queries per stratum, 104 in total)
neither correction produces a measurable biological benefit, nor clear harm. The intervals are wide and
cannot rule out moderate effects — "not detectable here" is not "does not work".

`mvp/jump_mvp.py` performs all three steps automatically: it detects and reports the nesting, runs the
condition-retention diagnostic, and evaluates within strata, writing `"status": "exploratory"` and an
explicit `caveat` field into the results JSON.

> A note on scientific integrity: an earlier version of this project reported the *opposite*
> conclusion about batch correction. That was traced to a data-quality bug — a winsorization step was
> used to select feature names but never written back to the analysis matrix, letting ~176 near-zero-MAD
> features with values up to 1e19 dominate every distance. After the fix, the conclusion reversed.
> The fix, the invalidation mechanism, and the corrected numbers are all in the repo.

---

## Quick start

```bash
git clone https://github.com/WayneChou-bot/AI_cytoprint.git
cd AI_cytoprint

# environment (either one)
conda env create -f environment.yml && conda activate morpho
# or: pip install -r requirements.txt

# (optional) image layer — needs the upstream image repo, ~200 MB
git clone https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted
python images/image_pipeline.py 2021_Chandrasekaran_submitted

# main analysis + page  (downloads ~50 MB of LINCS profiles on first run, cached in data/)
python pipeline.py
python build_page.py            # -> index.html
python test_consistency.py      # asserts the page's tables and headline metrics agree
```

Then open `index.html` in a browser. No server, no build step, no external requests at runtime.

**Note on the image data:** the TIFFs in `example_images/` come down with a normal `git clone` —
only `.gz` files in that repository use Git LFS.

### Deploying

`index.html` is fully self-contained (all CSS, JS, data and images inlined as base64), so any static
host works. On **Vercel**: import the repository and deploy — the root `index.html` is served
automatically, no framework preset or configuration needed.

---

## Repository layout

```
index.html                  Self-contained interactive app (generated — do not edit by hand)
index_template.html         Template with __DATA__ / __IMGDATA__ placeholders (edit this)
build_page.py               Injects results into the template -> index.html

pipeline.py                 Single source of truth for the profile analysis
test_consistency.py         Verifies neighbour tables == headline metrics
run_on_JUMP.py              Points the profile pipeline at full JUMP cpg0016 (S3, anonymous)

images/
  image_pipeline.py         Images -> segmentation -> 59 features -> fingerprints
  out/webimages.json        Channel images, segmentation overlays, fingerprints (generated)
  out/per_cell_features.csv 701 cells x 59 features (generated, inspectable)

mvp/
  rigorous_eval.py          Dual-metric evaluation engine + plate-block bootstrap
  jump_mvp.py               Same engine on real JUMP-Target (CPJUMP1)
  make_tradeoff.py          Renders the trade-off figure
  mvp_results.json          LINCS results
  jump_mvp_results.json     Real JUMP-Target results
  README_MVP.md             Method notes and limitations

demo/                       Two images for demoing the in-browser upload feature
web/webdata.json            Computed results for the main page (generated)
CITATION.md                 Full data sources and citations
```

---

## Method notes

**Feature sanitization (`pipeline.py`).** Non-finite values → NaN; features whose max |value| exceeds
100 are dropped as MAD-blowup artifacts; the rest are winsorized to ±15 **and written back**; then
variance and correlation filtering. Downstream code reads only the sanitized matrix, and
`data/manifest.json` records the commit, cleaning rules and a version tag so a stale cache cannot be
silently reused.

**Metrics.** Retrieval uses cosine similarity on standardized per-compound consensus profiles. The
random baseline **excludes the query itself** and is matched to eligible-query class sizes:
`mean_i (size(MoA_i) − 1) / (P − 1)`. The neighbour tables shown in the app and the headline metrics
are produced by the *same* computation — `test_consistency.py` asserts a third party can recompute
one from the other exactly.

**Candidate activity.** A heuristic screen (consensus L2 distance from the DMSO centroid exceeding the
95th percentile of a DMSO-only null), deliberately labelled *candidate active* rather than
*significant* — there is no plate-matched null, no whitening, and no per-compound FDR. The rigorous
upgrade is replicate mAP with permutation p-values and FDR.

**Channel mapping.** `ch5 = DNA` and `ch1 = Mito` were determined empirically (nuclear morphology and
segmentation for DNA; strongest nuclear exclusion and absent nucleolar contrast for Mito). `ch2/3/4 =
AGP/RNA/ER` are inferred from those anchors by acquisition wavelength order and are **not** confirmed
against the authoritative table, which lives in the LFS-tracked `load_data_csv/*.csv.gz`. Verify with:

```bash
git lfs pull --include="load_data_csv/2020_11_04_CPJUMP1/BR00117010/load_data.csv.gz"
```

---

## Limitations

- **One field per compound, no experimental replicates** in the image layer — the mechanism table is a
  sanity check on known compounds, not a statistical result. Cells within a field are not independent
  replicates; phenotype strength is not shown to be independent of cell density or cell death.
- **The JUMP-Target comparison is exploratory**: plate is perfectly confounded with `cell_line × timepoint`
  in CPJUMP1, so the pooled result cannot separate technical from condition-specific variation. Use the
  stratified output instead.
- **59 hand-built features** are a simplified re-implementation of CellProfiler (1,747 features),
  intended to demonstrate the end-to-end capability, not to replace it.
- The in-browser upload feature runs a **simplified analysis** (single grayscale channel, Otsu,
  connected components, three summary statistics) — it is not the 59-D Python fingerprint.
- LINCS is **single-centre**, so its batch effects are small by construction; results there should not
  be generalized to multi-source JUMP.
- Moving to full JUMP `cpg0016` is **not** a URL swap: each subset is pre-processed differently
  (compound = feature-select + Harmony; CRISPR = sphering + Harmony + PCA; ORF = sphering + Harmony),
  and all-pairs cosine over ~116k compounds is a >100 GB dense matrix — use FAISS or chunked top-k.
- mAP is currently query-weighted; macro-averaged-per-MoA should be reported alongside it.

## Roadmap

1. Replace the heuristic activity screen with replicate mAP + permutation p-values + FDR.
2. Extend the image layer beyond one field per compound so statistics become possible.
3. Self-supervised image embeddings (DINO-style) to replace hand-built features — needs GPU and far
   more images than the 9 fields used here.
4. Scale retrieval to full JUMP `cpg0016` with FAISS / chunked top-k.
5. An LLM layer that *explains retrieved evidence only* — never asserts an unretrieved mechanism.

---

## License & attribution

Source code is **MIT** licensed (see [LICENSE](LICENSE)).

Redistributed scientific data — the microscopy images in `demo/`, the images embedded in
`index.html`, and the derived feature tables — originates from the Broad Institute and the
JUMP-Cell Painting Consortium, whose upstream repositories dual-license **data, results and figures
under CC0 1.0**. That data remains CC0 here.

**Please cite the upstream datasets.** Full details in **[CITATION.md](CITATION.md)**; the primary ones are:

- Chandrasekaran, S. N. et al. (2024). *Three million images and morphological profiles of cells
  treated with matched chemical and genetic perturbations.* **Nature Methods** 21, 1114–1121.
  https://doi.org/10.1038/s41592-024-02241-6
- Natoli, T. et al. (2021). *broadinstitute/lincs-cell-painting: Full release of LINCS Cell Painting
  dataset.* Zenodo. https://doi.org/10.5281/zenodo.5008187
- Weisbart, E. et al. (2024). *Cell Painting Gallery: an open resource for image-based profiling.*
  **Nature Methods** 21, 1775–1778. https://doi.org/10.1038/s41592-024-02399-z
- Bray, M.-A. et al. (2016). *Cell Painting, a high-content image-based assay for morphological
  profiling.* **Nature Protocols** 11, 1757–1774. https://doi.org/10.1038/nprot.2016.105

Data access is free and requires no AWS account:
`aws s3 ls --no-sign-request s3://cellpainting-gallery/`
