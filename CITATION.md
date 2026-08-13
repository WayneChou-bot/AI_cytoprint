# Data Sources & Citations

Every number, image and figure in this project comes from public, openly licensed data.
This file lists exactly what was used, where it came from, and how to cite it.

If you use this repository, please cite the **upstream data sources** below — they did the
expensive work of generating and releasing these datasets.

---

## 1. Datasets actually used

### 1.1 LINCS Cell Painting (drives the main analysis page)

- **What is used here:** per-plate CellProfiler morphological profiles from batch
  `2016_04_01_a549_48hr_batch1` (A549 cells, 10 plate maps × 5 replicate plates = 50 plates,
  19,200 wells, 561 de-duplicated compounds).
- **Repository:** https://github.com/broadinstitute/lincs-cell-painting
- **Pinned commit:** `da8ae6a3bc103346095d61b4ee02f08fc85a5d98`
- **License:** source code BSD 3-Clause; **data, results and figures CC0 1.0**
- **Citation:**

  > Natoli, T., Way, G., Lu, X., Logan, D., Alimova, M., Hartland, K., Golub, T., Carpenter, A.,
  > Singh, S., Subramanian, A. (2021). *broadinstitute/lincs-cell-painting: Full release of LINCS
  > Cell Painting dataset* (Version v1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5008187

### 1.2 JUMP pilot — CPJUMP1 (drives the real-image section and the JUMP MVP)

- **What is used here:** (a) the raw 5-channel Cell Painting images in `example_images/`
  (9 compounds incl. DMSO control, U2OS cells, plates BR00117010–13, one field per compound);
  (b) well-level profiles + experimental metadata for the 16 compound plates (260 JUMP-Target
  compounds, 2 cell lines, 2 timepoints) used by `mvp/jump_mvp.py`.
- **Repository:** https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted
- **License:** source code BSD 3-Clause; **data, results and figures CC0 1.0**
- **Citation:**

  > Chandrasekaran, S. N., Ackerman, J., Alix, E., Ando, D. M., Arevalo, J., Bennion, M., et al.
  > (2024). *Three million images and morphological profiles of cells treated with matched chemical
  > and genetic perturbations.* **Nature Methods**, 21, 1114–1121.
  > https://doi.org/10.1038/s41592-024-02241-6
  > (preprint: https://doi.org/10.1101/2022.01.05.475090)

### 1.3 Mechanism-of-action (MoA) annotations

- **What is used here:** compound → MoA / target labels, joined by `InChIKey14`.
- **Source:** Broad Drug Repurposing Hub metadata, redistributed in the LINCS repository as
  `metadata/moa/repurposing_info_external_moa_map_resolved.tsv`
- **Citation:**

  > Corsello, S. M., Bittker, J. A., Liu, Z., Gould, J., McCarren, P., Hirschman, J. E., et al.
  > (2017). *The Drug Repurposing Hub: a next-generation drug library and information resource.*
  > **Nature Medicine**, 23, 405–408. https://doi.org/10.1038/nm.4306

### 1.4 Target dataset (referenced, not fully processed here)

- **JUMP-CP `cpg0016`** — ~116k compound and gene perturbations in U2OS
  (~15k unique genes across ORF + CRISPR).
- **Profile index:** https://github.com/jump-cellpainting/datasets
  (`manifests/profile_index.json`)
- **Hosting:** Cell Painting Gallery on the AWS Registry of Open Data —
  free, no AWS account required (`aws s3 ls --no-sign-request s3://cellpainting-gallery/`)
- **Citation:**

  > Weisbart, E., Kumar, A., Arevalo, J., Carpenter, A. E., Cimini, B. A., Singh, S. (2024).
  > *Cell Painting Gallery: an open resource for image-based profiling.*
  > **Nature Methods**, 21, 1775–1778. https://doi.org/10.1038/s41592-024-02399-z

---

## 2. Assay and method references

| Topic | Reference |
|---|---|
| Cell Painting assay | Bray, M.-A. et al. (2016). *Cell Painting, a high-content image-based assay for morphological profiling.* **Nature Protocols** 11, 1757–1774. https://doi.org/10.1038/nprot.2016.105 |
| Optimized protocol | Cimini, B. A. et al. (2023). *Optimizing the Cell Painting assay for image-based profiling.* **Nature Protocols** 18, 1981–2013. https://doi.org/10.1038/s41596-023-00840-9 |
| Profiling toolchain | Serrano, E. et al. (2025). *Reproducible image-based profiling with Pycytominer.* **Nature Methods**. https://github.com/cytomining/pycytominer |
| Retrieval metric (mAP) | Kalinin, A. A. et al. (2025). *A versatile information retrieval framework for evaluating profile strength and similarity.* **Nature Communications** 16, 5181. https://doi.org/10.1038/s41467-025-60306-2 |
| Percent replicating | Way, G. P. et al. (2022). *Morphology and gene expression profiling provide complementary information for mapping cell state.* **Cell Systems** 13, 911–923. |
| Batch-correction benchmark | Arevalo, J. et al. (2024). *Evaluating batch correction methods for image-based cell profiling.* **Nature Communications** 15, 6516. https://doi.org/10.1038/s41467-024-50613-5 |
| Harmony (batch integration) | Korsunsky, I. et al. (2019). *Fast, sensitive and accurate integration of single-cell data with Harmony.* **Nature Methods** 16, 1289–1296. |
| kBET (batch metric) | Büttner, M. et al. (2019). *A test metric for assessing single-cell RNA-seq batch correction.* **Nature Methods** 16, 43–49. |
| Self-supervised profiling | Doron, M. et al. (2023). *Unbiased single-cell morphology with self-supervised vision transformers.* bioRxiv. https://doi.org/10.1101/2023.06.16.545359 |
| Image segmentation tooling | van der Walt, S. et al. (2014). *scikit-image: image processing in Python.* **PeerJ** 2:e453. |

## 3. Biological interpretation cited in the image section

| Claim in the app | Reference |
|---|---|
| Aurora kinase inhibition → failed cytokinesis → polyploid / enlarged nuclei (AMG900) | Payton, M. et al. (2010). *Preclinical evaluation of AMG 900, a novel potent and highly selective pan-aurora kinase inhibitor.* **Cancer Research** 70, 9846–9854. |
| NAMPT inhibition → NAD⁺ depletion → apoptosis (FK-866) | Hasmann, M. & Schemainda, I. (2003). *FK866, a highly specific noncompetitive inhibitor of nicotinamide phosphoribosyltransferase, represents a novel mechanism for induction of tumor cell apoptosis.* **Cancer Research** 63, 7436–7442. |
| Glucocorticoid receptor acts on mitochondria (dexamethasone → reduced Mito signal) | Psarra, A.-M. G. & Sekeris, C. E. (2011) and Lapp, H. E. et al. (2019). See also *Mitochondrial Glucocorticoid Receptors and Their Actions.* **IJMS** 22, 6054. https://doi.org/10.3390/ijms22116054 |

> Note: the mechanism statements above are used to *sanity-check* the image pipeline on
> well-characterised compounds. With one field per compound and no experimental replicates,
> they are illustrative case evidence, not statistical claims.

---

## 4. How to cite this project

If this repository itself is useful to you:

```bibtex
@software{morphoprofile,
  author  = {Chou, Wayne},
  title   = {MorphoProfile: image-based morphological profiling and MoA retrieval
             on open Cell Painting data},
  year    = {2026},
  url     = {https://github.com/USERNAME/morphoprofile},
  note    = {Code MIT-licensed; redistributed data CC0 1.0 from LINCS Cell Painting
             and the JUMP-Cell Painting Consortium}
}
```

Please cite the upstream datasets (Section 1) as well — they are the primary scientific source.
