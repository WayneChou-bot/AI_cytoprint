# Demo images — for the in-browser upload feature

Two real Cell Painting DNA-channel images (1080×1080 PNG), chosen because they give the
strongest possible visual and numerical contrast for a live demo.

| File | Compound | What the viewer sees | What the tool should report |
|---|---|---|---|
| `01_control_DMSO_DNA.png` | **DMSO** (vehicle control) | Many small, evenly-sized nuclei, densely packed | **93 nuclei**, median area ≈ 984 px → nearest reference: **DMSO** ✅ |
| `02_treated_AMG900_DNA.png` | **AMG900** (Aurora kinase inhibitor) | Far fewer nuclei, visibly **enlarged and brighter** | **35 nuclei**, median area ≈ 1734 px → nearest reference: **AMG900** ✅ |

The in-browser detector lands on the correct compound for both images. Its counts differ slightly
from the Python pipeline (95 / 38 nuclei) because the browser uses a plain Otsu threshold without
watershed splitting of touching nuclei — worth mentioning on camera as an honest caveat.

Both are single fields from plates BR00117010–13 (U2OS cells), JUMP pilot CPJUMP1,
released under CC0 1.0. Source and citation: see [`../CITATION.md`](../CITATION.md).

## Suggested 30-second demo script

1. Open the app, scroll to **"Upload your own image"**.
2. Upload `01_control_DMSO_DNA.png` → the canvas outlines each nucleus in yellow.
   Point out the count and median nucleus area. *"This is the untreated control."*
3. Upload `02_treated_AMG900_DNA.png` → **the same detector, no retraining.**
   *"Half the cells are gone, and the surviving nuclei are ~50 % larger."*
4. Land the point: *"AMG900 is an Aurora kinase inhibitor — it blocks cytokinesis, so cells
   keep duplicating DNA without dividing. Enlarged, polyploid nuclei are exactly the expected
   phenotype. The tool recovered the mechanism from pixels alone."*
5. (Optional) Scroll up to the real-image viewer and switch channels
   (composite → segmentation → DNA/RNA/ER/AGP/Mito) to show it is genuine 5-channel data.

## Notes for recording

- The browser-side analysis is intentionally **simplified** (grayscale, Otsu threshold,
  connected components, three summary statistics). Say so on camera — the full 59-feature
  fingerprint is the Python pipeline (`images/image_pipeline.py`). Being upfront about this
  is more impressive than glossing over it.
- Counts are scaled assuming a full 1080×1080 field, which is what these files are. Cropping
  or resizing before upload will change the reported numbers.
- Nothing is uploaded anywhere — all processing happens in the browser on a `<canvas>`.
  Worth mentioning if privacy comes up.
