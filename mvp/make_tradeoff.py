#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Render the biology-vs-batch trade-off (with CIs) as a self-contained HTML page.

Usage:
    python make_tradeoff.py                      # LINCS results -> tradeoff.html (EN) + tradeoff.zh-TW.html
    python make_tradeoff.py jump_mvp_results.json   # same, from the real JUMP-Target run

The English page is the one embedded in README.md; the Chinese one in README.zh-TW.md.
PNGs are produced by screenshotting these pages (any headless browser).
"""
import json, locale, math, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE/(sys.argv[1] if len(sys.argv) > 1 else "mvp_results.json")
def _load(path):
    """Robust to JSON written under a Windows codepage (cp950/Big5) rather than UTF-8."""
    raw = Path(path).read_bytes()
    for enc in ("utf-8", "cp950", "big5", locale.getpreferredencoding(False), "latin-1"):
        try:
            return json.loads(raw.decode(enc))
        except (UnicodeDecodeError, LookupError):
            continue
    raise SystemExit(f"cannot decode {path}")
D = _load(SRC)
R = D["results"]
order = [k for k in ["raw", "sphered", "sphered+harmony"] if k in R]
COL = {"raw": "#5ad07a", "sphered": "#5a86ff", "sphered+harmony": "#ff9a46"}

T = {
 "en": dict(
   lab={"raw": "raw (standardized)", "sphered": "sphered (ZCA)", "sphered+harmony": "sphered + Harmony (by plate)"},
   title="Dual-metric trade-off of batch correction",
   sub="{wells:,} wells · {plates} plates (= block) · {dims}-D common space · plate-block bootstrap {draws} draws (biology axis).<br>"
       "Biology axis = MoA-retrieval mAP (higher is better); batch axis = plate-mixing enrichment (→1 is better; &gt;1 means batch structure remains).",
   ideal="↖ ideal: well mixed, biology preserved", unit="=1 fully mixed",
   xlab="batch mixing (→1 is better)", ylab="biology signal — MoA mAP (higher is better)",
   th=["", "Method", "Biology mAP (95% CI)", "Batch mixing (95% CI)"],
   note="<b>Honest reading (recomputed after fixing an extreme-value contamination bug):</b> the three methods have "
        "<b>almost identical biology signal</b> (mAP raw {r} / sphered {s} / +Harmony {h}; every pairwise paired-block Δ "
        "has a 95% CI that <b>crosses zero</b>, i.e. no significant difference). On the batch axis, "
        "<b>Harmony clearly improves mixing</b> ({br}× → {bh}×, toward 1) at <b>no measurable biology cost</b>; "
        "sphering alone is slightly worse.<br><br>"
        "Lesson: <b>de-biasing must be judged on both axes at once — biology retention and batch removal — using a "
        "plate-block bootstrap and paired Δ CIs.</b> This is a <b>single-centre</b> pilot where batch effects are small "
        "by construction, so most differences are not significant and <b>should not be generalized</b>. The same engine "
        "(jump_mvp.py) on multi-source JUMP, with <b>source</b> as the block, is where correction really earns its keep.<br><br>"
        "(For contrast: an earlier version reported \"Harmony loses on both axes\" — that was <b>an artefact of contaminated "
        "data</b>. After cleaning, the conclusion reversed, which is exactly why data quality and rigorous evaluation are "
        "both indispensable.)",
   lang="en", pagetitle="Batch-correction trade-off (real data, with CIs)"),
 "zh": dict(
   lab={"raw": "raw（標準化）", "sphered": "sphered（ZCA）", "sphered+harmony": "sphered + Harmony（依 plate）"},
   title="去偏方法的雙指標權衡",
   sub="{wells:,} wells · {plates} plates（= block）· {dims}-D 共同空間 · plate-block bootstrap {draws} draws（生物軸）。<br>"
       "生物軸 = MoA 檢索 mAP（越高越好）；批次軸 = plate-mixing enrichment（→1 越好，&gt;1 表仍有批次結構）。",
   ideal="↖ 理想：批次混合好、生物訊號高", unit="=1 完全混合",
   xlab="批次混合 batch mixing（→1 越好）", ylab="生物訊號 MoA mAP（越高越好）",
   th=["", "方法", "生物訊號 mAP（95% CI）", "批次混合（95% CI）"],
   note="<b>誠實解讀（已修掉極端值污染後重算）：</b>三種方法的<b>生物訊號幾乎相同</b>"
        "（mAP raw {r} / sphered {s} / +Harmony {h}；兩兩 paired block Δ 的 95% CI <b>都跨 0</b>，即無顯著差異）。"
        "批次軸上 <b>Harmony 明顯改善混合</b>（{br}× → {bh}×，往 1 靠），且<b>沒有付出可測的生物代價</b>；單獨球化則略微變差。<br><br>"
        "教訓：<b>去偏必須同時看「生物保存」與「批次移除」兩軸，並用 plate-block bootstrap 與 paired Δ CI。</b>"
        "本結果為<b>單一中心</b> pilot、批次效應本就小，故差異多不顯著，<b>不能一般化</b>；同一引擎（jump_mvp.py）"
        "在多來源 JUMP 以 <b>source</b> 為 block、批次效應大時，才是校正真正該發揮的場景。<br><br>"
        "（對照：先前「Harmony 兩軸皆輸」是<b>資料污染的假象</b>，清理後結論反轉——正說明資料品質與嚴謹評估缺一不可。）",
   lang="zh-Hant", pagetitle="去偏雙指標權衡（真實 · 含 CI）"),
}

W, H, PADL, PADB, PADT, PADR = 620, 430, 62, 54, 30, 22
bx = [R[k]["batch_mixing"] for k in order]
by = [R[k]["biology_mAP"] for k in order]
hi_ci = [R[k].get("biology_ci", [0, R[k]["biology_mAP"]])[1] for k in order]
lx = lambda v: math.log10(max(v, 1.0))
xmin, xmax = lx(1.0), lx(max(bx))*1.10 + 0.05
ymax = max(hi_ci + by)*1.22
sx = lambda x: PADL + (lx(x)-xmin)/(xmax-xmin)*(W-PADL-PADR)
sy = lambda y: H-PADB - (y/ymax)*(H-PADB-PADT)

def build(cfg):
    LAB = cfg["lab"]
    svg = []
    for xv in [1, 1.5, 2, 3, 5, 10, 20, 40]:
        if lx(xv) > xmax: continue
        gx = sx(xv)
        svg.append(f'<line x1="{gx:.0f}" y1="{PADT}" x2="{gx:.0f}" y2="{H-PADB}" stroke="rgba(150,160,210,.08)"/>')
        lbl = f"{xv:g}×"
        svg.append(f'<text x="{gx:.0f}" y="{H-PADB+16}" fill="rgba(150,160,190,.6)" font-size="9" text-anchor="middle">{lbl}</text>')
    for g in range(5):
        gy = H-PADB - g/4*(H-PADB-PADT); yv = g/4*ymax
        svg.append(f'<line x1="{PADL}" y1="{gy:.0f}" x2="{W-PADR}" y2="{gy:.0f}" stroke="rgba(150,160,210,.05)"/>')
        svg.append(f'<text x="{PADL-8}" y="{gy+3:.0f}" fill="rgba(150,160,190,.6)" font-size="9" text-anchor="end">{yv:.2f}</text>')
    svg.append(f'<text x="{PADL+4}" y="{PADT-8}" fill="rgba(90,208,122,.85)" font-size="10">{cfg["ideal"]}</text>')
    svg.append(f'<line x1="{sx(1):.0f}" y1="{PADT}" x2="{sx(1):.0f}" y2="{H-PADB}" stroke="rgba(90,208,122,.4)" stroke-dasharray="4 3"/>')
    svg.append(f'<text x="{sx(1)+4:.0f}" y="{H-PADB-6}" fill="rgba(90,208,122,.7)" font-size="8">{cfg["unit"]}</text>')

    # ---- markers + CI whiskers ----
    for k in order:
        x, y = sx(R[k]["batch_mixing"]), sy(R[k]["biology_mAP"])
        if "batch_ci" in R[k]:
            xl, xh = sx(R[k]["batch_ci"][0]), sx(R[k]["batch_ci"][1])
            svg.append(f'<line x1="{xl:.0f}" y1="{y:.0f}" x2="{xh:.0f}" y2="{y:.0f}" stroke="{COL[k]}" stroke-width="1.4" opacity=".65"/>')
        if "biology_ci" in R[k]:
            yl, yh = sy(R[k]["biology_ci"][0]), sy(R[k]["biology_ci"][1])
            svg.append(f'<line x1="{x:.0f}" y1="{yl:.0f}" x2="{x:.0f}" y2="{yh:.0f}" stroke="{COL[k]}" stroke-width="1.4" opacity=".65"/>')
        svg.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="7" fill="{COL[k]}"/>')

    # ---- labels with collision avoidance (placed above the point / its CI whisker) ----
    CHW, LH = 6.3, 15          # approx char width @ font-size 11 mono, line height
    placed = []                # (x0, x1, y)
    # top-most points first so their labels sit closest to their marker
    for k in sorted(order, key=lambda k: -R[k]["biology_mAP"]):
        x, y = sx(R[k]["batch_mixing"]), sy(R[k]["biology_mAP"])
        top = sy(R[k]["biology_ci"][1]) if "biology_ci" in R[k] else y
        w = len(LAB[k])*CHW
        ty = top - 10
        x0, x1 = x - w/2, x + w/2
        if x0 < PADL: x0, x1 = PADL, PADL + w                 # keep inside the plot
        if x1 > W-PADR: x0, x1 = W-PADR-w, W-PADR
        guard = 0
        while any(not (x1 < px0-6 or x0 > px1+6) and abs(ty-py) < LH for px0, px1, py in placed) and guard < 12:
            ty -= LH; guard += 1
        ty = max(ty, PADT+10)
        placed.append((x0, x1, ty))
        cx = (x0+x1)/2
        svg.append(f'<line x1="{cx:.0f}" y1="{ty+4:.0f}" x2="{x:.0f}" y2="{top-3:.0f}" stroke="{COL[k]}" stroke-width="1" opacity=".45"/>')
        svg.append(f'<text x="{cx:.0f}" y="{ty:.0f}" fill="#e9ecf8" font-size="11" text-anchor="middle">{LAB[k]}</text>')

    svg.append(f'<text x="{W/2:.0f}" y="{H-6}" fill="rgba(200,208,230,.75)" font-size="10" text-anchor="middle">{cfg["xlab"]}</text>')
    svg.append(f'<text x="14" y="{H/2:.0f}" fill="rgba(200,208,230,.75)" font-size="10" transform="rotate(-90 14 {H/2:.0f})" text-anchor="middle">{cfg["ylab"]}</text>')
    svg_s = (f'<svg viewBox="0 0 {W} {H}" width="100%" '
             'style="background:#0a0a12;border:1px solid rgba(150,160,210,.13);border-radius:12px">' + "".join(svg) + "</svg>")

    rows = "".join(
        f'<tr><td style="color:{COL[k]}">●</td><td>{LAB[k]}</td>'
        f'<td>{R[k]["biology_mAP"]} <span class="ci">[{R[k].get("biology_ci",["",""])[0]}, {R[k].get("biology_ci",["",""])[1]}]</span></td>'
        f'<td>{R[k]["batch_mixing"]} <span class="ci">[{R[k].get("batch_ci",["",""])[0]}, {R[k].get("batch_ci",["",""])[1]}]</span></td></tr>'
        for k in order)
    sub = cfg["sub"].format(wells=D.get("n_wells", 0), plates=D.get("n_plates", "—"),
                            dims=D.get("pca_dims", 50), draws=D.get("bootstrap_draws", "—"))
    g = lambda k, f: R[k][f] if k in R else "—"
    note = cfg["note"].format(r=g("raw", "biology_mAP"), s=g("sphered", "biology_mAP"),
                              h=g("sphered+harmony", "biology_mAP"),
                              br=g("raw", "batch_mixing"), bh=g("sphered+harmony", "batch_mixing"))
    return f"""<!DOCTYPE html><html lang="{cfg['lang']}"><head><meta charset="utf-8"/>
<title>{cfg['pagetitle']}</title><style>
body{{background:#06060c;color:#e9ecf8;font-family:ui-monospace,Menlo,Consolas,monospace;max-width:820px;margin:0 auto;padding:20px}}
h1{{font-size:17px;line-height:1.5}} .sub{{color:#828aa8;font-size:12px;line-height:1.7;margin:8px 0 14px}}
table{{width:100%;border-collapse:collapse;font-size:12px;margin-top:14px}}
th,td{{padding:7px 6px;border-bottom:1px solid rgba(150,160,210,.13);text-align:left}}
th{{color:#828aa8;font-size:9px;letter-spacing:.1em;text-transform:uppercase;font-weight:400}}
.ci{{color:#828aa8;font-size:10px}}
.note{{background:rgba(255,154,70,.07);border:1px solid rgba(255,154,70,.25);border-radius:10px;padding:12px;
  font-size:12px;color:#b9c0d8;line-height:1.8;margin-top:14px}}
b{{color:#e9ecf8}}</style></head><body>
<h1>{cfg['title']} · {D['dataset']}</h1>
<div class="sub">{sub}</div>
{svg_s}
<table><thead><tr>{''.join(f'<th>{h}</th>' for h in cfg['th'])}</tr></thead><tbody>{rows}</tbody></table>
<div class="note">{note}</div>
</body></html>"""

stem = "tradeoff" if SRC.name == "mvp_results.json" else SRC.stem.replace("_results", "")
for suffix, key in [("", "en"), (".zh-TW", "zh")]:
    out = HERE/f"{stem}{suffix}.html"
    out.write_text(build(T[key]), encoding="utf-8")
    print("wrote", out)
