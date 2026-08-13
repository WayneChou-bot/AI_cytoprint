#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render the biology-vs-batch trade-off (with CIs) from mvp_results.json -> tradeoff.html (+ optional PNG)."""
import json, locale
from pathlib import Path
_p = Path(__file__).resolve().parent/"mvp_results.json"
try:                                   # robust to files written under Windows' default codepage (cp950/Big5)
    D = json.loads(_p.read_text(encoding="utf-8"))
except UnicodeDecodeError:
    D = json.loads(_p.read_text(encoding=locale.getpreferredencoding(False), errors="replace"))
R = D["results"]
COL = {"raw": "#5ad07a", "sphered": "#5a86ff", "sphered+harmony": "#ff9a46"}
LAB = {"raw": "raw (標準化)", "sphered": "sphered (ZCA)", "sphered+harmony": "sphered+Harmony(依 plate)"}
order = [k for k in ["raw", "sphered", "sphered+harmony"] if k in R]

import math
bx = [R[k]["batch_mixing"] for k in order]; by = [R[k]["biology_mAP"] for k in order]
lx = lambda v: math.log10(max(v, 1.0))
xmin, xmax = lx(1.0), lx(max(bx))*1.08 + 0.05
ymax = max(by)*1.25
W, H, PADL, PADB, PADT, PADR = 620, 420, 60, 52, 26, 20
def sx(x): return PADL + (lx(x)-xmin)/(xmax-xmin)*(W-PADL-PADR)
def sy(y): return H-PADB - (y/ymax)*(H-PADB-PADT)
svg = []
# grid + axes (x is log10 of batch-mixing)
for xv in [1, 2, 3, 5, 10, 20, 40]:
    if lx(xv) > xmax: continue
    gx = sx(xv)
    svg.append(f'<line x1="{gx:.0f}" y1="{PADT}" x2="{gx:.0f}" y2="{H-PADB}" stroke="rgba(150,160,210,.08)"/>')
    svg.append(f'<text x="{gx:.0f}" y="{H-PADB+16}" fill="rgba(150,160,190,.6)" font-size="9" text-anchor="middle">{xv}×</text>')
for g in range(5):
    gy = H-PADB - g/4*(H-PADB-PADT); yv = g/4*ymax
    svg.append(f'<text x="{PADL-8}" y="{gy+3:.0f}" fill="rgba(150,160,190,.6)" font-size="9" text-anchor="end">{yv:.2f}</text>')
# ideal corner arrow (low batch, high biology = top-left)
svg.append(f'<text x="{PADL+4}" y="{PADT+4}" fill="rgba(90,208,122,.8)" font-size="10">↖ 理想：批次混合好、生物訊號高</text>')
svg.append(f'<line x1="{sx(1):.0f}" y1="{PADT}" x2="{sx(1):.0f}" y2="{H-PADB}" stroke="rgba(90,208,122,.4)" stroke-dasharray="4 3"/>')
svg.append(f'<text x="{sx(1)+4:.0f}" y="{H-PADB-6}" fill="rgba(90,208,122,.7)" font-size="8">=1 完全混合</text>')
for k in order:
    x, y = sx(R[k]["batch_mixing"]), sy(R[k]["biology_mAP"])
    xl, xh = sx(R[k]["batch_ci"][0]), sx(R[k]["batch_ci"][1])
    yl, yh = sy(R[k]["biology_ci"][0]), sy(R[k]["biology_ci"][1])
    svg.append(f'<line x1="{xl:.0f}" y1="{y:.0f}" x2="{xh:.0f}" y2="{y:.0f}" stroke="{COL[k]}" stroke-width="1.4" opacity=".7"/>')
    svg.append(f'<line x1="{x:.0f}" y1="{yl:.0f}" x2="{x:.0f}" y2="{yh:.0f}" stroke="{COL[k]}" stroke-width="1.4" opacity=".7"/>')
    svg.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="7" fill="{COL[k]}"/>')
    if x > PADL + 0.62*(W-PADL-PADR):   # near right edge → label to the left
        svg.append(f'<text x="{x-11:.0f}" y="{y+3:.0f}" fill="#e9ecf8" font-size="11" text-anchor="end">{LAB[k]}</text>')
    else:
        svg.append(f'<text x="{x+11:.0f}" y="{y+3:.0f}" fill="#e9ecf8" font-size="11">{LAB[k]}</text>')
svg.append(f'<text x="{(W)/2:.0f}" y="{H-6}" fill="rgba(200,208,230,.7)" font-size="10" text-anchor="middle">批次混合 batch-mixing（→1 越好）</text>')
svg.append(f'<text x="14" y="{H/2:.0f}" fill="rgba(200,208,230,.7)" font-size="10" transform="rotate(-90 14 {H/2:.0f})" text-anchor="middle">生物訊號 MoA mAP（越高越好）</text>')
svg_s = f'<svg viewBox="0 0 {W} {H}" width="100%" style="background:#0a0a12;border:1px solid rgba(150,160,210,.13);border-radius:12px">'+"".join(svg)+"</svg>"

rows = "".join(f'<tr><td style="color:{COL[k]}">●</td><td>{LAB[k]}</td>'
               f'<td>{R[k]["biology_mAP"]} <span class="ci">[{R[k]["biology_ci"][0]}, {R[k]["biology_ci"][1]}]</span></td>'
               f'<td>{R[k]["batch_mixing"]} <span class="ci">[{R[k]["batch_ci"][0]}, {R[k]["batch_ci"][1]}]</span></td></tr>' for k in order)
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"/>
<title>去偏雙指標權衡（真實 · 含 CI）</title><style>
body{{background:#06060c;color:#e9ecf8;font-family:ui-monospace,Menlo,Consolas,monospace;max-width:820px;margin:0 auto;padding:20px}}
h1{{font-size:18px}} .sub{{color:#828aa8;font-size:12px;line-height:1.7}}
table{{width:100%;border-collapse:collapse;font-size:12px;margin-top:14px}}
th,td{{padding:7px 6px;border-bottom:1px solid rgba(150,160,210,.13);text-align:left}}
th{{color:#828aa8;font-size:9px;letter-spacing:.1em;text-transform:uppercase;font-weight:400}}
.ci{{color:#828aa8;font-size:10px}} .note{{background:rgba(255,154,70,.07);border:1px solid rgba(255,154,70,.25);border-radius:10px;padding:12px;font-size:12px;color:#b9c0d8;line-height:1.8;margin-top:14px}}
b{{color:#e9ecf8}}</style></head><body>
<h1>去偏方法的雙指標權衡 · {D['dataset']}</h1>
<div class="sub">{D['n_wells']:,} wells · {D['n_plates']} plates(=block) · {D['pca_dims']}-D 共同空間 · plate-block bootstrap {D['bootstrap_draws']} draws（生物軸）。
生物軸=MoA 檢索 mAP（越高越好）；批次軸=plate-mixing enrichment（→1 越好，>1 表仍有批次結構）。</div>
{svg_s}
<table><thead><tr><th></th><th>方法</th><th>生物訊號 mAP（95% CI）</th><th>批次混合（95% CI）</th></tr></thead><tbody>{rows}</tbody></table>
<div class="note"><b>誠實解讀（已修掉極端值污染後重算）：</b>三種方法的<b>生物訊號幾乎相同</b>
（mAP raw {R['raw']['biology_mAP']} / sphered {R['sphered']['biology_mAP']} / +Harmony {R['sphered+harmony']['biology_mAP']}；
兩兩 paired block Δ 的 95% CI <b>都跨 0</b>，即無顯著差異）。批次軸上 <b>Harmony 明顯改善混合</b>
（{R['raw']['batch_mixing']}× → {R['sphered+harmony']['batch_mixing']}×，往 1 靠），且<b>沒有付出可測的生物代價</b>；
單獨球化則略微變差。<br><br>
教訓：<b>去偏必須同時看「生物保存」與「批次移除」兩軸、用 plate-block bootstrap 與 paired Δ CI</b>。
本結果為<b>單一中心</b> pilot、批次效應本就小,故差異多不顯著,<b>不能一般化</b>；同一引擎（jump_mvp.py）在多來源 JUMP 以 <b>source</b> 為 block、批次效應大時,才是校正真正該發揮的場景。
（對照:先前「Harmony 兩軸皆輸」是<b>資料污染的假象</b>,清理後結論反轉——正說明資料品質與嚴謹評估缺一不可。）</div>
</body></html>"""
out = Path(__file__).resolve().parent/"tradeoff.html"
out.write_text(html, encoding="utf-8"); print("wrote", out)
