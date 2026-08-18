#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Inject computed results into index_template.html -> index.html (repo root, for static hosting)
and web/MorphoProfile.html (local copy). Both are the same single self-contained file."""
from pathlib import Path
tpl = Path("index_template.html").read_text(encoding="utf-8")
data = Path("web/webdata.json").read_text(encoding="utf-8").replace("</", "<\\/")  # safe inside <script>
imgp = Path("images/out/webimages.json")
if imgp.exists():
    imgdata = imgp.read_text(encoding="utf-8").replace("</", "<\\/")
else:
    imgdata = '{"compounds":{},"source":"(not generated: run images/image_pipeline.py first)","channel_note":"","feature_names":[]}'
    print("WARNING: images/out/webimages.json missing -> image section will be empty. Run: python images/image_pipeline.py")
jp = Path("mvp/jump_mvp_results.json")
if jp.exists():
    jumpdata = jp.read_text(encoding="utf-8").replace("</", "<\\/")
else:
    jumpdata = '{}'
    print("WARNING: mvp/jump_mvp_results.json missing -> JUMP section will be empty. Run: python mvp/jump_mvp.py <repo>")
html = tpl.replace("__DATA__", data).replace("__IMGDATA__", imgdata).replace("__JUMPDATA__", jumpdata)
Path("web").mkdir(exist_ok=True)
Path("web/MorphoProfile.html").write_text(html, encoding="utf-8")
Path("index.html").write_text(html, encoding="utf-8")      # repo-root entry for Vercel / static hosting
print(f"wrote index.html + web/MorphoProfile.html ({len(html)//1024} KB)")
