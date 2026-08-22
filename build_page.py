#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Inject computed results into index_template.html -> index.html (repo root, for static hosting)
and web/MorphoProfile.html (local copy). Both are the same single self-contained file.

    python build_page.py            # write index.html + web/MorphoProfile.html
    python build_page.py --check    # verify index.html matches template + data, WITHOUT writing

--check exists so a reviewer can confirm the deployed page really is template+JSON without
mutating the working tree (a plain rebuild rewrites line endings on Windows and leaves git dirty).
Exit code 0 = in sync, 1 = out of date.
"""
import sys
from pathlib import Path

CHECK = "--check" in sys.argv[1:]

def esc(p, fallback=None):
    if not Path(p).exists():
        if fallback is None: raise SystemExit(f"missing required input: {p}")
        print(f"WARNING: {p} missing -> using empty placeholder")
        return fallback
    return Path(p).read_text(encoding="utf-8").replace("</", "<\\/")   # safe inside <script>

tpl  = Path("index_template.html").read_text(encoding="utf-8")
data = esc("web/webdata.json")
imgd = esc("images/out/webimages.json",
           '{"compounds":{},"source":"(not generated: run images/image_pipeline.py first)","channel_note":"","feature_names":[]}')
jumpd = esc("mvp/jump_mvp_results.json", "{}")
html = tpl.replace("__DATA__", data).replace("__IMGDATA__", imgd).replace("__JUMPDATA__", jumpd)

out = Path("index.html")
if CHECK:
    if not out.exists():
        print("index.html is missing"); raise SystemExit(1)
    cur = out.read_text(encoding="utf-8")
    if cur.replace("\r\n", "\n") == html.replace("\r\n", "\n"):
        print(f"index.html is in sync with index_template.html + the three JSON inputs ({len(html)//1024} KB)")
        raise SystemExit(0)
    print("index.html is OUT OF DATE — run: python build_page.py")
    raise SystemExit(1)

Path("web").mkdir(exist_ok=True)
Path("web/MorphoProfile.html").write_text(html, encoding="utf-8")
out.write_text(html, encoding="utf-8")          # repo-root entry for Vercel / static hosting
print(f"wrote index.html + web/MorphoProfile.html ({len(html)//1024} KB)")
