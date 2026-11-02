#!/usr/bin/env python3
"""Build paper/figures/src/fig1.svg (hilldiv3 package structure). Units: mm."""
import os
from xml.sax.saxutils import escape
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src", "fig1.svg")
W, H = 180, 90
NEU, PHY, FUN = "#66C2A5", "#FC8D62", "#8DA0CB"
INK, MID, LIGHT, EDGE = "#222222", "#555555", "#F0F0F0", "#8A8A8A"
FS = 2.2  # mm = 6.2 pt
o = []
def rect(x, y, w, h, fill="#fff", stroke=EDGE, sw=0.25, dash=None, r=0.8):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')
def text(x, y, s, size=FS, weight="normal", fill=INK, anchor="start", style="normal"):
    o.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" font-style="{style}" fill="{fill}" text-anchor="{anchor}">{escape(s)}</text>')
def lines(x, y, items, anchor="start", lh=2.75):
    for i, it in enumerate(items):
        s, wt, st = (it + ("normal", "normal"))[:3] if isinstance(it, tuple) else (it, "normal", "normal")
        text(x, y + i * lh, s, weight=wt, style=st, anchor=anchor)
def arrow(x1, y1, x2, y2, color=MID, sw=0.4):
    o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2-1.0}" y2="{y2}" stroke="{color}" stroke-width="{sw}"/>')
    o.append(f'<path d="M{x2} {y2} L{x2-1.5} {y2-0.9} L{x2-1.5} {y2+0.9} Z" fill="{color}"/>')

# column headers
for x, s in [(2, "1  Inputs"), (52, "2  Adapter and compute engine"), (136, "3  User-facing wrappers")]:
    text(x, 5.2, s, size=2.5, weight="bold")

# 1 inputs
ins = [
    (NEU, [("Count table", "bold"), "matrix, data frame, tibble,", "phyloseq, TreeSummarizedExperiment"]),
    (PHY, [("tree (optional)", "bold"), "phylo object; adds phylogenetic", "diversity"]),
    (FUN, [("dist (optional) + threshold τ", "bold"), "trait distances (see traits2dist());", "adds functional diversity"]),
    ("#9A9A9A", [("hierarchy (optional)", "bold"), "formula, e.g. ~ region / site", "(nested partition)"]),
]
ys = [10, 27, 44, 61]
for (c, ls), y in zip(ins, ys):
    rect(2, y, 44, 14.5, fill="#fff")
    o.append(f'<rect x="2" y="{y}" width="2.4" height="14.5" rx="0.8" fill="{c}"/>')
    lines(6.2, y + 4.3, [l if isinstance(l, tuple) else (l,) for l in ls])
    arrow(46, y + 7.2, 52, y + 7.2, color=c if c != "#9A9A9A" else MID, sw=0.6)

# adapter
rect(52, 10, 31, 79, fill=LIGHT)
text(54, 15, "Input adapter", weight="bold")
text(54, 17.8, "and validation", weight="bold")
lines(54, 23.5, ["Converts every input to", "one internal object.", "",
                 "Taxon names in counts,", "tree and dist are", "compared as sets;",
                 "counts are reordered", "to match (a mismatch", "is an error).", "",
                 "match_data() applies the", "same matching by hand."])
# engine
rect(88, 10, 44, 79, fill=LIGHT)
text(90, 15, "Compute engine (internal)", weight="bold")
text(90, 20.2, "Diversity type inferred from inputs", style="italic")
chips = [(NEU, "counts only: neutral"), (PHY, "+ tree: phylogenetic"), (FUN, "+ dist, τ: functional")]
for i, (c, s) in enumerate(chips):
    y = 22.2 + i * 6
    rect(90, y, 40, 5, fill=c, stroke="none")
    text(91.5, y + 3.4, s)
rect(90, 41.5, 40, 11, fill="#fff")
lines(91.5, 45.4, [("Shared structure, once", "bold"), "branch abundances: one", "post-order traversal"])
text(90, 57.5, "Engine pieces", style="italic")
pieces = ["per-sample alpha", "alpha / gamma / beta partition", "nested partition",
          "beta to (dis)similarity"]
for i, s in enumerate(pieces):
    y = 59.5 + i * 7
    rect(90, y, 40, 5.6, fill="#fff")
    text(91.5, y + 3.7, s)
# colour thread: thin bars along the engine pieces
for c, dx in [(NEU, 0), (PHY, 1), (FUN, 2)]:
    pass
# flow adapter -> engine -> wrappers
arrow(83, 49.5, 88, 49.5, sw=0.6)
arrow(132, 49.5, 136, 49.5, sw=0.6)

# 3 wrappers
wr = [
    ("hilldiv()", "per-sample Hill numbers"),
    ("hillprof(), hilleven()", "profiles across q; evenness"),
    ("hillpart()", "alpha/beta/gamma, flat or nested"),
    ("hilldiss(), hillsim()", "bounded (dis)similarity"),
    ("hillpair()", "pairwise dissimilarity"),
    ("hillred()", "redundancy"),
    ("hillshare(), hillcontrib()", "share and dissimilarity by taxon"),
]
for i, (f, d) in enumerate(wr):
    y = 10 + i * 8.2
    rect(136, y, 42, 6.6)
    text(137.5, y + 2.7, f, weight="bold")
    text(137.5, y + 5.4, d, fill=MID)
rect(2, 78, 44, 11, dash="1.2 0.8")
text(4, 82.2, "Helpers (optional)", weight="bold")
text(4, 85.0, "traits2dist(): trait distances", fill=MID)
text(4, 87.6, "tss(): total-sum scaling", fill=MID)
# outputs
rect(136, 70, 42, 19, fill=LIGHT)
text(137.5, 73.6, "Outputs", weight="bold")
lines(137.5, 76.6, ["tidy long data frame", "  (print and plot methods)", "dist object (hillpair)",
                    'matrix (out = "matrix")'])
o.append(f'<line x1="157" y1="65.8" x2="157" y2="68.6" stroke="{MID}" stroke-width="0.6"/>')
o.append(f'<path d="M157 70 L156.1 68.5 L157.9 68.5 Z" fill="{MID}"/>')

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm" viewBox="0 0 {W} {H}" '
       f'font-family="Helvetica, Arial, sans-serif">\n<rect width="{W}" height="{H}" fill="#fff"/>\n' + "\n".join(o) + "\n</svg>\n")
open(OUT, "w").write(svg)
print("wrote", OUT)
