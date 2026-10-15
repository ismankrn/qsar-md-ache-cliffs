#!/usr/bin/env python
"""Compose PyMOL renders into publication-ready figures.

All figure text is in English: these panels go into the manuscript verbatim."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
F = ROOT / "figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

def trim(img, pad=6):
    """Trim excess white margin."""
    import numpy as np
    a = img[..., :3] if img.ndim == 3 else img
    mask = (a < 0.985).any(axis=2) if a.ndim == 3 else a < 0.985
    ys, xs = mask.any(axis=1), mask.any(axis=0)
    if not ys.any():
        return img
    y0, y1 = ys.argmax(), len(ys) - ys[::-1].argmax()
    x0, x1 = xs.argmax(), len(xs) - xs[::-1].argmax()
    y0, x0 = max(0, y0-pad), max(0, x0-pad)
    y1, x1 = min(img.shape[0], y1+pad), min(img.shape[1], x1+pad)
    return img[y0:y1, x0:x1]

# ---------------- Gambar 1 ----------------
a = trim(mpimg.imread(F/"_p1a.png"))
b = trim(mpimg.imread(F/"_p1b.png"))
fig = plt.figure(figsize=(11.0, 6.9), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.02,
                      left=0.012, right=0.988, top=0.935, bottom=0.215)
panels = [(a, "A", "Overall fold of hAChE"),
          (b, "B", "Architecture of the catalytic gorge")]
for i, (img, lab, ttl) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i]); ax.imshow(img); ax.axis("off")
    ax.set_anchor("N")                      # top-align both panels
    bb = ax.get_position()
    # titles placed in FIGURE coords so they line up despite differing aspects
    fig.text(bb.x0 + 0.004, 0.945, lab, fontsize=19, fontweight="bold",
             ha="left", va="top")

handles = [
    Line2D([], [], color="#1f6fd0", lw=6, label="Trp86 — catalytic anionic site (CAS)"),
    Line2D([], [], color="#b22222", lw=6, label="Trp286 — peripheral anionic site (PAS)"),
    Line2D([], [], color="#ff8c00", lw=6, label="Tyr337 — swinging gate"),
    Line2D([], [], color="#228b22", lw=6, label="Catalytic triad Ser203–His447–Glu334"),
    Line2D([], [], color="#e8e800", lw=6, label="Donepezil (crystal pose)"),
]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
           fontsize=10.5, bbox_to_anchor=(0.5, 0.055), handlelength=1.5,
           columnspacing=1.8, handletextpad=0.6)
fig.savefig(F/"Figure1_hAChE_structure.png", dpi=300, facecolor="white")
fig.savefig(F/"Figure1_hAChE_structure.pdf", facecolor="white")
plt.close(fig)

# ---------------- Figure 2 (one panel) ----------------
# Panels B and C showed huperzine A self-docked into 4EY5 against the same ligand
# cross-docked into 4EY7, which was the evidence that receptor rigidity and not the
# protocol limits the docking. That diagnosis was taken out of the manuscript, so the
# figure keeps only the redocking it validates. The renders for B and C are still in
# figures/_p2b.png and _p2c.png if it is ever wanted back.
img = trim(mpimg.imread(F/"_p2.png"))
fig = plt.figure(figsize=(6.2, 5.9), dpi=300)
ax = fig.add_axes([0.01, 0.20, 0.98, 0.76]); ax.imshow(img); ax.axis("off")
ax.text(0.97, 0.97, "RMSD = 0.75 \u00c5", transform=ax.transAxes, fontsize=12.5,
        ha="right", va="top", fontweight="bold", color="#1a7a1a",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#1a7a1a", lw=1.2))

handles = [
    Line2D([], [], color="#22c022", lw=6, label="Crystal pose (reference)"),
    Line2D([], [], color="#e020e0", lw=6, label="Docked pose (top-ranked)"),
    Patch(facecolor="0.78", edgecolor="0.5", label="Binding-cavity surface"),
    Line2D([], [], color="#1f6fd0", lw=6, label="Trp86 (CAS)"),
    Line2D([], [], color="#ff8c00", lw=6, label="Tyr337"),
]
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
           fontsize=9.5, bbox_to_anchor=(0.5, 0.01), handlelength=1.5,
           columnspacing=2.0, handletextpad=0.6)
fig.savefig(F/"Figure2_redocking_validation.png", dpi=300, facecolor="white")
fig.savefig(F/"Figure2_redocking_validation.pdf", facecolor="white")
plt.close(fig)
print("done:")
for f in sorted(F.glob("Figure*")):
    print(f"  {f.name}  {f.stat().st_size/1e6:.2f} MB")
