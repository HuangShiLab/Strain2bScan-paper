#!/usr/bin/env python3
"""Produce the draft-backed numbered figure set.

Outputs are Figures 1-12 and S1-S4 in ``figures/numbered`` and match
``figures/numbered/MANIFEST.tsv``. Single-source figures copy the vector PDF;
multi-source figures are montaged vertically from their source PNGs. Re-run after
regenerating any panel figure.
"""
import os, shutil
from pathlib import Path
from PIL import Image

PAPER = Path(__file__).resolve().parent.parent
FIG = PAPER / "figures"
OUT = FIG / "numbered"; os.makedirs(OUT, exist_ok=True)

MAIN = {
    "01": ["overview"],
    "02": ["mash_2brad_vs_16s"],
    "03": ["cross_species", "depth_sensitivity"],
    "04": ["containment_mechanism", "refqual_figure"],
    "05": ["enzyme_sweep"],
    "06": ["fig6_2brad"],
    "07": ["saliva_individual_discrimination", "saliva_temporal_ml"],
    "08": ["saliva_concordance"],
    "09": ["performance", "scalability", "community_throughput"],
    "10": ["species_expansion"],
    "11": ["sim_headtohead"],
    "12": ["fig12_wms_toolcompare"],
}
SUPP = {
    # (former S1 rank-rank scatter is now Fig 2 panel B — see scripts/plot_fig2_combined.py)
    "S1": ["mock_msa1002_titration"],
    "S2": ["gate_calibration"],
    "S3": ["figS_tree_expansion"],    # DB-expansion cost: 20- vs 28-species combined tree
    "S4": ["figS_instrain_derep"],   # non-dereplicated vs dereplicated inStrain reference
}

# Human-readable output labels. By default the source stems are joined with '+';
# these overrides avoid redundant names such as Fig12_fig12_wms_toolcompare.
LABELS = {
    "06": "native_2brad_mocks",
    "10": "species_expansion_strainscan",
    "12": "wms_toolcompare",
    "S3": "tree_expansion",
    "S4": "instrain_derep",
}

def montage(tag, sources):
    missing = [s for s in sources if not os.path.exists(f"{FIG}/{s}.png")]
    if missing:
        raise FileNotFoundError(f"Fig {tag}: missing sources {missing}")
    imgs = [Image.open(f"{FIG}/{s}.png") for s in sources]
    # normalize each row to a common width, then stack vertically
    W = max(im.width for im in imgs)
    resized = []
    for im in imgs:
        if im.width != W:
            h = int(round(im.height * W / im.width))
            resized.append(im.resize((W, h), Image.LANCZOS))
        else:
            resized.append(im)
    total_h = sum(im.height for im in resized)
    canvas = Image.new("RGB", (W, total_h), (255, 255, 255))
    y = 0
    for im in resized:
        canvas.paste(im, (0, y))
        y += im.height
    label = LABELS.get(tag, "+".join(sources))
    base = f"{OUT}/Fig{tag}_{label}"
    canvas.save(base + ".png", "PNG")
    # PDF: keep the original vector PDF for single-panel figures; raster montage for multi-panel
    if len(sources) == 1 and os.path.exists(f"{FIG}/{sources[0]}.pdf"):
        shutil.copy(f"{FIG}/{sources[0]}.pdf", base + ".pdf")
    else:
        canvas.save(base + ".pdf", "PDF", resolution=150)
    print(f"  Fig {tag}: {len(sources)} panel(s) -> {os.path.basename(base)}.png/.pdf")

print("Main figures:")
for tag, srcs in MAIN.items():
    montage(tag, srcs)
print("Supplementary:")
for tag, srcs in SUPP.items():
    montage(tag, srcs)
print(f"\nnumbered figures -> {OUT}")
