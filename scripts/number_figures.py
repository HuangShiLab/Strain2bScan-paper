#!/usr/bin/env python3
"""Produce the draft-backed numbered figure set.

Outputs are Figures 1-11 and S1-S5 in ``figures/numbered`` and match
``figures/numbered/MANIFEST.tsv``. Single-source figures copy the vector PDF;
multi-source figures are montaged vertically from their source PNGs. Re-run after
regenerating any panel figure.
"""
import os, shutil
from pathlib import Path
from PIL import Image

PAPER = Path(__file__).resolve().parent.parent
FIG = PAPER / "figures/source_panels"
OUT = PAPER / "figures/numbered"; os.makedirs(OUT, exist_ok=True)

MAIN = {
    "01": ["overview"],
    "02": ["mash_2brad_vs_16s"],
    "03": ["fig06_native_2brad_mocks"],
    "04": ["saliva_fig07"],
    "05": ["sim2025_summary"],
    "06": ["enzyme_sweep"],
    "07": ["sim_headtohead"],
    "08": ["computational_performance"],
    "09": ["fig11_wms_toolcompare"],
    "10": ["saliva_fig08"],
}

SUPP = {
    # (former S1 rank-rank scatter is now Fig 2 panel B — see scripts/plot_fig2_combined.py)
    "S1": ["mock_msa1002_titration"],
    "S2": ["gate_calibration"],
    "S3": ["figS_tree_expansion"],    # DB-expansion cost: 20- vs 28-species combined tree
    "S4": ["figS_instrain_derep"],   # non-dereplicated vs dereplicated inStrain reference
    "S5": ["species_expansion"],     # exploratory StrainScan-curated-database comparison
    "S6": ["containment_mechanism", "refqual_figure"], # reference-completeness / containment optimization
    "S7": ["rankrank_2brad_16s"],
}

# Human-readable output labels. By default the source stems are joined with '+';
# these overrides avoid redundant names such as Fig12_fig11_wms_toolcompare.
LABELS = {
    "03": "native_2brad_mocks",
    "04": "real_saliva_timecourse",
    "05": "sim2025_summary",
    "06": "enzyme_ladder",
    "07": "sim_headtohead",
    "08": "computational_performance",
    "09": "wms_toolcompare",
    "10": "paired_saliva_wms_reanalysis",
    "S6": "reference_quality_containment",
    "S3": "tree_expansion",
    "S4": "instrain_derep",
    "S5": "species_expansion_strainscan",
    "S7": "rankrank_2brad_16s",
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
