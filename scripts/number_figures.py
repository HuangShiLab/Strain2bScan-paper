#!/usr/bin/env python3
"""Produce a numbered figure set (figures/numbered/Fig01.png … Fig10.png, FigS1–S3) that matches the
manuscript numbering in manuscript/figures.md. Single-source figures are copied; multi-source figures
are montaged vertically (one source figure per row). Re-run after any panel figure is regenerated."""
import os, shutil, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg

FIG = os.path.expanduser("~/Downloads/Strain2bScan-paper/figures")
OUT = f"{FIG}/numbered"; os.makedirs(OUT, exist_ok=True)

MAIN = {
    "01": ["overview"],
    "02": ["mash_2brad_vs_16s"],
    "03": ["cross_species", "depth_sensitivity"],
    "04": ["containment_mechanism", "refqual_figure"],
    "05": ["enzyme_sweep"],
    "06": ["mock_native2brad"],
    "07": ["saliva_individual_discrimination", "saliva_temporal_ml"],
    "08": ["saliva_concordance"],
    "09": ["performance", "scalability", "community_throughput"],
    "10": ["species_expansion"],
    "11": ["sim_headtohead"],
    "12": ["mock_shotgun"],
}
SUPP = {
    # (former S1 rank-rank scatter is now Fig 2 panel B — see scripts/plot_fig2_combined.py)
    "S1": ["mock_msa1002_titration"],
    "S2": ["gate_calibration"],
    "S3": ["mock_hostcontam"],   # was Fig 6 (62-species DNA mock, host-contamination series)
}

def montage(tag, sources):
    imgs = [mpimg.imread(f"{FIG}/{s}.png") for s in sources if os.path.exists(f"{FIG}/{s}.png")]
    if not imgs:
        print(f"  Fig {tag}: MISSING sources {sources}"); return
    # normalize each row to a common width; stack vertically
    W = max(im.shape[1] for im in imgs)
    heights = [im.shape[0] * W / im.shape[1] for im in imgs]
    total_h = sum(heights)
    fig = plt.figure(figsize=(W / 150, total_h / 150), dpi=150)
    gs = fig.add_gridspec(len(imgs), 1, height_ratios=heights, hspace=0.02)
    for i, im in enumerate(imgs):
        ax = fig.add_subplot(gs[i, 0]); ax.imshow(im); ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    base = f"{OUT}/Fig{tag}_" + "+".join(sources)
    fig.savefig(base + ".png", dpi=150, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    # PDF: copy the original vector PDF for single-panel figures; montage-PDF for multi-panel
    if len(sources) == 1 and os.path.exists(f"{FIG}/{sources[0]}.pdf"):
        shutil.copy(f"{FIG}/{sources[0]}.pdf", base + ".pdf")
    else:
        fig2 = plt.figure(figsize=(W / 150, sum(heights) / 150), dpi=150)
        gs2 = fig2.add_gridspec(len(imgs), 1, height_ratios=heights, hspace=0.02)
        for i, im in enumerate(imgs):
            ax = fig2.add_subplot(gs2[i, 0]); ax.imshow(im); ax.axis("off")
        fig2.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig2.savefig(base + ".pdf", bbox_inches="tight", pad_inches=0.02); plt.close(fig2)
    print(f"  Fig {tag}: {len(imgs)} panel(s) -> {os.path.basename(base)}.png/.pdf")

print("Main figures:")
for tag, srcs in MAIN.items():
    montage(tag, srcs)
print("Supplementary:")
for tag, srcs in SUPP.items():
    montage(tag, srcs)
print(f"\nnumbered figures -> {OUT}")
