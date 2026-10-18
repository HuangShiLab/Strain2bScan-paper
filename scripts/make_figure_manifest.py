#!/usr/bin/env python3
"""Regenerate figures/numbered/MANIFEST.tsv from the current numbered figure set."""
import csv, hashlib
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "figures/numbered/MANIFEST.tsv"
SOURCES = {
 "Fig01_overview": "figures/source_panels/overview.{png,pdf}",
 "Fig02_mash_2brad_vs_16s": "figures/source_panels/mash_2brad_vs_16s.{png,pdf}",
 "Fig03_enzyme_sweep": "figures/source_panels/enzyme_sweep.{png,pdf}",
 "Fig04_native_2brad_mocks": "figures/source_panels/fig06_native_2brad_mocks.{png,pdf}",
 "Fig05_saliva_fig07": "figures/source_panels/saliva_fig07.{png,pdf} (five-panel A-E)",
 "Fig06_saliva_fig08": "figures/source_panels/saliva_fig08.{png,pdf} (directional paired comparison)",
 "Fig07_sim2025_summary": "figures/source_panels/sim2025_summary.{png,pdf} (complete 2,025-sample simulation aggregate)",
 "Fig08_sim_headtohead": "figures/source_panels/sim_headtohead.{png,pdf} (Strain2bScan/StrainScan/StrainGST 225-sample simulation)",
 "Fig09_performance+scalability+community_throughput": "figures/source_panels/performance.{png,pdf} + scalability.{png,pdf} + community_throughput.{png,pdf} (combined)",
 "Fig10_wms_toolcompare": "figures/source_panels/fig11_wms_toolcompare.{png,pdf} (retitled for manuscript numbering)",
 "FigS1_mock_msa1002_titration": "figures/source_panels/mock_msa1002_titration.{png,pdf}",
 "FigS2_gate_calibration": "figures/source_panels/gate_calibration.{png,pdf}",
 "FigS3_tree_expansion": "figures/source_panels/figS_tree_expansion.{png,pdf}",
 "FigS4_instrain_derep": "figures/source_panels/figS_instrain_derep.{png,pdf}",
 "FigS5_species_expansion_strainscan": "figures/source_panels/species_expansion.{png,pdf}",
 "FigS6_reference_quality_containment": "figures/source_panels/containment_mechanism.{png,pdf} + refqual_figure.{png,pdf} (combined)",
}
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for stem,source in SOURCES.items():
    png=ROOT/"figures/numbered"/f"{stem}.png"; pdf=ROOT/"figures/numbered"/f"{stem}.pdf"
    if not png.exists() or not pdf.exists(): raise FileNotFoundError(stem)
    rows.append({"figure_number":stem,"sources":source,"png":png.name,"pdf":pdf.name,
                 "png_sha256":sha(png),"pdf_sha256":sha(pdf)})
with path.open("w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=["figure_number","sources","png","pdf","png_sha256","pdf_sha256"],delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
print(f"wrote {len(rows)} rows to {path}")
