#!/usr/bin/env python3
"""Fetch assembly QC (level, CheckM completeness/contamination, contigs, length) for the downloaded
16S-panel genomes and classify each as high-quality (complete/near-complete) or not."""
import urllib.request, urllib.parse, urllib.error, json, os, glob, time

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
GEN = f"{WORK}/genomes_16s"
UA = {"User-Agent": "s2b"}

def get_retry(url, tries=6):
    for i in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90).read()
        except urllib.error.HTTPError as e:
            if e.code == 429: time.sleep(3*(i+1)); continue
            raise
    raise RuntimeError("429")

def report(sp):
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{urllib.parse.quote(sp)}/"
           f"dataset_report?filters.assembly_version=current&page_size=800")
    return json.loads(get_retry(url).decode()).get("reports", [])

rows = [["species","accession","assembly_level","checkm_completeness","checkm_contamination","n_contigs","total_len","high_quality"]]
for sd in sorted(glob.glob(f"{GEN}/*/")):
    sp = os.path.basename(sd.rstrip("/")); spname = sp.replace("_", " ")
    have = {os.path.basename(g)[:-4].split(".")[0] for g in glob.glob(f"{sd}/*.fna")}  # GCA numeric stems
    meta = {}
    for r in report(spname):
        for acc in (r.get("accession",""), r.get("paired_accession","")):
            num = acc.replace("GCF_","").replace("GCA_","").split(".")[0]
            if f"GCA_{num}" in have or num in {h.replace('GCA_','') for h in have}:
                ck = r.get("checkm_info",{}) or {}
                st = r.get("assembly_stats",{}) or {}
                meta[f"GCA_{num}"] = (
                    r.get("assembly_info",{}).get("assembly_level",""),
                    ck.get("completeness"), ck.get("contamination"),
                    st.get("number_of_contigs"), st.get("total_sequence_length"))
    for g in sorted(have):
        key = g if g.startswith("GCA_") else f"GCA_{g}"
        lvl, comp, cont, nc, tl = meta.get(key, ("",None,None,None,None))
        hq = (lvl in ("Complete Genome","Chromosome")
              and (comp is None or comp >= 97) and (cont is None or cont <= 5))
        rows.append([sp, key, lvl, comp, cont, nc, tl, int(bool(hq))])
    print(f"{sp}: {sum(1 for x in rows[1:] if x[0]==sp and x[7]==1)}/{len(have)} high-quality", flush=True)

open(f"{WORK}/genome_qc.tsv","w").write("\n".join("\t".join(str(c) for c in r) for r in rows)+"\n")
print(f"\nwrote {WORK}/genome_qc.tsv")
