#!/usr/bin/env python3
"""Download ~50 genomes/species for the 15 slide-deck species (16S-vs-2bRAD motivation figure).
NCBI accession lists (RefSeq-preferred, complete/chromosome first) -> ENA FASTA (fast, unthrottled)."""
import urllib.request, urllib.parse, urllib.error, json, os, gzip, time
from concurrent.futures import ThreadPoolExecutor

OUT = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad/genomes_16s"
N_PER, WORKERS = 50, 8
UA = {"User-Agent": "s2b-bench"}
SPECIES = [
    "Fusobacterium nucleatum", "Porphyromonas gingivalis", "Lactiplantibacillus plantarum",
    "Klebsiella pneumoniae", "Escherichia coli", "Prevotella copri", "Salmonella enterica",
    "Akkermansia muciniphila", "Phocaeicola dorei", "Clostridioides difficile",
    "Mycobacterium tuberculosis", "Streptococcus pneumoniae", "Staphylococcus aureus",
    "Cutibacterium acnes", "Staphylococcus epidermidis",
]

def get(url, t=120):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=t).read()

def get_retry(url, t=90, tries=6):
    for i in range(tries):
        try: return get(url, t)
        except urllib.error.HTTPError as e:
            if e.code == 429: time.sleep(3 * (i + 1)); continue
            raise
    raise RuntimeError("429 backoff exhausted")

def list_accs(sp, n):
    enc = urllib.parse.quote(sp)
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{enc}/dataset_report"
           f"?filters.assembly_version=current&page_size=600")
    reports = json.loads(get_retry(url, 90).decode()).get("reports", [])
    order = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    def key(r):
        lvl = r.get("assembly_info", {}).get("assembly_level", "")
        refseq = 0 if r.get("source_database") == "SOURCE_DATABASE_REFSEQ" else 1
        return (order.get(lvl, 4), refseq)
    reports.sort(key=key)
    seen, uniq = set(), []
    for r in reports:
        a = r["accession"]
        if a not in seen: seen.add(a); uniq.append(a)
    return uniq[:n]

def dl_ena(gcf, d):
    gca = gcf.replace("GCF", "GCA").split(".")[0]
    out = f"{d}/{gca}.fna"
    if os.path.exists(out) and os.path.getsize(out) > 0: return 1
    try:
        data = get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{gca}?download=true&gzip=true", 120)
        if data[:2] == b"\x1f\x8b": data = gzip.decompress(data)
        if not data.startswith(b">"): return 0
        open(out, "wb").write(data); return 1
    except Exception:
        return 0

tot = 0
for sp in SPECIES:
    d = f"{OUT}/{sp.replace(' ', '_')}"; os.makedirs(d, exist_ok=True)
    have = len([x for x in os.listdir(d) if x.endswith(".fna")])
    if have >= N_PER: print(f"{sp}: {have} (cached)", flush=True); tot += have; continue
    try: accs = list_accs(sp, N_PER)
    except Exception as e: print(f"{sp}: LIST FAIL {e}", flush=True); continue
    ok = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for r in ex.map(lambda a: dl_ena(a, d), accs): ok += r
    tot += ok
    print(f"{sp}: {ok}/{len(accs)} genomes", flush=True)
print(f"DONE  {tot} genomes -> {OUT}", flush=True)
