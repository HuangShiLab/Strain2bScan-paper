#!/usr/bin/env python3
"""Download the 62-species MSA-1002 mock DB genome panel.
55 species with pinned GCF accessions (multispecies_55x4.tsv) + 7 extra mock species
(fetched via ENA taxonomy). NCBI datasets primary, ENA-GCA fallback. Idempotent."""
import urllib.request, urllib.parse, json, os, gzip, io, zipfile, sys, time

PAPER = "/Users/macstudio/Downloads/Strain2bScan-paper"
OUT   = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad/genomes"
EXTRA = ["Bifidobacterium adolescentis", "Deinococcus radiodurans", "Phocaeicola vulgatus",
         "Lactobacillus gasseri", "Bacillus cereus", "Clostridium beijerinckii", "Cereibacter sphaeroides"]
UA = {"User-Agent": "s2b-bench"}

def get(url, t=180):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=t).read()

def dl_ncbi(acc, out):
    """NCBI datasets zip -> extract the .fna."""
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
    data = get(url, 240)
    z = zipfile.ZipFile(io.BytesIO(data))
    names = [n for n in z.namelist() if n.endswith(".fna")]
    if not names: return False
    open(out, "wb").write(z.read(names[0]))
    return os.path.getsize(out) > 0

def dl_ena(acc, out):
    """ENA GCA fallback (GCF->GCA, unversioned)."""
    gca = acc.replace("GCF", "GCA").split(".")[0]
    data = get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{gca}?download=true&gzip=true", 240)
    if data[:2] == b"\x1f\x8b": data = gzip.decompress(data)
    if not data.startswith(b">"): return False
    open(out, "wb").write(data)
    return os.path.getsize(out) > 0

def fetch(acc, out):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "cached"
    for fn, tag in ((dl_ncbi, "ncbi"), (dl_ena, "ena")):
        for attempt in range(3):
            try:
                if fn(acc, out): return tag
            except Exception as e:
                time.sleep(2)
                last = e
    return None

def ena_taxid(sp):
    d = json.loads(get(f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{urllib.parse.quote(sp)}", 30).decode())
    return d[0]["taxId"]

def ena_assemblies(tid, n=4):
    url = (f"https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&"
           f"query={urllib.parse.quote('tax_tree(%s)' % tid)}&fields=accession,assembly_level&format=tsv&limit=100")
    rows = [l.split("\t") for l in get(url, 60).decode().splitlines()[1:] if l.strip()]
    order = {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}
    rows.sort(key=lambda r: order.get(r[1].lower() if len(r) > 1 else "", 4))
    return [r[0] for r in rows[:n]]

# ---- 55 pinned species ----
pinned = {}
for line in open(f"{PAPER}/data/accessions/multispecies_55x4.tsv"):
    if line.startswith("#") or not line.strip(): continue
    sp, acc = line.rstrip("\n").split("\t")
    pinned.setdefault(sp, []).append(acc)

total_ok = 0; total_try = 0
for sp, accs in pinned.items():
    d = f"{OUT}/{sp.replace(' ', '_')}"; os.makedirs(d, exist_ok=True)
    ok = 0
    for a in accs:
        total_try += 1
        r = fetch(a, f"{d}/{a}.fna")
        if r: ok += 1; total_ok += 1
        else: print(f"  FAIL {sp} {a}", flush=True)
    print(f"[pinned] {sp}: {ok}/{len(accs)}", flush=True)

# ---- 7 extra mock species (ENA taxonomy discovery) ----
for sp in EXTRA:
    d = f"{OUT}/{sp.replace(' ', '_')}"; os.makedirs(d, exist_ok=True)
    try:
        accs = ena_assemblies(ena_taxid(sp), 4)
    except Exception as e:
        print(f"[extra]  {sp}: taxonomy FAIL {e}", flush=True); continue
    ok = 0
    for a in accs:
        total_try += 1
        # extra species come as GCA already -> ENA route works directly
        out = f"{d}/{a}.fna"
        if os.path.exists(out) and os.path.getsize(out) > 0: ok += 1; total_ok += 1; continue
        try:
            data = get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{a}?download=true&gzip=true", 240)
            if data[:2] == b"\x1f\x8b": data = gzip.decompress(data)
            if data.startswith(b">"):
                open(out, "wb").write(data); ok += 1; total_ok += 1
            else:
                print(f"  FAIL {sp} {a} (no fasta)", flush=True)
        except Exception as e:
            print(f"  FAIL {sp} {a} {e}", flush=True)
    print(f"[extra]  {sp}: {ok}/{len(accs)} ({','.join(accs)})", flush=True)

print(f"DONE  {total_ok}/{total_try} genomes -> {OUT}", flush=True)
