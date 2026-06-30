import urllib.request, json, io, zipfile, os, time, sys
BASE="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT="multispecies/genomes"; os.makedirs(OUT, exist_ok=True)
K=4
SPECIES=["Staphylococcus aureus","Staphylococcus epidermidis","Escherichia coli","Klebsiella pneumoniae",
"Pseudomonas aeruginosa","Enterococcus faecalis","Enterococcus faecium","Streptococcus pneumoniae",
"Streptococcus pyogenes","Streptococcus agalactiae","Listeria monocytogenes","Salmonella enterica",
"Acinetobacter baumannii","Clostridioides difficile","Cutibacterium acnes","Helicobacter pylori",
"Campylobacter jejuni","Neisseria meningitidis","Haemophilus influenzae","Bacteroides fragilis",
"Bifidobacterium longum","Lacticaseibacillus rhamnosus","Vibrio cholerae","Proteus mirabilis",
"Serratia marcescens","Enterobacter cloacae","Citrobacter freundii","Streptococcus mutans",
"Bifidobacterium breve","Fusobacterium nucleatum","Porphyromonas gingivalis","Klebsiella oxytoca",
"Moraxella catarrhalis","Stenotrophomonas maltophilia","Mycobacteroides abscessus","Bordetella pertussis",
"Legionella pneumophila","Shigella sonnei","Yersinia enterocolitica","Pseudomonas putida"]
def get(url):
    req=urllib.request.Request(url, headers={"User-Agent":"strain2bscan-bench"})
    return urllib.request.urlopen(req, timeout=60).read()
def accessions(sp):
    url=f"{BASE}/genome/taxon/{urllib.parse.quote(sp)}/dataset_report?filters.assembly_level=complete_genome&filters.exclude_atypical=true&page_size=20"
    d=json.loads(get(url)); accs=[]
    for r in d.get("reports",[]):
        a=r.get("accession","")
        if a.startswith("GCF_"): accs.append(a)
        if len(accs)>=K: break
    return accs
def download(acc, d):
    url=f"{BASE}/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
    z=zipfile.ZipFile(io.BytesIO(get(url)))
    fna=[n for n in z.namelist() if n.endswith("_genomic.fna")]
    if not fna: return False
    with z.open(fna[0]) as fh, open(os.path.join(d,f"{acc}.fna"),"wb") as o: o.write(fh.read())
    return True
import urllib.parse
tot=0
for sp in SPECIES:
    safe=sp.replace(" ","_"); d=os.path.join(OUT,safe); os.makedirs(d,exist_ok=True)
    have=len([f for f in os.listdir(d) if f.endswith(".fna")])
    if have>=2: tot+=have; print(f"{sp}: have {have}, skip",flush=True); continue
    try: accs=accessions(sp)
    except Exception as e: print(f"{sp}: report FAIL {e}",flush=True); continue
    ok=0
    for a in accs:
        if os.path.exists(os.path.join(d,f"{a}.fna")): ok+=1; continue
        try:
            if download(a,d): ok+=1
        except Exception as e: print(f"  {a} dl FAIL {e}",flush=True)
        time.sleep(0.3)
    tot+=ok; print(f"{sp}: {ok} genomes",flush=True)
print(f"TOTAL genomes: {tot} across {len(SPECIES)} species")
