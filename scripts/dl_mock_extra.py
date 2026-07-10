import urllib.request, urllib.parse, json, os, gzip, time
OUT="/private/tmp/claude-501/-Users-shihuang-Downloads-StrainScan-rust-claude/9f8a698a-0227-4dd2-8236-2536be092ed9/scratchpad/mock/extra_genomes"
SPECIES=["Bifidobacterium adolescentis","Deinococcus radiodurans","Phocaeicola vulgatus",
         "Lactobacillus gasseri","Bacillus cereus","Clostridium beijerinckii","Cereibacter sphaeroides"]
N=6
def get(url,t=120):
    return urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"s2b-bench"}),timeout=t).read()
def taxid(sp):
    d=json.loads(get(f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{urllib.parse.quote(sp)}",30).decode())
    return d[0]["taxId"]
def assemblies(tid,n):
    url=f"https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query={urllib.parse.quote('tax_tree(%s)'%tid)}&fields=accession,assembly_level&format=tsv&limit=80"
    rows=[l.split('\t') for l in get(url,60).decode().splitlines()[1:] if l.strip()]
    order={"complete genome":0,"chromosome":1,"scaffold":2,"contig":3}
    rows.sort(key=lambda r: order.get(r[1] if len(r)>1 else "",4))
    return [r[0] for r in rows[:n]]
def dl(gca,out):
    if os.path.exists(out) and os.path.getsize(out)>0: return True
    data=get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{gca}?download=true&gzip=true",180)
    if data[:2]==b'\x1f\x8b': data=gzip.decompress(data)
    if not data.startswith(b'>'): return False
    open(out,"wb").write(data); return True
for sp in SPECIES:
    d=f"{OUT}/{sp.replace(' ','_')}"; os.makedirs(d,exist_ok=True)
    try:
        accs=assemblies(taxid(sp),N)
    except Exception as e:
        print(f"{sp}: FAIL {e}",flush=True); continue
    ok=0
    for a in accs:
        try:
            if dl(a,f"{d}/{a}.fna"): ok+=1
        except Exception as e:
            print(f"  {a} FAIL {e}")
    print(f"{sp}: {ok}/{len(accs)} genomes",flush=True)
print("DONE")
