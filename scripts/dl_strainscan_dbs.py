#!/usr/bin/env python3
"""Download 3 of StrainScan's own pre-built species databases (Akkermansia muciniphila,
Prevotella copri, M. tuberculosis) from Google Drive, into ss_db3/<short>_db/ in the workdir.
Needs `gdown` (pip install gdown) and a cloned StrainScan repo (for StrainScan.py + the
jellyfish-linux fix -- see docs/headtohead_strainscan.md) to actually profile against them.
"""
import os, subprocess, tarfile

DBS = [
    ("akk", "1BAoi5u4JuPTapULjbZRgaBBW5JrNYd7P", "DB_Akk"),
    ("prev", "1qhc17ZSRop0hp5lrM2sAQiffv0PXuQ6r", "DB_Pre"),
    ("mtb", "18pPGyHODRwYV_d-l-xOKiLkbj7WcHCR_", "DB_Mtb"),
]

os.makedirs("ss_db3", exist_ok=True)
for short, gdrive_id, db_name in DBS:
    tgz = f"ss_db3/{short}_db.tar.gz"
    outdir = f"ss_db3/{short}_db"
    if os.path.isdir(f"{outdir}/{db_name}"):
        print(f"{short}: already extracted, skip")
        continue
    if not os.path.exists(tgz):
        subprocess.run(["gdown", gdrive_id, "-O", tgz], check=True)
    os.makedirs(outdir, exist_ok=True)
    with tarfile.open(tgz) as t:
        t.extractall(outdir)
    print(f"{short}: extracted -> {outdir}/{db_name}")
print("done. Point StrainScan.py -d at ss_db3/<short>_db/<DB_name>, e.g. ss_db3/akk_db/DB_Akk")
