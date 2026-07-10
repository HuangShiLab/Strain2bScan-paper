# Strain2bScan — session handoff (→ Mac Studio, 2026-07-10)

Paste this as the opening message on Mac Studio to continue. The two git repos are the durable
state; scratchpad intermediates (genomes, mock reads, DBs) are ephemeral but reproducible from the
committed scripts + accession/manifest files.

## Repos (clone both on Mac Studio)
- **Software:** `github.com/HuangShiLab/Strain2bScan` (Rust profiler; formerly "Strain2bfunc").
- **Paper:** `github.com/HuangShiLab/Strain2bScan-paper` (manuscript, benchmarks, figures).
- Build: `cargo build --release` (zero deps). Binary: `target/release/strain2bscan`.

## What Strain2bScan is
Fast/light Rust strain profiler on **2bRAD-reduced markers** (32–38 bp type-IIB tags). Ports
StrainScan's clustering + unique-marker scoring onto 2bRAD; tags are **interoperable with
Fast2bRAD-M / 2bRADExtraction.pl**. Two input modes: native BcgI 2bRAD libraries, or in-silico
digestion of shotgun. Two-layer: Layer-1 species gate → Layer-2 within-species strain resolution.

## Software state (all committed + pushed)
- Latest: **Layer-1 gate = breadth-fraction term + three-tier reporting** (`species_tier` in
  `src/main.rs`). Flags: `--min-species-markers` (200), `--min-species-marker-frac` (0.0),
  `--min-species-detect` (10). 27 tests pass. Defaults validated by calibration (keep them).
- Three outcomes per species: strain-resolved / detected-not-resolvable / absent.
- Key files: `src/main.rs` (multi-profile), `src/enzymes.rs` (Fast2bRAD-M table),
  `src/identify.rs` (Layer-2), `src/markers.rs` (digest+canonical), `src/db.rs`, `src/cst.rs`.

## Paper state
- **Figure plan:** `docs/results_figure_plan.md` (10-figure sequence).
- **Done:** Fig 1 schematic (`figures/overview.*`); all simulated benchmarks (cross-species, depth,
  enzyme knob, robustness=panelsize+refqual, efficiency+55-species community scale, head-to-head vs
  StrainScan); **Fig 9 real-DNA mock** (`figures/mock_msa1002_titration.*`, `docs/mock_msa1002.md`);
  **2bRAD-vs-whole-genome Mash recompute** 3 species (`results/mash_2brad_vs_wgs_recomputed.tsv`,
  Spearman 0.89–0.996); gate calibration (`docs/gate_calibration.md`).
- **Prose:** `manuscript/{abstract,introduction,results,discussion,methods}.md` rewritten around the
  corrected story: **accurate (precision 1.0) + fast/light + community-scale + 2bRAD-native +
  Fast2bRAD-M-interoperable**. (The old "resolvability/operating-envelope" framing was a digestion
  bug, now fixed — don't reintroduce it.)
- ⚠ **Check `git status`/`git push`** first: the SRA-manifest commit may be unpushed (a GitHub SSL
  hiccup at handoff time). Push it.

## Environment quirks (carry these forward)
- **NCBI API/FTP was IP-blocked** on the old MacBook (fast SSL EOF). Test NCBI on Mac Studio; if
  blocked, use **ENA**: genomes `ebi.ac.uk/ena/browser/api/fasta/<GCA>?download=true&gzip=true`
  (our pins are GCF → ENA serves GCA, same numeric part, unversioned auto-resolves); SRA reads via
  the `fastq_ftp` field; metadata via ENA portal `filereport`. Figshare + ENA + GitHub all reachable.
- **No bioinformatics tools installed.** Install on Mac Studio as needed: `barrnap`+`mash` (16S
  figure), `inStrain` (WMS contrast), `unar` (rar, if `bsdtar` fails). `.rar` from Figshare extract
  fine with `bsdtar`.
- **HKU HPC** `hpc2021.hku.hk` (user `shihuang`): needs HKUVPN + key/ControlMaster (password-based
  SSH won't work; never handle the password). Not set up; optional.

## Data sources
- **Figshare 12272360 v8** (2B-RAD-M unique tags DB): MSA-1002 mock `2BRAD-M`(id 25933027, used for
  Fig 9), `shotgun_1/2`(25933009/25933021), `16S`(25933024), `FFPE_samples`(25933072). URL
  `ndownloader.figshare.com/files/<id>`.
- **SRA PRJNA1131785** (npj Biofilms Microbiomes 2025, s41522-025-00851-2): MSA-1002 mock **+90/99%
  human** and real HoC **saliva (diurnal 8×4) / oral cancer / ECC**; 2bRAD/WMS/16S/PacBio. Full
  manifest: `data/sra_PRJNA1131785_manifest.tsv`. **~1.78 TB total but the 2bRAD subset is only
  ~23 GB** (mock 4.2 GB; real HoC 18.9 GB); WMS is 1.75 TB (grab only a few for inStrain).
  Run→condition mapping via ENA sample `source_material_id` (e.g. saliva "S11-9") + the paper's
  supplementary table.

## Three prior Strain2bfunc results being folded in (PPTX `Strain2bFunc_20260710--Shi.pptx`)
1. **2bRAD-vs-16S Mash concordance** (motivation fig): 2bRAD tracks whole-genome (recomputed:
   Spearman 0.89–0.996 ✅); **16S half still TODO** (needs barrnap). `results/mash_2brad_vs_16s.tsv`
   has the prior 15-species table.
2. **DNA mock** (Fig 9): ✅ done on real MSA-1002 2bRAD (precision 1.0, recall 20/20 to 0.1 ng).
   Extensions TODO: 90/99%-host axis (SRA), WMS/inStrain contrast, cross-species abundance.
3. **Saliva individual discrimination** (Fig 10, TODO): target = strain Adonis R² up to 0.87 >
   species 0.72; ML 100% at strain level. Data = saliva 2bRAD in PRJNA1131785.

## Next steps on Mac Studio (priority)
1. **Saliva (Fig 10)** — download diurnal-saliva 2bRAD (subset of PRJNA1131785), build oral-species
   DBs, `multi-profile`, PERMANOVA R² strain vs species. Highest-value biological result.
2. **Mock host-contamination axis** — grab mock 90/99%-human 2bRAD (4.2 GB), extend Fig 9.
3. **16S contrast** — install barrnap+mash, extract 16S for the 15 species, compute WGS/2bRAD/16S
   distances → complete the motivation figure.
4. **WMS/inStrain contrast** — install inStrain, a few mock/saliva WMS runs.
5. **Oral cancer + ECC** — new application directions (2bRAD in PRJNA1131785; FFPE on Figshare).

## Memory (may not auto-transfer; key facts captured above + in committed docs)
`~/.claude/.../memory/`: user-huanglab, project-2brad-strainscan, project-strain2bfunc-prior-data,
reference-sra-prjna1131785, reference-ena-genome-fallback, feedback-no-claude-coauthor (omit the
Co-Authored-By Claude trailer in commits). Consider copying the memory dir, or just rely on this doc.

## To resume
Clone both repos, read `docs/results_figure_plan.md` + this file, then: "continue Strain2bScan on
Mac Studio — start with the saliva 2bRAD (Fig 10)" and point me at the data paths.
