# Strain2bScan — Introduction (draft)

Microbial strains of the same species can differ sharply in phenotype — virulence,
antibiotic resistance, metabolic capacity and host interaction — so resolving *which* strains
are present, and at what abundance, is central to microbiome and clinical metagenomics.
Shotgun metagenomes carry this information, and a family of tools now recover it: reference
k-mer methods such as StrainScan and StrainGE, marker-gene methods such as StrainPhlAn, and
sketch-based abundance estimators such as sylph. These methods have made strain-level
profiling routine for individual samples and individual species.

The bottleneck is now **scale**. Modern studies span hundreds of samples and aim to resolve
strains across the dozens–to–hundreds of species that make up a community. Full k-mer methods
index and query the entire k-mer content of every reference genome; the per-sample cost is
dominated by counting the sample's k-mers, and it is paid again for every species database
queried. For a community of *S* species profiled across *N* samples this scales as *N × S ×
(k-mer count + search)*, which becomes hours-to-days of compute and hundreds of megabytes of
memory per run — a practical ceiling for population-scale, multi-species strain surveillance.

Reduced-representation sequencing offers a way out. Type-IIB restriction ("2bRAD") digestion
releases a sparse, reproducible set of fixed-length tags — roughly 1–2 % of the genome — that
has been used for accurate species-level profiling (2bRAD-M). Crucially, strain-resolution
methods like StrainScan do not need the whole genome: they score samples on *unique* markers
that distinguish a strain or a cluster of near-identical strains. A 2bRAD tag set is exactly
this kind of low-redundancy, taxonomically informative marker set — 50–100× smaller than the
full k-mer set — yet, to our knowledge, 2bRAD has not been used for strain-level profiling.

We present **Strain2bScan**, which ports the StrainScan resolution framework — within-species
clustering into a search tree, followed by unique-marker scoring and abundance estimation —
onto 2bRAD tags, implemented in dependency-free Rust. Strain2bScan accepts both BcgI 2bRAD
experimental libraries (a capability no full k-mer tool provides) and in-silico digestion of
conventional shotgun reads with up to 16 enzymes. Operating on the sparse tag set makes the
database ~50–100× smaller and lets a sample be digested **once** and matched against every
per-species database, so the per-sample cost becomes independent of the number of species.
We show that this yields ~8× faster, ~11× lighter per-sample profiling and near-flat scaling
in species count on real *Cutibacterium acnes* data — reaching ~132× throughput over a
per-species k-mer tool on a 55-species community — while holding **precision 1.0** across
species, matching StrainScan's detection onset at 0.5× coverage and matching or exceeding its
recall on StrainScan's own databases; and that the number of enzymes is a tunable knob trading
efficiency against resolution, with single-enzyme BcgI operation enabling native 2bRAD data.
Where recall is limited it reflects genuine near-clonality — the intrinsic short-read resolution
limit — which Strain2bScan reports honestly, and at which StrainScan failed to complete on
*Mycobacterium tuberculosis* while Strain2bScan finished in ~1 s.
