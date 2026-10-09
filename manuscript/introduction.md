# Introduction

Microbial strains of the same species can differ sharply in phenotype; virulence, antibiotic
resistance, metabolic capacity and host interaction; so resolving *which* strains are present, and at
what abundance, is central to microbiome and clinical metagenomics. The 16S rRNA gene, the workhorse of
community surveys, cannot see this variation: its between-strain distances are essentially uncorrelated
with genome-wide divergence, so 16S resolves species but not strains (below, and Fig 2). Shotgun
metagenomes do carry strain information, and a family of tools now recovers it; reference k-mer methods
such as StrainScan [1] and StrainGE [2], marker-gene methods such as StrainPhlAn [3], and sketch-based
abundance estimators such as sylph [4]; making strain-level profiling routine for individual samples and
species.

Two obstacles keep shotgun-based strain profiling from two increasingly important settings.

**Scale.** Modern studies span hundreds of samples and aim to resolve strains across the
dozens-to-hundreds of species in a community. Full k-mer methods index and query the entire k-mer
content of every reference genome; per-sample cost is dominated by counting the sample's k-mers, and is
paid again for every species database queried. For *S* species across *N* samples this scales as
*N × S × (k-mer count + search)*; hours-to-days of compute and hundreds of megabytes to gigabytes of
memory; a practical ceiling for population-scale, multi-species strain surveillance.

**Low biomass and host contamination.** Many host-rich, low-input specimen types; saliva and other
oral sites, tumour and FFPE tissue, skin; yield little microbial DNA against a large human background
(often ≥90–99 % host). Shotgun sequencing spends most of its reads on the host, so the informative
bacterial fraction, and with it the rarer strains, is lost; strain resolution collapses exactly where it
would be most valuable.

sparse, reproducible set of fixed-length tags (roughly 1–2 % of the genome) [5]. Applied as a wet-lab
protocol (2bRAD-M / Fast2bRAD-M) [6,7], it enriches these informative markers *before* host DNA can swamp
the library, giving accurate profiles from picogram inputs, heavily host-contaminated samples and degraded
material. To date, however, 2bRAD has mainly been used for **species**-level profiling. Yet strain-resolution
methods do not need the whole genome; they score samples on *unique* markers that distinguish a strain or
a cluster of near-identical strains; and a 2bRAD tag set is exactly this kind of low-redundancy,
taxonomically informative marker set, 50–100× smaller than the full k-mer set.

We present **Strain2bScan**, which ports the StrainScan [1] resolution framework; within-species clustering
into a search structure, followed by unique-marker scoring and abundance estimation; onto 2bRAD tags, in
Rust. Its tag lengths and recognition patterns match Fast2bRAD-M / `2bRADExtraction.pl` [7]
exactly, so the tags are interoperable with the Fast2bRAD-M species layer, and Strain2bScan accepts **two input modes** that map onto the two obstacles above:

1. **Native BcgI 2bRAD experimental libraries**, enabling strain-level analysis of
   low-biomass and high-host microbiomes. In a controlled native BcgI mock, the reduction happens at the
   bench and Strain2bScan retained precision 1.0 and full strain recall at 99% host DNA; the primary
   comparator test under the same host ladder used in-silico-digested shotgun reads. On real saliva, it
   resolved individual-specific, within-day strain signatures and produced candidate low-abundance calls
   that require independent validation.

2. **In-silico digestion of conventional shotgun metagenomes**, enabling community-scale strain
   profiling. The sample is digested **once** and matched against every per-species database, so
   per-sample cost is independent of the number of species and linear in the number of samples: ~8×
   faster and ~11× lighter per sample than StrainScan, and projected to be 121–146× lower in cost on a 55-species community,
   while matching StrainScan's precision, its 0.5× detection onset and its recall on its own databases;
   and completing near-clonal *Mycobacterium tuberculosis* in ~1 s where StrainScan does not finish.

The two modes are directionally concordant: in three usable paired samples, 65/65 shotgun strain-cluster calls were contained in native calls. Thus a single
tool, on one 2bRAD tag space, spans both host-rich, low-input profiling and cohort-scale shotgun analysis.
