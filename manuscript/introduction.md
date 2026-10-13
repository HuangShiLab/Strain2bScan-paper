# Introduction

Microbial strains within a species can differ in virulence, antibiotic resistance, metabolic capacity and
host interaction. Resolving which strains are present, and at what abundance, is therefore central to
microbiome and clinical metagenomics. The 16S rRNA gene, the workhorse of community surveys, resolves
species but not strains because its between-strain distances are essentially uncorrelated with genome-wide
divergence (below, and Fig 2). Shotgun metagenomes retain this signal, and reference k-mer methods such as
StrainScan [1] and StrainGE [2], marker-gene methods such as StrainPhlAn [3], and sketch-based abundance
estimators such as sylph [4] now make strain-level profiling routine for individual samples and species.

Two constraints restrict broad use of shotgun-based strain profiling. The first is scale: full k-mer
methods index and query the k-mer content of each reference genome, and a sample's k-mers are counted again
for every species database. For *S* species across *N* samples, cost therefore grows approximately as
*N × S × (k-mer count + search)*, creating a practical barrier to population-scale, multi-species
surveillance. The second is host-rich, low-biomass material. Saliva and other oral sites, tumour and FFPE
tissue, and skin often contain 90–99% host DNA, so shotgun sequencing spends most reads on the host and
the informative bacterial fraction is depleted precisely where strain resolution would be most valuable.

Reduced-representation sequencing offers a route through both constraints. Type-IIB restriction ("2bRAD")
digestion releases a sparse, reproducible set of fixed-length tags that represent roughly 1–2% of the
genome [5]. As a wet-lab protocol (2bRAD-M / Fast2bRAD-M) [6,7], it enriches informative markers before
host DNA dominates the library and supports profiling of picogram inputs, heavily host-contaminated samples
and degraded material. This marker space has mainly supported species-level profiling, but strain-resolved
methods do not require the whole genome. They score samples on unique markers that distinguish a strain or
a cluster of near-identical strains, and a 2bRAD tag set is exactly such a low-redundancy, taxonomically
informative marker set, 50–100× smaller than the full k-mer set.

We present **Strain2bScan**, which ports the StrainScan [1] resolution framework of within-species
clustering, unique-marker scoring and abundance estimation onto 2bRAD tags in Rust. Its tag lengths and
recognition patterns match Fast2bRAD-M / `2bRADExtraction.pl` [7] exactly, so the species and strain layers
share one marker space. The implementation supports two entry points. Native BcgI 2bRAD experimental
libraries enable strain-level analysis of low-biomass, high-host microbiomes. In-silico digestion of
conventional shotgun metagenomes enables digest-once, cohort-scale analysis, in which the marginal cost of
another species is a search against a compact panel rather than a second pass over sample k-mers. We
evaluate marker information content, strain detection and abundance accuracy, comparator performance,
computational scaling and application to real oral, vaginal, meconium and infant-gut metagenomes.
