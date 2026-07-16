Bacteria of a single named species are not interchangeable. Isolates that a 16S survey — or even a
species-level shotgun profile — would report as one taxon can differ by tens of thousands of single-
nucleotide variants, by the presence or absence of whole genomic islands, plasmids and prophages, and,
consequently, in phenotypes that matter directly to health: virulence, antibiotic-resistance carriage,
toxin production, metabolic capacity and the ability to colonise a particular host. *Escherichia coli*
ranges from harmless commensal to enterohaemorrhagic pathogen; *Cutibacterium acnes* phylotypes partition
with healthy versus acne-associated skin; *Klebsiella pneumoniae* lineages differ sharply in
carbapenem resistance and hypervirulence. Resolving *which* strains are present in a sample, and at what
relative abundance, is therefore a central problem for microbiome science and for clinical metagenomics —
and it is a fundamentally different problem from cataloguing which species are present.

### The limits of 16S and the promise of shotgun strain typing

The 16S rRNA gene, the historical workhorse of community surveys, is structurally unable to resolve this
variation. The gene is short, highly conserved, and often present in multiple, non-identical copies per
genome; its between-strain distances are essentially uncorrelated with genome-wide divergence. This
chapter quantifies that failure directly (Fig 2): across fifteen species, the rank correlation between
16S distance and whole-genome distance has a median of only 0.36, with several species indistinguishable
from zero. 16S resolves species, not strains.

Whole-genome shotgun metagenomics does carry strain information, and over the past decade a family of
tools has learned to extract it. They fall into a few families. *Reference k-mer methods* — StrainScan,
StrainGE — index the k-mer content of a curated set of reference genomes per species and assign a sample's
k-mers to the best-supported combination of reference strains (or clusters of near-identical strains).
*Marker-gene methods* — StrainPhlAn — call SNPs within a fixed set of clade-specific marker genes and
reconstruct per-sample haplotypes. *Sketch/containment estimators* — sylph, and the Mash/sourmash
lineage they build on — estimate the containment of reference genomes in a sample from subsampled k-mer
sketches, trading per-SNP resolution for speed. These tools have made strain-level profiling routine for
individual samples and individual species, and each embodies a different point on the accuracy/speed/
generality trade-off surface.

### Two settings where shotgun strain typing struggles

Two obstacles keep shotgun-based strain profiling out of two settings that are becoming central to the
field.

**Scale.** Modern microbiome studies span hundreds to thousands of samples and increasingly aim to
resolve strains across the dozens-to-hundreds of species that co-occur in a community, not one species at
a time. Full k-mer methods index and query the entire k-mer content of every reference genome; the
dominant per-sample cost is counting the k-mers of the sample, and — critically — that cost is *paid
again for every species database queried*, because there is no shared per-sample representation across
species. For *S* species across *N* samples the work scales as *N × S ×* (k-mer count + search): the
species dimension multiplies rather than amortises. In practice this means hours-to-days of compute and
hundreds of megabytes to gigabytes of resident memory for a community-scale, multi-species survey — a
ceiling that this chapter shows is not intrinsic to the strain-typing problem but to the full-k-mer
representation.

**Low biomass and host contamination.** Many of the clinical niches where strain resolution would be most
valuable yield very little microbial DNA against an overwhelming human background. Saliva and other oral
sites, tumour and formalin-fixed paraffin-embedded (FFPE) tissue, and skin routinely present ≥90–99 %
host DNA. Ordinary shotgun sequencing spends its reads in proportion to DNA abundance, so at 99 % host a
sequencing run devotes ~99 % of its reads to the human genome; the informative bacterial fraction, and
with it the rarer strains, is simply not sequenced deeply enough to be resolved. Strain resolution
collapses in exactly the samples where it would matter most, and no amount of downstream computation can
recover reads that were never generated.

### Reduced-representation 2bRAD sequencing

Reduced-representation sequencing offers a route around both obstacles. Type-IIB restriction enzymes
(the basis of the "2bRAD" method) cut on *both* sides of a short, degenerate recognition site, excising a
fixed-length fragment — the 2bRAD tag, here 32–38 bp — at every occurrence of the site in a genome. The
result is a sparse, reproducible, genome-wide sample of roughly 1–2 % of the genome. Because the tag set
is defined by the recognition sequence rather than by abundance, the *same* loci are recovered from any
genome that contains them, making tags directly comparable across samples and reference genomes.

Two properties make this attractive for the two problem settings above. First, the tag set is 50–100×
smaller than the full k-mer set of a genome while remaining genome-wide and taxonomically structured —
precisely the low-redundancy, informative marker set that a strain-resolution framework needs, since such
frameworks never use the whole genome but only the *unique* markers that distinguish a strain or a cluster
of near-identical strains. Second, and decisively for the low-biomass problem, 2bRAD can be realised as a
*wet-lab* protocol (2bRAD-M and its faster variant Fast2bRAD-M): the reduction to informative tags happens
at the bench, during library preparation, *before* host DNA can swamp the sample. A native 2bRAD library
of a 99 %-host sample concentrates sequencing on bacterial tags rather than on the human genome, giving
accurate profiles from picogram inputs and heavily contaminated or degraded material.

To date, however, 2bRAD has been used only for **species**-level profiling. The tags carry strain-level
signal — this chapter demonstrates that they track whole-genome divergence where 16S does not — but no
method had yet exploited that signal for strain resolution, and it was not obvious *a priori* that a marker
set 50–100× sparser than a full k-mer index would retain enough discriminating, single-copy loci to
separate near-identical strains at realistic sequencing depths.

### Contribution of this chapter

This chapter presents **Strain2bScan**, a tool that ports the two-layer StrainScan resolution framework —
within-species clustering into a search structure, followed by unique-marker detection and abundance
estimation — onto 2bRAD tags, implemented in dependency-free Rust for speed and portability. Its tag
lengths and recognition patterns match Fast2bRAD-M / `2bRADExtraction.pl` exactly, so its markers are
interoperable with the Fast2bRAD-M species layer and the two operate in one shared tag space. Uniquely,
Strain2bScan accepts **two input modes**, one for each of the obstacles above:

1. **Native BcgI 2bRAD experimental libraries**, whose reads *are* the tags — enabling, for the first
   time, strain-level analysis of low-biomass, high-host microbiomes. Because the reduction happens at the
   bench, Strain2bScan holds precision 1.0 and full strain recall at 99 % host DNA, where in-silico
   digestion of ordinary shotgun of the same material loses most of its strains; on real saliva it
   resolves individual-specific, temporally stable strain signatures and recovers low-abundance strains
   host-limited shotgun cannot reach.

2. **In-silico digestion of conventional shotgun metagenomes** — enabling community-scale strain
   profiling. The sample is digested **once** into a shared tag representation and matched against every
   per-species database, so the marginal cost of an additional species is a hash-set lookup rather than a
   re-count, and per-sample cost is essentially independent of the number of species. This turns the
   *N × S ×* (count + search) scaling of full-k-mer methods into *N ×* (digest + *S·ε*), an *S*-fold
   structural saving that this chapter measures directly.

The two modes are shown to agree — in-silico and native digestion of the same material recover the same
strains — so a single tool, operating on one 2bRAD tag space, spans both the low-biomass clinical regime
and the cohort-scale regime. The chapter first establishes the shared foundation (the 2bRAD-versus-16S
motivation, core accuracy, sensitivity, and the effect of reference-genome quality), then develops the two
input-mode pillars in turn, and closes with a systematic head-to-head against StrainScan on a common
simulated benchmark that isolates the accuracy and the cost of the two approaches under identical inputs.