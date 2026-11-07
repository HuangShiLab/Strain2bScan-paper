# Discussion

Strain2bScan couples a sparse 2bRAD marker space to StrainScan-style [1] strain-cluster resolution in a
compact Rust implementation. The resulting capability is not another full-k-mer profiler with a faster
engine. It is a single marker framework in which native reduction precedes host depletion for low-biomass
libraries and one in-silico digest serves every species panel for cohort-scale shotgun analysis. Across the
benchmarks, 2bRAD markers preserved genome-wide strain relationships, precision matched StrainScan in
completed paired simulations, and the strongest operational gains appeared in database construction,
multi-species runtime and retention at 99% host DNA (Figs 2, 5–9).

Native 2bRAD-M addresses the setting in which shotgun coverage is most constrained. Because reduction occurs
before host DNA dominates the library, the primary MSA-1002 comparison retained both detection and abundance
at 99% human DNA, whereas the primary shotgun comparators lost one or both quantities in their own call
spaces (Fig 9). The saliva analysis extended this observation from a synthetic host ladder to person-specific
profiles: strain features identified the subject in every leave-one-timepoint-out test, compared with 78.1%
for species features, although the subject-aware improvement interval included zero (Fig 4). The full-depth
32-pair comparison showed both modality-specific sensitivity and modality-specific omissions, reinforcing that
these real-data comparisons describe assay scope rather than validate either modality (Fig 10).

A rhesus macaque GI spatial cohort extended this logic from person-specific stability to spatial ecology.
A cohort-specific panel built from 315 study MAGs recovered strain calls in 113/128 GI samples and 56/66
luminal samples, whereas only two body controls produced calls. Eighty-nine clusters were LI-enriched and
thirteen remained luminal-enriched after matched-position Bonferroni correction (Supplementary Tables
S5–S10). In shotgun/HUMAnN/bile-acid subsets, strain composition was associated with bile-acid classes,
butyrate and propionate programs and N-acetyl-sugar degradation. These are cohort-level ecological
associations rather than assignments of individual pathway genes to strains, but they show that
reduced-representation strain calls can support spatial and metabolic hypotheses in real animal cohorts.
Dominant-cluster sharing followed the same ecological logic: it was higher within monkeys between GI
positions than between monkeys, and old/young sharing was higher for LI luminal than LI mucosal calls
(Fig 11). This cluster-sharing summary is a population-resolution ecological measure rather than
SNV/popANI-level transmission proof.

In shotgun mode, the structural advantage is a change in search economics. A sample is digested once and the
resulting tags index against all species, so adding species increases panel searches rather than sample
k-mer counting. This design produced the 121–146× projected cohort-level advantage on the 55-species
analysis and made the marginal cost of community profiling approximately flat in species count (Fig 8).
The head-to-head timing comparison measures the combined effect of the sparse marker space, Rust
implementation and compact data structures; it was not designed to apportion these contributions causally.
Clean-mock precision depended on input composition and threshold semantics, so the scale claim should be
read together with the accuracy and cluster-space conditions in Tables 1–3. The framework's strongest
current position is therefore scalable strain resolution across mixed cohorts, rather than uniform
higher precision on every mock.

The biological unit of the claim is the strain cluster. At a 0.95 similarity cut, diverse species can
resolve close to individual genomes, whereas near-clonal populations resolve to fewer clusters. This is an
inherent property of short-read evidence shared by StrainScan [1] and Strain2bScan, not a loss introduced by
reduced-representation markers. Precision remained high in low-diversity species, while near-clonal recall
defined the most important algorithmic boundary (Fig 7). A within-cluster overlap or regression layer that
uses read placement after occurrence-based uniqueness would directly test whether near-clonal recall can
increase without sacrificing precision.

Panel and reference quality control the transportability of all strain calls. A niche-appropriate,
genome-rich oral panel recovered individual-specific structure that a generic pathogen panel did not, and
body-site-complete pangenomes are needed before sparse vaginal or meconium marker evidence can support
biological absence calls (Tables 4–7). The containment mode addresses incomplete genomes that have a more
complete relative, but markers absent from every reference remain unavailable. Broader pangenome panels,
completeness-aware gating and explicit absent-strain models are therefore the next tests for open-world
surveillance rather than optional engineering refinements.

The comparator scope is also explicit. StrainScan [1] and inStrain [8] were executed on the shotgun
benchmarks, whereas StrainGST [2] was analyzed in StrainGE's 0.90-reference space and reported separately
because those calls cannot be pooled with 0.95-space Strain2bScan or StrainScan metrics (Tables 8–9). The
other comparators likewise used their own call units, so cross-tool accuracy panels are descriptive rather
than pooled estimates. This
StrainGST analysis covered the primary simulations and ATCC shotgun mocks, but native BcgI input for StrainGE, StrainGR
and the supplementary curated-database subset remain untested. The comparison therefore establishes cost and
feasibility boundaries for the tested StrainGST workflow without claiming equivalence to the full StrainGE
family.

Two applications follow directly from the evidence. Native 2bRAD-M is appropriate when host content, input
mass or degradation limits shotgun recovery, with tumour, FFPE and skin as priority targets that remain to
be tested. In-silico digestion is appropriate for large shotgun cohorts in which repeated per-species
k-mer searches dominate cost. The rhesus GI analysis illustrates the second use: one panel spans multiple
species and links strain-cluster structure to regional, mucosal/luminal and metabolic contrasts. Resolving
those regimes in one marker space turns reduced-representation sequencing from a species-profile protocol
into a practical route to strain-resolved microbiome analysis.
