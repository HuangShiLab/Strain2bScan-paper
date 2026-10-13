# Discussion

Strain2bScan couples a sparse 2bRAD marker space to StrainScan-style [1] strain-cluster resolution in a
compact Rust implementation. The resulting capability is not another full-k-mer profiler with a faster
engine. It is a single marker framework in which native reduction precedes host depletion for low-biomass
libraries and one in-silico digest serves every species panel for cohort-scale shotgun analysis. Across the
benchmarks, 2bRAD markers preserved genome-wide strain relationships, precision matched StrainScan in
completed paired simulations, and the strongest operational gains appeared in database construction,
multi-species runtime and retention at 99% host DNA (Figs 2, 3, 9–11).

Native 2bRAD-M addresses the setting in which shotgun coverage is most constrained. Because reduction occurs
before host DNA dominates the library, the primary MSA-1002 comparison retained both detection and abundance
at 99% human DNA, whereas the tested shotgun profilers lost one or both quantities (Fig 11). The saliva
analysis extended this observation from synthetic host ladders to person-specific profiles: strain features
identified the subject in every leave-one-timepoint-out test, compared with 78.1% for species features
(Fig 7). This is a small within-subject classification result rather than independent strain-identity
validation. Its value lies in showing that native libraries carry enough strain-specific structure to
support subject-level discrimination, a property that 65/65 concordant shotgun calls and 128–163 additional
candidate strain-cluster calls per sample further support on usable paired libraries (Fig 8).

In shotgun mode, the structural advantage is a change in search economics. A sample is digested once and the
resulting tags index against all species, so adding species increases panel searches rather than sample
k-mer counting. This design produced the 121–146× projected cohort-level advantage on the 55-species
analysis and made the marginal cost of community profiling approximately flat in species count (Fig 9).
Clean-mock precision depended on input composition and threshold semantics, so the scale claim should be
read together with the accuracy and cluster-space conditions in Tables 1–3. The framework's strongest
current position is therefore scalable strain resolution across mixed cohorts, rather than uniform
precision superiority on every mock.

The biological unit of the claim is the strain cluster. At a 0.95 similarity cut, diverse species can
resolve close to individual genomes, whereas near-clonal populations resolve to fewer clusters. This is an
inherent property of short-read evidence shared by StrainScan [1] and Strain2bScan, not a loss introduced by
reduced-representation markers. Precision remained high in low-diversity species, while near-clonal recall
defined the most important algorithmic boundary (Fig 3). A within-cluster overlap or regression layer that
uses read placement after occurrence-based uniqueness would directly test whether near-clonal recall can
increase without sacrificing precision.

Panel and reference quality control the transportability of all strain calls. A niche-appropriate,
genome-rich oral panel recovered individual-specific structure that a generic pathogen panel did not, and
body-site-complete pangenomes are needed before sparse vaginal or meconium marker evidence can support
biological absence calls (Tables 4–7). The containment mode addresses incomplete genomes that have a more
complete relative, but markers absent from every reference remain unavailable. Broader pangenome panels,
completeness-aware gating and explicit absent-strain models are therefore the next tests for open-world
surveillance rather than optional engineering refinements.

The comparator scope is also explicit. StrainScan [1] and inStrain [8] were executed on the frozen shotgun
benchmarks, whereas StrainGST [2] was rerun in StrainGE's 0.90-reference space and reported separately
because those calls cannot be pooled with 0.95-space Strain2bScan or StrainScan metrics (Tables 8–9). This
rerun covered the primary simulations and ATCC shotgun mocks, but native BcgI input for StrainGE, StrainGR
and the supplementary curated-database subset remain untested. The comparison therefore establishes cost and
feasibility boundaries for the tested StrainGST workflow without claiming equivalence to the full StrainGE
family.

Two applications follow directly from the evidence. Native 2bRAD-M is appropriate when host content, input
mass or degradation limits shotgun recovery, with tumour, FFPE and skin as priority targets that remain to
be tested. In-silico digestion is appropriate for large shotgun cohorts in which repeated per-species
k-mer searches dominate cost. Resolving those two regimes in one marker space turns reduced-representation
sequencing from a species-profile protocol into a practical route to strain-resolved microbiome analysis.
