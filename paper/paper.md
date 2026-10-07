---
title: 'hilldiv3: a unified Hill numbers framework for neutral, phylogenetic and functional diversity analysis in R'
tags:
  - R
  - ecology
  - biodiversity
  - Hill numbers
  - phylogenetic diversity
  - functional diversity
  - microbiome
  - metabarcoding
authors:
  - name: Antton Alberdi
    orcid: 0000-0002-2875-6446
    corresponding: true
    affiliation: 1
affiliations:
  - name: Center for Evolutionary Hologenomics, Globe Institute, University of Copenhagen, Denmark
    index: 1
    ror: 035b05819
date: 1 October 2026
bibliography: paper.bib
---

# Summary

Hill numbers are a mathematically coherent family of diversity metrics that
express the diversity of a community as an *effective number* of equally
abundant, equally distinct units [@Hill1973; @Jost2006]. A single order
parameter, $q$, sets how strongly rare and common units are weighted, and the
same construction extends from counting taxa (neutral diversity) to weighting
them by evolutionary relatedness (phylogenetic diversity) or by trait
dissimilarity (functional diversity) [@Chao2014]. `hilldiv3` is an R package
that computes all of these within one engine: per-sample Hill numbers, diversity
profiles and evenness, multiplicative alpha/beta/gamma partitioning
[@Whittaker1960; @Jost2007],
(dis)similarity and pairwise distance matrices, and phylogenetic and functional
redundancy. Every function takes the same inputs and returns results on the
same scale whichever diversity type is analysed, so neutral, phylogenetic and
functional results can be compared directly. `hilldiv3` also provides nested
multiplicative partitioning across arbitrarily deep sampling hierarchies for all
three diversity types.

# Statement of need

High-throughput sequencing now routinely produces feature tables of operational
taxonomic units (OTUs), amplicon sequence variants (ASVs) or
metagenome-assembled genomes (MAGs), holding hundreds to thousands of features
across many samples, from which researchers compute a range of diversity
metrics [@Gotelli2001; @Ji2013; @Neu2021]. Classical indices are grounded in
different mathematical bases and sit on different, non-linear scales, so values
from one study often cannot be placed alongside those from another
[@Chao2014; @Alberdi2019]. Hill numbers resolve much of this, and
@Alberdi2019 argued for their routine use with DNA-based data. In practice,
however, a complete Hill-number analysis still means combining several
packages, each with its own input conventions, alignment assumptions, output
formats and normalisations. Phylogenetic and functional analyses of large
datasets can also take a long time with existing implementations.

`hilldiv3` is aimed at ecologists and microbiologists analysing metabarcoding,
amplicon and shotgun metagenomic data, and more generally at anyone who needs
comparable alpha, beta and gamma diversity across taxonomic, phylogenetic and
functional dimensions. One worked example from the package documentation shows
why this matters. In a simulated gut microbiome of 24 MAGs in 12 hosts, an
intervention swaps which of two deep clades dominates. The two clades are
functional mirrors of each other. Alpha diversity is almost identical between
control and treatment in all three diversity types. Between-group beta
diversity, by contrast, grows with $q$ for neutral and phylogenetic diversity
but stays at $\beta \approx 1$ for functional diversity (\autoref{fig:beta}):
community membership changes while function does not. Detecting this
decoupling requires all three diversity types on the same beta scale. A second
example, dietary metabarcoding of bats nested within roosts and habitats, uses
nested partitioning to show that neutral turnover is spread between individuals
and habitats, whereas phylogenetic turnover is concentrated between habitats.

![Between-group beta diversity across diversity orders for neutral,
phylogenetic and functional Hill numbers in the simulated gut-microbiome
example. The dashed line marks $\beta = 1$ (no turnover).\label{fig:beta}](figure1.png){ width=80% }

# State of the field

Several R packages compute Hill numbers or related measures. `vegan`
[@vegan] provides neutral Hill numbers through Rényi diversity and many
non-Hill dissimilarities. `BAT` [@Cardoso2015] covers taxonomic, phylogenetic
and functional alpha and beta diversity, but outside the Hill-number
partitioning framework. `iNEXT` [@Hsieh2016] focuses on rarefaction and
extrapolation, which is complementary to `hilldiv3` rather than overlapping
with it. `entropart` [@Marcon2015] partitions metacommunity diversity,
including similarity-based and phylogenetic diversity, but uses different
normalisations and requires ultrametric trees for phylogenetic diversity.
`hierDiversity` [@Marion2015] provides hierarchical partitioning of Hill
numbers, but not for phylogenetic or functional diversity. The closest
comparator, `hillR` [@Li2018], covers all three diversity types and partitions
them, but follows different conventions. For functional diversity it reports
attribute diversity $FD = {}^qD^2 Q$ at a fixed threshold $\tau = Q$ (Rao's
quadratic entropy), so its functional beta lies on a squared scale that is not
comparable with its neutral and phylogenetic betas. None of these packages
offers nested partitioning across more than two levels for phylogenetic or
functional diversity.

`hilldiv3` is a ground-up redesign of the author's `hilldiv` [@hilldiv] and
`hilldiv2` [@hilldiv2] packages, not a contribution to an existing package.
Nested partitioning and the common-scale guarantee depend on a single shared
construction of alpha, beta and gamma for all three diversity types. This
could not be added to packages with different normalisations without breaking
their existing results. The package documentation lists each convention
difference from `hillR` and `hilldiv2`, together with the arguments that
reproduce their results exactly where possible.

We benchmarked `hilldiv3` against `hilldiv2`, `hillR`, `entropart`, `vegan`
and `BAT` on simulated data (200 taxa, 50 samples, 10 iterations). The tests
covered four operations (alpha diversity, flat partitioning, nested
partitioning and pairwise beta matrices), each in neutral, phylogenetic and
functional form. `hilldiv3` supported all twelve combinations directly;
`hillR` and `hilldiv2` supported nine, lacking nested partitioning. For
neutral alpha diversity and neutral partitioning, other packages were faster,
because `hilldiv3`'s fixed per-call cost of validation and output assembly
dominates when the computation is trivial. With a phylogeny or trait distance
the ordering largely reversed. The largest differences were in pairwise
matrices: `hilldiv3` built them in 54–70 ms using 45–98 MB of allocated memory,
whereas `hilldiv2` needed 2.1–50.6 s and up to 20 GB, and `hillR` needed
0.8–16.2 s and up to 10.6 GB. The benchmark scripts and results ship with the
package.

# Software design

`hilldiv3` keeps the mathematics separate from the user interface. Per-sample
alpha, the alpha/gamma/beta partition, the nested partition and the
beta-to-(dis)similarity conversions all live in an internal compute engine. The
user-facing functions (\autoref{tab:functions}) are thin wrappers around it.
This makes the mathematics auditable in one place and keeps the public API
stable.

: Main user-facing functions. \label{tab:functions}

| Function | Purpose |
|---|---|
| `hilldiv()` | Per-sample Hill numbers |
| `hillprof()`, `hilleven()` | Diversity profiles across $q$; evenness |
| `hillpart()` | Alpha/beta/gamma partitioning, flat or nested |
| `hilldiss()`, `hillsim()` | Bounded (dis)similarity [@Chiu2014] |
| `hillpair()` | Pairwise dissimilarity as `dist` objects |
| `hillred()` | Phylogenetic or functional redundancy |
| `traits2dist()`, `tss()` | Trait distances [@Gower1971]; total-sum scaling |

**One input layer.** All inputs pass through one adapter and validation layer.
It accepts matrices, data frames, tibbles and `phyloseq` [@McMurdie2013] and
`TreeSummarizedExperiment` [@Huang2021] objects. It compares taxon names
between the count table and the tree or distance matrix as sets, and reorders
the counts to match. This replaces `hilldiv2`'s order-dependent comparison,
which could silently misalign data. Validation adds a small fixed cost to every
call; we accepted this cost, which is visible in the neutral-alpha benchmark,
in exchange for safety.

**The diversity type follows from the inputs.** Counts alone give neutral
diversity, adding a `tree` gives phylogenetic diversity [@Chao2010], and adding
a trait `dist` gives functional diversity at a distance threshold $\tau$
[@ChiuChao2014]. An optional `type` argument asserts and validates the
intended type. For phylogenetic diversity on non-ultrametric trees, the
`reference` argument makes explicit whether samples are read at one common tree
depth (the default) or each at its own. This is a convention on which existing
packages silently differ.

**Shared structure is computed once.** Branch abundances are obtained in a
single post-order traversal with `ape` [@Paradis2019], instead of a
descendant look-up for each tip. `hillpair()` computes the type-specific
structure once over all samples and reuses it for every pair. Pairs can
optionally be computed in parallel through `future`.

**Nested partitioning.** Supplying a formula such as `~ region / site` to
`hillpart()` returns the finest-scale alpha, one beta per level of the
hierarchy and the overall gamma [@Jost2007; @Chao2016; @Gaggiotti2018]. Each
sample receives a single global weight, and each beta is the ratio of
consecutive scale-level Hill numbers. Each beta is therefore independent of
alpha for all $q$ and bounded by the number of child units it aggregates, and
the chain telescopes exactly ($\gamma = \alpha \prod_k \beta_k$). The three
diversity types share this construction and differ only in the per-unit measure
and its normaliser. The phylogenetic version holds one tree depth fixed across
all levels, and the functional version holds one $\tau$ fixed. The current
release assumes equal sample weights, and the documentation states this.

**Testing.** A `testthat` suite tests the engine independently of the
wrappers. It includes golden-value tests against `vegan` and analytical
identities: Faith's PD as the $q = 0$ phylogenetic gamma, continuity at the
$q = 1$ limit, and exact telescoping for every diversity type, including on
non-ultrametric trees. `R CMD check` runs on Linux, macOS and Windows through
continuous integration. Results are returned as tidy long-format data frames
with `print` and `plot` methods, or as plain matrices with `out = "matrix"`.

# Research impact statement

`hilldiv3` is [available on CRAN](https://CRAN.R-project.org/package=hilldiv3),
where version 3.0.0 was first published on 6 October 2026. It continues a
lineage that is already in research use. Its predecessor `hilldiv` was
distributed through CRAN from October 2019 until it
was archived in March 2025, and was downloaded about 24,900 times from the
RStudio CRAN mirror. The accompanying
guide, @Alberdi2019, has been cited 273 times (OpenAlex, October 2026).
`hilldiv3` keeps `hilldiv`'s function names, so existing users can migrate
directly. Its documentation explains how to reproduce `hilldiv2` and `hillR`
results where the two agree.

**TODO (author):** add concrete, verifiable evidence of `hilldiv3`'s own use,
for example: studies or preprints that used it (with DOIs); groups or projects
that have adopted it (e.g. Earth Hologenome Initiative analyses); pipelines
that integrate it. JOSS does not accept
statements of intended future use as evidence.

# AI usage disclosure

Generative AI was used in developing `hilldiv3`. Claude Code (Anthropic;
Claude Opus 5 and Claude Opus 5.5) assisted with code generation and
refactoring of the compute engine and wrappers, test scaffolding,
documentation and vignette drafting, and benchmark scripts. It was also used to
adapt the author's manuscript into this paper's format. The author made the
core design decisions, including the unified engine architecture and the
global-weight construction for nested partitioning, and reviewed, edited and
validated all AI-assisted outputs.
The numerical results were checked against independent references through the
test suite.

# Acknowledgements

The author acknowledges the Danish National Research Foundation (DNRF), grant
DNRF143.

# References
