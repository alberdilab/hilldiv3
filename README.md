# hilldiv3

<!-- badges: start -->
[![R-CMD-check](https://github.com/alberdilab/hilldiv3/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alberdilab/hilldiv3/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/alberdilab/hilldiv3/graph/badge.svg)](https://app.codecov.io/gh/alberdilab/hilldiv3)
[![lint](https://github.com/alberdilab/hilldiv3/actions/workflows/lint.yaml/badge.svg)](https://github.com/alberdilab/hilldiv3/actions/workflows/lint.yaml)
[![CRAN version](https://www.r-pkg.org/badges/version/hilldiv3)](https://CRAN.R-project.org/package=hilldiv3)
[![CRAN checks](https://badges.cranchecks.info/summary/hilldiv3.svg)](https://cran.r-project.org/web/checks/check_results_hilldiv3.html)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/hilldiv3)](https://CRAN.R-project.org/package=hilldiv3)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`hilldiv3` measures and compares the diversity of biological communities —
OTU, ASV or MAG count tables — using **Hill numbers**. Hill numbers are a
single, intuitive family of diversity metrics: each one is an *effective number
of taxa* ("how many equally-abundant taxa would give this much diversity"), and
one parameter, the *diversity order* `q`, slides smoothly between counting all
taxa equally (richness), weighting them by abundance (Shannon) and focusing on
the common ones (Simpson). Because everything is expressed in the same currency,
results are directly comparable across samples, studies and methods.

From that one foundation, `hilldiv3` provides a unified toolkit for **neutral**
diversity (abundances only), **phylogenetic** diversity (accounting for how
related taxa are) and **functional** diversity (accounting for how different
their traits are), covering measurement, partitioning, (dis)similarity,
profiles, evenness and redundancy. You call the same functions for all three —
the diversity type is chosen by whether you supply a tree or a distance matrix.

`hilldiv3` is aimed at ecologists and microbiologists analysing DNA-based
community data (metabarcoding, amplicon and shotgun metagenomics), and more
generally anyone who needs comparable alpha, beta and gamma diversity across
taxonomic, phylogenetic and functional dimensions — including nested
multi-scale designs such as individuals within sites within regions.

## Installation

Install the released version from CRAN:

```r
install.packages("hilldiv3")
```

To install the version in this repository, which may contain changes not yet
available on CRAN:

```r
# install.packages("devtools")
devtools::install_github("alberdilab/hilldiv3")
```

## Quick start

```r
library(hilldiv3)

# Bundled simulated gut-microbiome MAG data.
hilldiv(gut_counts)                    # neutral Hill numbers q = 0, 1, 2
hilldiv(gut_counts, tree = gut_tree)   # neutral + phylogenetic

dist <- traits2dist(gut_traits)
hilldiv(gut_counts, dist = dist)                   # neutral + functional
hilldiv(gut_counts, tree = gut_tree, dist = dist)  # all three types at once

# Attribute Hill power-sum weight to a MAG and a set of MAGs.
sets <- list(C = rownames(gut_counts)[3],
             ABC = rownames(gut_counts)[1:3])
hillshare(gut_counts, q = c(0, 1, 2), tree = gut_tree, sets = sets)
hillshare(gut_counts, q = c(0, 1, 2), dist = dist, sets = sets)

# Rank taxa by their contribution to pairwise or collective dissimilarity.
head(hillcontrib(gut_counts, q = 1, metric = "C"))
head(hillcontrib(gut_counts, q = 1, metric = "C", by = "collective"))

# Core diversity results are tidy by default and plot directly.
plot(hillprof(gut_counts))             # diversity profile
hilldiv(gut_counts, out = "matrix")    # matrix: samples x q orders
```

## Documentation

Full documentation lives on the package website:
**<https://alberdilab.github.io/hilldiv3/>**

* **Get started** — `vignette("hilldiv3")`, a gentle introduction for anyone
  using Hill numbers for the first time.
* **Articles** — step-by-step guides to diversity types, taxon and set shares,
  partitioning & (dis)similarity, profiles/evenness/redundancy, and preparing
  your data.
* **Examples** — complete worked analyses (bat diets, gut microbiomes).
* **Reference** — every exported function, grouped by task.

## What's new in v3

If you have used [hilldiv2](https://github.com/anttonalberdi/hilldiv2), v3 is a
full redesign that keeps the familiar function names (`hilldiv()`,
`hillpart()`, `hilldiss()`, `hillsim()`, `hillpair()`, `hillred()`, `tss()`,
`traits2dist()`) while changing how they work underneath:

* A tested, isolated **compute engine** — the diversity maths lives in one
  place and is unit-tested independently of the user-facing functions.
* A single **validation/alignment layer** that reorders data to match the tree
  or distance matrix (fixing silent misalignment in v2), plus a real
  `match_data()` helper.
* **Broad input support**: matrices, data frames, tibbles, `phyloseq` and
  `TreeSummarizedExperiment` objects.
* Faster phylogenetic computation using an `ape` post-order traversal in place
  of `geiger::tips()`; `hillpair()` computes the shared structure once and
  reuses it across all sample pairs.
* **Tidy by default**: the core diversity and comparison functions return a
  long-format `data.frame` with `print()`/`plot()`/`autoplot()` methods; pass
  `out = "matrix"` for a plain matrix. `hillshare()` and `hillcontrib()` return
  data frames of taxon and set shares or contributions.
* An explicit `type = c("auto", "neutral", "phylogenetic", "functional")`
  argument that asserts and validates the diversity type (auto-detected by
  default).
* New functions `hillprof()` (diversity profiles), `hilleven()` (evenness)
  `hillshare()` (taxon and set shares) and `hillcontrib()` (taxon turnover
  contributions), hierarchical multi-scale
  partitioning in `hillpart()`, plus bundled example data (`gut_counts`,
  `gut_tree`, `gut_traits`).

See [`NEWS.md`](NEWS.md) for the full changelog.

## Testing

The diversity engine is covered by a `testthat` suite, including golden-value
tests against independent references (`vegan`, hand-computed values and
analytical identities such as Faith's PD and exact telescoping of nested
partitions). Run it with:

```r
devtools::test()
```

`R CMD check` runs on Linux, macOS and Windows on every push via GitHub
Actions.

## Getting help and contributing

* Questions and bug reports: open an
  [issue](https://github.com/alberdilab/hilldiv3/issues).
* Contributions are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) for how
  to report problems, propose features and submit pull requests, and for the
  support we aim to provide.

Please note that the hilldiv3 project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

## Citation

If you use `hilldiv3`, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button), and
the methodological guide:

* Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
  numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*, 19, 804-817.
  <https://doi.org/10.1111/1755-0998.13014>

## References

* Hill, M.O. (1973). Diversity and evenness. *Ecology*, 54, 427-432.
* Jost, L. (2007). Partitioning diversity into independent alpha and beta
  components. *Ecology*, 88, 2427-2439.
* Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures based
  on Hill numbers. *Phil. Trans. R. Soc. B*, 365, 3599-3609.
* Chiu, C.-H., Jost, L. & Chao, A. (2014). Phylogenetic beta diversity,
  similarity, and differentiation measures based on Hill numbers. *Ecological
  Monographs*, 84, 21-44.
* Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
  numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*, 19, 804-817.
