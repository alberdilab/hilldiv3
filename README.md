# hilldiv3

<!-- badges: start -->
[![R-CMD-check](https://github.com/alberdilab/hilldiv3/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alberdilab/hilldiv3/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`hilldiv3` is a redesign of [hilldiv2](https://github.com/anttonalberdi/hilldiv2)
for analysing the diversity of biological communities (OTU/ASV/MAG tables)
based on Hill numbers. It provides a unified framework for **neutral**,
**phylogenetic** and **functional** diversity: measurement, partitioning,
(dis)similarity, profiles, evenness and redundancy.

## What's new in v3

* A tested, isolated **compute engine** — the diversity maths lives in one
  place and is unit-tested independently of the user-facing functions.
* A single **validation/alignment layer** that reorders data to match the tree
  or distance matrix (fixing silent-misalignment in v2), plus a real
  `match_data()` helper.
* **Dual input support**: matrices, data frames, tibbles, `phyloseq` and
  `TreeSummarizedExperiment` objects via `as_hill_input()`.
* Faster phylogenetic computation using an `ape` post-order traversal in place
  of `geiger::tips()`.
* New: `hillprof()` (diversity profiles) and `hilleven()` (evenness).
* The familiar `hilldiv()`, `hillpart()`, `hilldiss()`, `hillsim()`,
  `hillpair()`, `hillred()`, `tss()` and `traits2dist()` names are kept.

## Installation

```r
# install.packages("devtools")
devtools::install_github("alberdilab/hilldiv3")
```

## Documentation

Full documentation lives on the package website:
**<https://alberdilab.github.io/hilldiv3/>**

* **Get started** — `vignette("hilldiv3")`.
* **Articles** — in-depth guides to diversity types, partitioning &
  (dis)similarity, profiles/evenness/redundancy, and preparing your data.
* **Reference** — every exported function, grouped by task.

This package is a complete redesign of hilldiv2; see [`NEWS.md`](NEWS.md) for
what changed.

## Quick start

```r
library(hilldiv3)

counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))

hilldiv(counts)              # neutral Hill numbers q = 0, 1, 2
hilldiv(counts, tree = tree) # phylogenetic
hilldiv(counts, dist = dist) # functional
```

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
