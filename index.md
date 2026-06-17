# hilldiv3

`hilldiv3` is a redesign of
[hilldiv2](https://github.com/anttonalberdi/hilldiv2) for analysing the
diversity of biological communities (OTU/ASV/MAG tables) based on Hill
numbers. It provides a unified framework for **neutral**,
**phylogenetic** and **functional** diversity: measurement,
partitioning, (dis)similarity, profiles, evenness and redundancy.

## What’s new in v3

- A tested, isolated **compute engine** — the diversity maths lives in
  one place and is unit-tested independently of the user-facing
  functions.
- A single **validation/alignment layer** that reorders data to match
  the tree or distance matrix (fixing silent-misalignment in v2), plus a
  real
  [`match_data()`](https://alberdilab.github.io/hilldiv3/reference/match_data.md)
  helper.
- **Dual input support**: matrices, data frames, tibbles, `phyloseq` and
  `TreeSummarizedExperiment` objects via `as_hill_input()`.
- Faster phylogenetic computation using an `ape` post-order traversal in
  place of `geiger::tips()`;
  [`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md)
  computes the shared structure once and reuses it across all sample
  pairs.
- **Tidy by default**: every `hill*` function returns a long-format
  `data.frame` with
  [`print()`](https://rdrr.io/r/base/print.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/`autoplot()`
  methods; pass `out = "matrix"` for the legacy shape.
- An explicit
  `type = c("auto", "neutral", "phylogenetic", "functional")` argument
  that asserts and validates the diversity type (auto-detected by
  default).
- New:
  [`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md)
  (diversity profiles) and
  [`hilleven()`](https://alberdilab.github.io/hilldiv3/reference/hilleven.md)
  (evenness), plus bundled example data (`gut_counts`, `gut_tree`,
  `gut_traits`).
- The familiar
  [`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md),
  [`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md),
  [`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md),
  [`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md),
  [`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md),
  [`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md),
  [`tss()`](https://alberdilab.github.io/hilldiv3/reference/tss.md) and
  [`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md)
  names are kept.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("alberdilab/hilldiv3")
```

## Documentation

Full documentation lives on the package website:
**<https://alberdilab.github.io/hilldiv3/>**

- **Get started** —
  [`vignette("hilldiv3")`](https://alberdilab.github.io/hilldiv3/articles/hilldiv3.md).
- **Articles** — in-depth guides to diversity types, partitioning &
  (dis)similarity, profiles/evenness/redundancy, and preparing your
  data.
- **Reference** — every exported function, grouped by task.

This package is a complete redesign of hilldiv2; see
[`NEWS.md`](https://alberdilab.github.io/hilldiv3/NEWS.md) for what
changed.

## Quick start

``` r

library(hilldiv3)

# Bundled simulated gut-microbiome MAG data.
hilldiv(gut_counts)                    # neutral Hill numbers q = 0, 1, 2
hilldiv(gut_counts, tree = gut_tree)   # phylogenetic

dist <- traits2dist(gut_traits)
hilldiv(gut_counts, dist = dist)       # functional

# Results are tidy by default and plot directly.
plot(hillprof(gut_counts))             # diversity profile
hilldiv(gut_counts, out = "matrix")    # legacy matrix shape
```

## References

- Hill, M.O. (1973). Diversity and evenness. *Ecology*, 54, 427-432.
- Jost, L. (2007). Partitioning diversity into independent alpha and
  beta components. *Ecology*, 88, 2427-2439.
- Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity
  measures based on Hill numbers. *Phil. Trans. R. Soc. B*, 365,
  3599-3609.
- Chiu, C.-H., Jost, L. & Chao, A. (2014). Phylogenetic beta diversity,
  similarity, and differentiation measures based on Hill numbers.
  *Ecological Monographs*, 84, 21-44.
- Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of
  Hill numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*,
  19, 804-817.
