# Attribute Hill dissimilarity to taxa or sets

Show which taxa contribute most to dissimilarity between every pair of
samples or among samples in a group. Default singleton-taxon
contributions are nonnegative and add up to the corresponding value from
[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md)
or
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md).
A taxon's `share` is its fraction of that dissimilarity; it is `NA` when
the dissimilarity is zero.

## Usage

``` r
hillcontrib(
  data,
  q = c(0, 1, 2),
  metric = c("S", "C", "U", "V"),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  type = c("auto", "neutral", "phylogenetic", "functional"),
  by = c("pairwise", "collective"),
  groups = NULL,
  sets = NULL
)
```

## Arguments

- data:

  A count table (taxa x samples) or a supported object; a single sample
  is not meaningful for partitioning.

- q:

  Numeric vector of diversity orders (\>= 0). Defaults to `c(0, 1, 2)`
  (richness, Shannon, Simpson).

- metric:

  Dissimilarity metric(s) to return, any of `"S"`, `"C"`, `"U"`, `"V"`.
  Defaults to all four.

- tree:

  A phylogenetic tree of class `phylo` whose tip labels match the taxa
  in `data`.

- dist:

  A functional distance matrix (or `dist`) over the taxa.

- tau:

  Optional functional distance threshold. Defaults to `max(dist)`.

- type:

  Diversity type: `"auto"` (default) infers it from the inputs (counts
  only -\> neutral, `+tree` -\> phylogenetic, `+dist` -\> functional);
  an explicit `"neutral"`, `"phylogenetic"` or `"functional"` asserts
  the type and is validated against the inputs (e.g. `"phylogenetic"`
  requires a `tree`; `"neutral"` ignores any tree/dist carried by the
  object).

- by:

  `"pairwise"` (default) compares every pair of samples; `"collective"`
  compares all samples together or each specified group.

- groups:

  For `by = "collective"`, an optional named list of character vectors
  of sample names. Each group must contain at least two distinct
  samples. Defaults to one group named `"all"` containing every sample.

- sets:

  An optional named list of character vectors of taxon names, as in
  [`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md).
  Defaults to one set per taxon. Sets may overlap.

## Value

A data frame with columns `first`, `second`, `group`, `q`, `type`,
`metric`, `set`, `contribution`, `share`, and `rank`. Pairwise rows use
`first` and `second`; collective rows use `group`. Rows are ranked from
largest to smallest contribution within each comparison, order, and
metric. Sets that omit taxa sum to less than the total dissimilarity;
overlapping sets can sum to more.

## Details

For each taxon (or phylogenetic branch), the contribution kernel is the
Jensen gap between its pooled abundance and its abundances across
samples. The gap is zero if its abundance is constant. At `q = 0` the
gap measures uneven occurrence; at `q = 1` it is the entropy gap; at
other orders it is the power-sum gap. For neutral and functional
diversity the kernel belongs directly to a taxon. For phylogenetic
diversity each branch's gap is assigned to descendant taxa in proportion
to their own neutral Jensen gaps at the same order. This is one explicit
allocation of shared branches, rather than a unique taxon-level
decomposition of phylogenetic turnover.

The kernel fractions are multiplied by the selected dissimilarity
metric. Thus singleton-taxon contributions add exactly to that metric,
while their relative ranking is the same for `S`, `C`, `U`, and `V` at a
given order and comparison. Zero-length branches and taxa absent from
every selected sample receive zero. Every selected sample must have a
positive total. Functional beta uses the package's raw-count
partitioning, so unequal sample totals can contribute to functional
dissimilarity even when relative composition is unchanged.

## See also

[`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md),
[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md),
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md)

## Examples

``` r
counts <- matrix(c(5, 0, 2, 3, 4, 0, 0, 2, 6), nrow = 3,
                 dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
hillcontrib(counts, q = 1, metric = "C")
#>   first second group q    type metric set contribution      share rank
#> 1    s1     s2  <NA> 1 neutral      C   B   0.28571429 0.62848337    1
#> 2    s1     s2  <NA> 1 neutral      C   C   0.14285714 0.31424169    2
#> 3    s1     s2  <NA> 1 neutral      C   A   0.02603771 0.05727494    3
#> 4    s1     s3  <NA> 1 neutral      C   A   0.35714286 0.63781170    1
#> 5    s1     s3  <NA> 1 neutral      C   B   0.12500000 0.22323410    2
#> 6    s1     s3  <NA> 1 neutral      C   C   0.07780745 0.13895420    3
#> 7    s2     s3  <NA> 1 neutral      C   C   0.37500000 0.58972928    1
#> 8    s2     s3  <NA> 1 neutral      C   A   0.21428571 0.33698816    2
#> 9    s2     s3  <NA> 1 neutral      C   B   0.04659928 0.07328255    3
hillcontrib(counts, q = 1, metric = "C", by = "collective")
#>   first second group q    type metric set contribution     share rank
#> 1  <NA>   <NA>   all 1 neutral      C   C    0.1601445 0.3704045    1
#> 2  <NA>   <NA>   all 1 neutral      C   A    0.1515502 0.3505265    2
#> 3  <NA>   <NA>   all 1 neutral      C   B    0.1206555 0.2790690    3
```
