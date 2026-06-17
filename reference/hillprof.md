# Diversity profile across a range of orders

Compute a diversity profile: Hill numbers evaluated over a fine sweep of
diversity orders `q`. Profiles are the standard diagnostic for comparing
the diversity of assemblages, since the ranking of samples can change
with `q`.

## Usage

``` r
hillprof(
  data,
  q = seq(0, 3, by = 0.1),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  out = c("tibble", "matrix")
)
```

## Arguments

- data:

  Counts: a numeric vector (one sample), a matrix/data.frame (taxa x
  samples), a `phyloseq` object or a `TreeSummarizedExperiment`.

- q:

  Numeric vector of diversity orders to evaluate. Defaults to a fine
  sweep from 0 to 3.

- tree:

  A phylogenetic tree of class `phylo` whose tip labels match the taxa
  in `data`.

- dist:

  A functional distance matrix (or `dist`) over the taxa.

- tau:

  Optional functional distance threshold. Defaults to `max(dist)`.

- out:

  Output type: `"tibble"` (default, long format ready for plotting) or
  `"matrix"`.

## Value

A long-format `data.frame` of class `hill_profile` (columns `q`,
`sample`, `value`) with a
[plot()](https://alberdilab.github.io/hilldiv3/reference/plot.hill_profile.md)
method, or a matrix (orders in rows, samples in columns) when
`out = "matrix"`.

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
prof <- hillprof(counts)
#> Computing neutral diversity profile over 31 orders.
plot(prof)
```
