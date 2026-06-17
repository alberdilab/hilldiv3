# Hill numbers-based dissimilarity

Compute overall (multi-sample) dissimilarity metrics from the
Hill-number beta diversity following Chiu et al. (2014). These are the
complements of the similarities returned by
[`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md).

## Usage

``` r
hilldiss(
  data,
  q = c(0, 1, 2),
  metric = c("S", "C", "U", "V"),
  tree = NULL,
  dist = NULL,
  tau = NULL
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

## Value

A matrix of dissimilarities (diversity orders in rows, metrics in
columns), or a vector if a single metric is requested.

## See also

[`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md),
[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md),
[`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
hilldiss(counts)
#> dissimilarity from neutral Hill numbers of "q0", "q1", and "q2".
#>            S         C         U         V
#> q0 0.3333333 0.2000000 0.3333333 0.2000000
#> q1 0.6081390 0.5229848 0.5229848 0.4369251
#> q2 0.7308320 0.7308320 0.5758355 0.5758355
```
