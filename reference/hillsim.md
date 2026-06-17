# Hill numbers-based similarity

Compute overall similarity metrics from the Hill-number beta diversity
(Chiu et al. 2014). These are `1 -` the dissimilarities from
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md).

## Usage

``` r
hillsim(
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

A matrix of similarities (diversity orders in rows, metrics in columns),
or a vector if a single metric is requested.

## See also

[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
hillsim(counts)
#> similarity from neutral Hill numbers of "q0", "q1", and "q2".
#>            S         C         U         V
#> q0 0.6666667 0.8000000 0.6666667 0.8000000
#> q1 0.3918610 0.4770152 0.4770152 0.5630749
#> q2 0.2691680 0.2691680 0.4241645 0.4241645
```
