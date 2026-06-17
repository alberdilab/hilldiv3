# Hill numbers diversity partitioning

Partition neutral, phylogenetic or functional Hill-number diversity into
alpha, gamma and beta components across a set of samples.

## Usage

``` r
hillpart(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL)
```

## Arguments

- data:

  A count table (taxa x samples) or a supported object; a single sample
  is not meaningful for partitioning.

- q:

  Numeric vector of diversity orders (\>= 0). Defaults to `c(0, 1, 2)`
  (richness, Shannon, Simpson).

- tree:

  A phylogenetic tree of class `phylo` whose tip labels match the taxa
  in `data`.

- dist:

  A functional distance matrix (or `dist`) over the taxa.

- tau:

  Optional functional distance threshold. Defaults to `max(dist)`.

## Value

A matrix with columns `alpha`, `gamma`, `beta` and diversity orders in
rows.

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md),
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md),
[`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
hillpart(counts)
#> Partitioning neutral Hill numbers of "q0", "q1", and "q2".
#>       alpha    gamma     beta
#> q0 2.500000 3.000000 1.200000
#> q1 2.009791 2.887919 1.436925
#> q2 1.776509 2.799486 1.575835
```
