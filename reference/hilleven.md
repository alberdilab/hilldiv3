# Hill-number evenness

Evenness expressed through Hill numbers as the ratio of diversity of
order `q` to richness (`qD / 0D`), which ranges from 0 to 1.

## Usage

``` r
hilleven(data, q = c(1, 2), tree = NULL, dist = NULL, tau = NULL)
```

## Arguments

- data:

  Counts: a numeric vector (one sample), a matrix/data.frame (taxa x
  samples), a `phyloseq` object or a `TreeSummarizedExperiment`.

- q:

  Numeric vector of diversity orders (\> 0 are meaningful for evenness).
  Defaults to `c(1, 2)`.

- tree:

  A phylogenetic tree of class `phylo` whose tip labels match the taxa
  in `data`.

- dist:

  A functional distance matrix (or `dist`) over the taxa.

- tau:

  Optional functional distance threshold. Defaults to `max(dist)`.

## Value

A matrix of evenness values (orders in rows, samples in columns).

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
hilleven(counts)
#> Computing neutral evenness of "q1" and "q2".
#>           s1        s2
#> q1 0.9449408 0.7124362
#> q2 0.9000000 0.5845411
```
