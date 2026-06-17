# Hill numbers redundancy

Estimate phylogenetic or functional redundancy by fitting the saturating
relationship between neutral diversity and phylogenetic/functional
diversity across samples: `y = -a * 2^(-x / b) + c`. Redundancy is
summarised as `1 - b / max(x)`.

## Usage

``` r
hillred(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL)
```

## Arguments

- data:

  A count table (taxa x samples); requires either `tree` or `dist`.

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

A matrix with columns `redundancy`, `a`, `b`, `c`, one row per `q`.

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)
