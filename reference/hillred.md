# Hill numbers redundancy

Estimate phylogenetic or functional redundancy by fitting the saturating
relationship between neutral diversity and phylogenetic/functional
diversity across samples: `y = -a * 2^(-x / b) + c`. Redundancy is
summarised as `1 - b / max(x)`.

## Usage

``` r
hillred(
  data,
  q = c(0, 1, 2),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  type = c("auto", "phylogenetic", "functional"),
  out = c("tibble", "matrix")
)
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

- type:

  Diversity type: `"auto"` (default) infers it from the inputs (counts
  only -\> neutral, `+tree` -\> phylogenetic, `+dist` -\> functional);
  an explicit `"neutral"`, `"phylogenetic"` or `"functional"` asserts
  the type and is validated against the inputs (e.g. `"phylogenetic"`
  requires a `tree`; `"neutral"` ignores any tree/dist carried by the
  object).

- out:

  Output shape: `"tibble"` (default) returns a `data.frame` with one row
  per `q` and columns `q`, `redundancy`, `a`, `b`, `c`; `"matrix"`
  returns the legacy matrix (orders in rows).

## Value

A `data.frame` of class `hill_redundancy` (default), or a matrix with
columns `redundancy`, `a`, `b`, `c` (one row per `q`) when
`out = "matrix"`.

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)
