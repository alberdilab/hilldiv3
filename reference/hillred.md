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
  reference = c("pool", "sample"),
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

- reference:

  Reference tree depth for *phylogenetic* Hill numbers (ignored for
  neutral and functional types). `"pool"` (default) reads every sample
  at one common depth `T = mean(T_j)`, so values share a comparable axis
  across samples; `"sample"` reads each sample at its own depth `T_j`
  (effective lineages at that sample's depth). The two coincide on
  ultrametric trees. This reference depth is intentionally *not* offered
  by
  [`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md):
  in a partition `T` is fixed at the mean per-sample depth of Chiu et
  al. (2014), the unique value for which `gamma / alpha` is a valid
  decomposition with `beta` in `[1, N]`.

- out:

  Output shape: `"tibble"` (default) returns a `data.frame` with one row
  per `q` and columns `q`, `redundancy`, `a`, `b`, `c`; `"matrix"`
  returns the legacy matrix (orders in rows).

## Value

A `data.frame` of class `hill_redundancy` (default) with a
[plot()](https://alberdilab.github.io/hilldiv3/reference/plot.hill_redundancy.md)
method, or a matrix with columns `redundancy`, `a`, `b`, `c` (one row
per `q`) when `out = "matrix"`. The tibble carries the per-sample
neutral and phylogenetic/functional diversity used for the fit as a
`"hill_fit"` attribute, which the plot method draws.

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md),
[`plot.hill_redundancy()`](https://alberdilab.github.io/hilldiv3/reference/plot.hill_redundancy.md)

## Examples

``` r
d <- traits2dist(gut_traits)
red <- hillred(gut_counts, dist = d)
#> Warning: Redundancy for "q0" could not be estimated: singular gradient matrix at initial
#> parameter estimates
red
#> <hilldiv3 result: functional>
#> 3 rows x 5 cols
#> 
#>   q redundancy         a        b        c
#> 1 0         NA        NA       NA       NA
#> 2 1  0.6957595 1.1191121 4.770952 2.083805
#> 3 2  0.6024721 0.5929264 5.143371 2.057040
plot(red)
```
