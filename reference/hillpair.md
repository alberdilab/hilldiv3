# Pairwise Hill numbers-based dissimilarity

Compute dissimilarity metrics for every pair of samples, returning
distance objects suitable for ordination (e.g. NMDS, PCoA).

## Usage

``` r
hillpair(
  data,
  q = c(0, 1, 2),
  metric = c("S", "C", "U", "V"),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  type = c("auto", "neutral", "phylogenetic", "functional"),
  out = c("dist", "tibble"),
  parallel = FALSE
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

- out:

  Output type: `"dist"` (default) returns a `dist` object per requested
  metric/order combination; `"tibble"` returns a long-format table.

- parallel:

  Logical; if `TRUE` and `furrr` is installed, compute pairs in
  parallel.

## Value

For `out = "dist"`, a named list of `dist` objects (one per
order/metric, named e.g. `"q0S"`), collapsed to a single `dist` when
only one combination is requested. For `out = "tibble"`, a long-format
`data.frame` with columns `first`, `second`, `q`, `metric`, `value`.

## Details

The type-specific structure (per-sample normalisation, the tree
traversal or the functional similarity product) is computed **once**
over all samples via the partitioning engine; each pair then only
combines its two precomputed columns into beta, which is turned into the
requested overlap metrics. The maths are therefore identical to
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md)
on two samples, without re-running the full engine per pair. When
`parallel = TRUE` and the `furrr` package is installed, pairs are
computed in parallel via the active `future` plan. A `progressr`
progress bar is reported when that package is installed and a handler is
active.

## See also

[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md),
[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1, 3, 4, 0, 6, 2, 7), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"),
                                 c("s1", "s2", "s3", "s4")))
hillpair(counts, q = 1, metric = "C")
#> Computing neutral pairwise dissimilarity for 6 sample pairs.
#>            s1         s2         s3
#> s2 0.52298484                      
#> s3 0.47119926 0.08924381           
#> s4 0.09902103 0.29281463 0.33948330
```
