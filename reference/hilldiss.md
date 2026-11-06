# Hill numbers-based dissimilarity

Compute overall (multi-sample) dissimilarity metrics from the
Hill-number beta diversity following Chiu et al. (2014). These are the
complements of the similarities returned by
[`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md).
With `hierarchy`, each transition's beta is converted using its
data-dependent effective ceiling rather than the number of child units,
keeping the metrics in `[0, 1]`.

## Usage

``` r
hilldiss(
  data,
  q = c(0, 1, 2),
  metric = c("S", "C", "U", "V"),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  type = c("auto", "neutral", "phylogenetic", "functional"),
  out = c("tibble", "matrix"),
  hierarchy = NULL,
  metadata = NULL,
  weights = "equal",
  depth = NULL
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
  Defaults to all four. `"S"` (Jaccard-type) and `"V"` (Sorensen-type)
  are turnover measures, relative to the pooled assemblage and to a
  single assemblage respectively; `"U"` (Jaccard-type) and `"C"`
  (Sorensen-type; one minus the Morisita-Horn similarity at `q = 2`) are
  the corresponding overlap complements (proportions of non-shared
  units). At `q = 0`, `"U"` equals `"S"` and `"C"` equals `"V"`; at
  `q = 1`, `"U"` equals `"C"`.

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

  Output shape: `"tibble"` (default) returns a long-format `data.frame`
  with columns `q`, `metric`, `value`; `"matrix"` returns the legacy
  matrix (orders in rows, metrics in columns, dropped to a vector for a
  single metric). With `hierarchy`, the table adds a `scale` column for
  each transition above the sample scale. A `tau` grid adds a `tau`
  column. Matrix output has one row per `q` and transition, with a named
  list of matrices for a `tau` grid.

- hierarchy:

  Optional one-sided nesting formula, coarsest to finest, e.g.
  `~ region / site`. Every transition above `sample`, including `total`,
  gets its own dissimilarity metrics. Each run uses one shared tree
  `depth` or functional `tau`; a functional `tau` grid runs the full
  hierarchy once per threshold.

- metadata:

  Optional per-sample `data.frame` supplying the variables named in
  `hierarchy`; rows are matched to the count-table columns by name when
  possible, otherwise by position.

- weights:

  Sample weights for hierarchical partitioning: `"equal"` (default),
  `"abundance"` (proportional to each sample's total count), or a
  positive numeric vector with one value per sample. A named vector is
  matched to count-table columns; an unnamed vector uses column order.
  Weights are normalised to sum to one. Requires `hierarchy`.

- depth:

  Optional positive phylogenetic reference depth for hierarchical
  partitioning. On an ultrametric tree, use only branch segments within
  the last `depth` units before the tips. Must not exceed the tree
  height. `NULL` (default) uses the current shared mean-depth
  normaliser. Requires `hierarchy` and a phylogenetic `tree`.

## Value

A long-format `data.frame` of class `hill_dissimilarity` (default, with
a [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method), or a
matrix/vector of dissimilarities when `out = "matrix"`.

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
#> <hilldiv3 result: neutral>
#> 12 rows x 3 cols
#> 
#>    q metric     value
#> 1  0      S 0.3333333
#> 2  1      S 0.6081390
#> 3  2      S 0.7308320
#> 4  0      C 0.2000000
#> 5  1      C 0.5229848
#> 6  2      C 0.7308320
#> 7  0      U 0.3333333
#> 8  1      U 0.5229848
#> 9  2      U 0.5758355
#> 10 0      V 0.2000000
#> 11 1      V 0.4369251
#> 12 2      V 0.5758355
plot(hilldiss(counts))
#> dissimilarity from neutral Hill numbers of "q0", "q1", and "q2".
```
