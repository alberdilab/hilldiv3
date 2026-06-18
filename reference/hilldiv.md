# Hill numbers computation

Compute neutral, phylogenetic or functional Hill numbers (alpha
diversity) from a single sample or a count table. The diversity type is
inferred from the inputs: counts only -\> neutral; counts + `tree` -\>
phylogenetic; counts + `dist` -\> functional.

## Usage

``` r
hilldiv(
  data,
  q = c(0, 1, 2),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  type = c("auto", "neutral", "phylogenetic", "functional"),
  reference = c("pool", "sample"),
  out = c("tibble", "matrix")
)
```

## Arguments

- data:

  Counts: a numeric vector (one sample), a matrix/data.frame (taxa x
  samples), a `phyloseq` object or a `TreeSummarizedExperiment`.

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
  (matching hilldiv2's `multi` behaviour); `"sample"` reads each sample
  at its own depth `T_j` (effective lineages at that sample's depth).
  The two coincide on ultrametric trees. This reference depth is
  intentionally *not* offered by
  [`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md):
  in a partition `T` is fixed at the mean per-sample depth of Chiu et
  al. (2014), the unique value for which `gamma / alpha` is a valid
  decomposition with `beta` in `[1, N]`.

- out:

  Output shape: `"tibble"` (default) returns a long-format `data.frame`
  with columns `q`, `sample`, `value` and
  [`print()`](https://rdrr.io/r/base/print.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods; `"matrix"` returns a matrix with samples in rows and
  diversity orders (`q0`, `q1`, ...) in columns.

## Value

A long-format `data.frame` of class `hill_diversity` (default), or a
matrix of Hill numbers with samples in rows and diversity orders (`q0`,
`q1`, ...) in columns when `out = "matrix"`.

## References

Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
based on Hill numbers. Phil. Trans. R. Soc. B, 365, 3599-3609.  
  
Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
numbers to DNA-based diversity analyses. Mol. Ecol. Resour., 19,
804-817.

## See also

[`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md),
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md),
[`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md)

## Examples

``` r
counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
hilldiv(counts)
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 2 samples.
#> <hilldiv3 result: neutral>
#> 6 rows x 3 cols
#> 
#>   q sample    value
#> 1 0     s1 2.000000
#> 2 1     s1 1.889882
#> 3 2     s1 1.800000
#> 4 0     s2 3.000000
#> 5 1     s2 2.137309
#> 6 2     s2 1.753623
hilldiv(counts, q = c(0, 1, 2))
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 2 samples.
#> <hilldiv3 result: neutral>
#> 6 rows x 3 cols
#> 
#>   q sample    value
#> 1 0     s1 2.000000
#> 2 1     s1 1.889882
#> 3 2     s1 1.800000
#> 4 0     s2 3.000000
#> 5 1     s2 2.137309
#> 6 2     s2 1.753623
plot(hilldiv(counts, q = c(0, 1, 2)))
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 2 samples.
```
