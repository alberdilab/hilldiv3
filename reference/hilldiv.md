# Hill numbers computation

Compute neutral, phylogenetic or functional Hill numbers (alpha
diversity) from a single sample or a count table. The diversity type is
inferred from the inputs: counts only -\> neutral; counts + `tree` -\>
phylogenetic; counts + `dist` -\> functional.

## Usage

``` r
hilldiv(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL)
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

## Value

A matrix of Hill numbers with diversity orders in rows (`q0`, `q1`, ...)
and samples in columns.

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
#>          s1       s2
#> q0 2.000000 3.000000
#> q1 1.889882 2.137309
#> q2 1.800000 1.753623
hilldiv(counts, q = c(0, 1, 2))
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#>          s1       s2
#> q0 2.000000 3.000000
#> q1 1.889882 2.137309
#> q2 1.800000 1.753623
```
