# Convert a trait table into a distance matrix

Build a pairwise functional distance matrix from a table of taxon
traits, suitable as the `dist` argument of
[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)
and friends.

## Usage

``` r
traits2dist(traits, method = c("gower", "euclidean", "manhattan"))
```

## Arguments

- traits:

  A table with taxa (OTUs/ASVs/MAGs) in rows and traits in columns.
  Traits may be continuous, binary or proportional.

- method:

  Distance metric passed to
  [`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html):
  `"gower"` (default), `"euclidean"` or `"manhattan"`.

## Value

A numeric distance matrix.

## Examples

``` r
traits <- data.frame(body = c(1, 0.2, 0.9), diet = c(0L, 1L, 1L),
                     row.names = c("t1", "t2", "t3"))
traits2dist(traits)
#>        t1     t2     t3
#> t1 0.0000 1.0000 0.5625
#> t2 1.0000 0.0000 0.4375
#> t3 0.5625 0.4375 0.0000
```
