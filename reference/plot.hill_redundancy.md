# Plot a redundancy fit

Base-graphics plot of a
[`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md)
result. For each diversity order `q` it shows the per-sample neutral
diversity (x) against phylogenetic/functional diversity (y), overlaid
with the fitted saturating curve `y = -a * 2^(-x / b) + c`. A curve that
bends sharply and plateaus well below the points' spread indicates high
redundancy; a near-linear fit indicates low redundancy. This mirrors the
profile plot of
[`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md).

## Usage

``` r
# S3 method for class 'hill_redundancy'
plot(x, ...)
```

## Arguments

- x:

  A `hill_redundancy` object from
  [`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md).

- ...:

  Further arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The `hill_redundancy` object, invisibly.

## See also

[`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md),
[`plot.hill_profile()`](https://alberdilab.github.io/hilldiv3/reference/plot.hill_profile.md)
