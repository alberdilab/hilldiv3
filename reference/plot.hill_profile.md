# Plot a diversity profile

Base-graphics plot of a
[`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md)
result: one line per sample showing the Hill number against the
diversity order `q`.

## Usage

``` r
# S3 method for class 'hill_profile'
plot(x, ...)
```

## Arguments

- x:

  A `hill_profile` object from
  [`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md).

- ...:

  Further arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The `hill_profile` object, invisibly.
