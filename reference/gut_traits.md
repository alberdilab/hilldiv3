# Functional traits for the simulated gut MAGs

A trait table for the 24 MAGs in
[gut_counts](https://alberdilab.github.io/hilldiv3/reference/gut_counts.md),
mixing continuous, categorical and binary traits. Convert it to a
functional distance with
[`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md)
for the functional-diversity paths.

## Usage

``` r
gut_traits
```

## Format

A data frame with 24 rows (MAGs) and 4 columns:

- genome_size:

  Approximate genome size in Mbp (numeric).

- gc_content:

  GC content as a proportion (numeric).

- oxygen:

  Oxygen tolerance: aerobe, anaerobe or facultative (factor).

- motility:

  Motility indicator, 0/1 (integer).

## See also

[gut_counts](https://alberdilab.github.io/hilldiv3/reference/gut_counts.md),
[gut_tree](https://alberdilab.github.io/hilldiv3/reference/gut_tree.md),
[`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md)

## Examples

``` r
d <- traits2dist(gut_traits)
hilldiv(gut_counts, q = c(0, 1), dist = d)
#> Computing functional Hill numbers of "q0" and "q1".
#> <hilldiv3 result: functional>
#> 24 rows x 3 cols
#> 
#>    q sample    value
#> 1  0 ctrl01 1.768437
#> 2  1 ctrl01 1.725206
#> 3  0 ctrl02 1.794542
#> 4  1 ctrl02 1.751175
#> 5  0 ctrl03 1.785493
#> 6  1 ctrl03 1.742742
#> 7  0 ctrl04 1.798031
#> 8  1 ctrl04 1.754829
#> 9  0 ctrl05 1.798443
#> 10 1 ctrl05 1.754831
#> 11 0 ctrl06 1.780635
#> 12 1 ctrl06 1.736590
#> 13 0  trt01 1.642686
#> 14 1  trt01 1.577995
#> 15 0  trt02 1.670348
#> 16 1  trt02 1.606515
#> 17 0  trt03 1.655314
#> 18 1  trt03 1.591533
#> 19 0  trt04 1.665013
#> 20 1  trt04 1.600605
#> 21 0  trt05 1.661506
#> 22 1  trt05 1.595776
#> 23 0  trt06 1.668690
#> 24 1  trt06 1.604980
```
