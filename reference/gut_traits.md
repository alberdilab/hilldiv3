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
#> 1  0 ctrl01 1.903890
#> 2  1 ctrl01 1.864119
#> 3  0 ctrl02 2.011087
#> 4  1 ctrl02 1.997614
#> 5  0 ctrl03 1.962929
#> 6  1 ctrl03 1.928385
#> 7  0 ctrl04 1.947593
#> 8  1 ctrl04 1.913489
#> 9  0 ctrl05 1.864104
#> 10 1 ctrl05 1.816147
#> 11 0 ctrl06 1.914662
#> 12 1 ctrl06 1.874010
#> 13 0  trt01 1.974545
#> 14 1  trt01 1.954243
#> 15 0  trt02 1.841115
#> 16 1  trt02 1.796213
#> 17 0  trt03 1.954853
#> 18 1  trt03 1.941623
#> 19 0  trt04 1.874825
#> 20 1  trt04 1.854832
#> 21 0  trt05 1.763973
#> 22 1  trt05 1.722022
#> 23 0  trt06 1.831259
#> 24 1  trt06 1.802398
```
