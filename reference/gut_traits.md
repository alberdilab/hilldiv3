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
#> Computing "neutral" and "functional" Hill numbers of "q0" and "q1".
#> ℹ 24 taxa across 12 samples.
#> <hilldiv3 result: neutral, functional>
#> 48 rows x 4 cols
#> 
#>    q sample       type     value
#> 1  0 ctrl01    neutral 24.000000
#> 2  1 ctrl01    neutral 12.329592
#> 3  0 ctrl02    neutral 24.000000
#> 4  1 ctrl02    neutral 14.784342
#> 5  0 ctrl03    neutral 24.000000
#> 6  1 ctrl03    neutral 11.832847
#> 7  0 ctrl04    neutral 23.000000
#> 8  1 ctrl04    neutral 11.266414
#> 9  0 ctrl05    neutral 24.000000
#> 10 1 ctrl05    neutral 10.035886
#> 11 0 ctrl06    neutral 24.000000
#> 12 1 ctrl06    neutral 13.277945
#> 13 0  trt01    neutral 24.000000
#> 14 1  trt01    neutral 11.888942
#> 15 0  trt02    neutral 24.000000
#> 16 1  trt02    neutral  8.350790
#> 17 0  trt03    neutral 24.000000
#> 18 1  trt03    neutral 15.681517
#> 19 0  trt04    neutral 23.000000
#> 20 1  trt04    neutral 11.671685
#> 21 0  trt05    neutral 23.000000
#> 22 1  trt05    neutral  9.043005
#> 23 0  trt06    neutral 23.000000
#> 24 1  trt06    neutral 10.901319
#> 25 0 ctrl01 functional  1.903890
#> 26 1 ctrl01 functional  1.864119
#> 27 0 ctrl02 functional  2.011087
#> 28 1 ctrl02 functional  1.997614
#> 29 0 ctrl03 functional  1.962929
#> 30 1 ctrl03 functional  1.928385
#> 31 0 ctrl04 functional  1.947593
#> 32 1 ctrl04 functional  1.913489
#> 33 0 ctrl05 functional  1.864104
#> 34 1 ctrl05 functional  1.816147
#> 35 0 ctrl06 functional  1.914662
#> 36 1 ctrl06 functional  1.874010
#> 37 0  trt01 functional  1.974545
#> 38 1  trt01 functional  1.954243
#> 39 0  trt02 functional  1.841115
#> 40 1  trt02 functional  1.796213
#> 41 0  trt03 functional  1.954853
#> 42 1  trt03 functional  1.941623
#> 43 0  trt04 functional  1.874825
#> 44 1  trt04 functional  1.854832
#> 45 0  trt05 functional  1.763973
#> 46 1  trt05 functional  1.722022
#> 47 0  trt06 functional  1.831259
#> 48 1  trt06 functional  1.802398
```
