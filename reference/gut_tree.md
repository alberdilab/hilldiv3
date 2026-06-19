# Phylogeny for the simulated gut MAGs

An ultrametric coalescent tree over the 24 MAGs in
[gut_counts](https://alberdilab.github.io/hilldiv3/reference/gut_counts.md),
scaled to unit depth. Use it for the phylogenetic-diversity paths.

## Usage

``` r
gut_tree
```

## Format

A `phylo` object (see the ape package) with 24 tips whose labels match
the rows of
[gut_counts](https://alberdilab.github.io/hilldiv3/reference/gut_counts.md).

## See also

[gut_counts](https://alberdilab.github.io/hilldiv3/reference/gut_counts.md),
[gut_traits](https://alberdilab.github.io/hilldiv3/reference/gut_traits.md)

## Examples

``` r
hilldiv(gut_counts, q = c(0, 1), tree = gut_tree)
#> Computing "neutral" and "phylogenetic" Hill numbers of "q0" and "q1".
#> ℹ 24 taxa across 12 samples.
#> <hilldiv3 result: neutral, phylogenetic>
#> 48 rows x 4 cols
#> 
#>    q sample         type     value
#> 1  0 ctrl01      neutral 24.000000
#> 2  1 ctrl01      neutral 12.329592
#> 3  0 ctrl02      neutral 24.000000
#> 4  1 ctrl02      neutral 14.784342
#> 5  0 ctrl03      neutral 24.000000
#> 6  1 ctrl03      neutral 11.832847
#> 7  0 ctrl04      neutral 23.000000
#> 8  1 ctrl04      neutral 11.266414
#> 9  0 ctrl05      neutral 24.000000
#> 10 1 ctrl05      neutral 10.035886
#> 11 0 ctrl06      neutral 24.000000
#> 12 1 ctrl06      neutral 13.277945
#> 13 0  trt01      neutral 24.000000
#> 14 1  trt01      neutral 11.888942
#> 15 0  trt02      neutral 24.000000
#> 16 1  trt02      neutral  8.350790
#> 17 0  trt03      neutral 24.000000
#> 18 1  trt03      neutral 15.681517
#> 19 0  trt04      neutral 23.000000
#> 20 1  trt04      neutral 11.671685
#> 21 0  trt05      neutral 23.000000
#> 22 1  trt05      neutral  9.043005
#> 23 0  trt06      neutral 23.000000
#> 24 1  trt06      neutral 10.901319
#> 25 0 ctrl01 phylogenetic  3.757442
#> 26 1 ctrl01 phylogenetic  2.156361
#> 27 0 ctrl02 phylogenetic  3.757442
#> 28 1 ctrl02 phylogenetic  2.105705
#> 29 0 ctrl03 phylogenetic  3.757442
#> 30 1 ctrl03 phylogenetic  2.035022
#> 31 0 ctrl04 phylogenetic  3.750466
#> 32 1 ctrl04 phylogenetic  2.263603
#> 33 0 ctrl05 phylogenetic  3.757442
#> 34 1 ctrl05 phylogenetic  1.881365
#> 35 0 ctrl06 phylogenetic  3.757442
#> 36 1 ctrl06 phylogenetic  2.159882
#> 37 0  trt01 phylogenetic  3.757442
#> 38 1  trt01 phylogenetic  2.135261
#> 39 0  trt02 phylogenetic  3.757442
#> 40 1  trt02 phylogenetic  1.649469
#> 41 0  trt03 phylogenetic  3.757442
#> 42 1  trt03 phylogenetic  2.212881
#> 43 0  trt04 phylogenetic  3.742870
#> 44 1  trt04 phylogenetic  1.853561
#> 45 0  trt05 phylogenetic  3.750466
#> 46 1  trt05 phylogenetic  1.906460
#> 47 0  trt06 phylogenetic  3.754523
#> 48 1  trt06 phylogenetic  1.996142
```
