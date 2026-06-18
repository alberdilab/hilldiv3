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
#> Computing phylogenetic Hill numbers of "q0" and "q1".
#> ℹ 24 taxa across 12 samples.
#> <hilldiv3 result: phylogenetic>
#> 24 rows x 3 cols
#> 
#>    q sample    value
#> 1  0 ctrl01 3.757442
#> 2  1 ctrl01 2.156361
#> 3  0 ctrl02 3.757442
#> 4  1 ctrl02 2.105705
#> 5  0 ctrl03 3.757442
#> 6  1 ctrl03 2.035022
#> 7  0 ctrl04 3.750466
#> 8  1 ctrl04 2.263603
#> 9  0 ctrl05 3.757442
#> 10 1 ctrl05 1.881365
#> 11 0 ctrl06 3.757442
#> 12 1 ctrl06 2.159882
#> 13 0  trt01 3.757442
#> 14 1  trt01 2.135261
#> 15 0  trt02 3.757442
#> 16 1  trt02 1.649469
#> 17 0  trt03 3.757442
#> 18 1  trt03 2.212881
#> 19 0  trt04 3.742870
#> 20 1  trt04 1.853561
#> 21 0  trt05 3.750466
#> 22 1  trt05 1.906460
#> 23 0  trt06 3.754523
#> 24 1  trt06 1.996142
```
