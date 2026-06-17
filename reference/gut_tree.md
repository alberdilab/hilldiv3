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
#> <hilldiv3 result: phylogenetic>
#> 24 rows x 3 cols
#> 
#>    q sample     value
#> 1  0 ctrl01 46.000000
#> 2  1 ctrl01  6.267857
#> 3  0 ctrl02 46.000000
#> 4  1 ctrl02  6.240741
#> 5  0 ctrl03 46.000000
#> 6  1 ctrl03  6.295456
#> 7  0 ctrl04 45.000000
#> 8  1 ctrl04  6.295552
#> 9  0 ctrl05 45.000000
#> 10 1 ctrl05  6.318648
#> 11 0 ctrl06 45.000000
#> 12 1 ctrl06  6.188841
#> 13 0  trt01 45.000000
#> 14 1  trt01  5.600712
#> 15 0  trt02 46.000000
#> 16 1  trt02  5.675155
#> 17 0  trt03 46.000000
#> 18 1  trt03  5.645288
#> 19 0  trt04 45.000000
#> 20 1  trt04  5.667823
#> 21 0  trt05 46.000000
#> 22 1  trt05  5.659486
#> 23 0  trt06 45.000000
#> 24 1  trt06  5.680068
```
