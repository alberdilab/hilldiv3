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
#>    q sample    value
#> 1  0 ctrl01 4.206600
#> 2  1 ctrl01 2.361459
#> 3  0 ctrl02 4.206600
#> 4  1 ctrl02 2.357467
#> 5  0 ctrl03 4.206600
#> 6  1 ctrl03 2.376611
#> 7  0 ctrl04 4.143380
#> 8  1 ctrl04 2.381757
#> 9  0 ctrl05 4.143380
#> 10 1 ctrl05 2.392420
#> 11 0 ctrl06 4.143380
#> 12 1 ctrl06 2.325091
#> 13 0  trt01 4.157958
#> 14 1  trt01 2.049575
#> 15 0  trt02 4.206600
#> 16 1  trt02 2.085843
#> 17 0  trt03 4.206600
#> 18 1  trt03 2.073179
#> 19 0  trt04 4.157958
#> 20 1  trt04 2.079583
#> 21 0  trt05 4.206600
#> 22 1  trt05 2.080046
#> 23 0  trt06 4.157958
#> 24 1  trt06 2.085294
```
