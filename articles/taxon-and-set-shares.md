# Taxon and set shares of Hill power sums

``` r

library(hilldiv3)
```

[`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md)
assigns a fraction of a community’s Hill **power sum** to each taxon, or
to a user-defined set of taxa. The fraction ranges from zero to one
within a sample, diversity type and order `q`. It is a descriptive
allocation: it does not provide a test of enrichment between groups. It
is also not a fraction of the effective Hill number returned by
[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md).

The first allocation method is `"qpower"` (the default). You can use it
with neutral, phylogenetic and functional Hill numbers. With
`type = "auto"`, adding a tree or functional distance matrix returns the
corresponding shares alongside neutral shares.

## A small community

Consider three taxa with counts 5, 3 and 2. We will report the share
assigned to taxon C and to the set A+B+C. Define sets as a named list;
the name becomes the `set` label in the output.

``` r

counts <- c(A = 5, B = 3, C = 2)
sets <- list(C = "C", ABC = c("A", "B", "C"))
hillshare(counts, q = c(0, 1, 2), sets = sets)
#>   q  sample    type set     share allocation
#> 1 0 sample1 neutral   C 0.3333333     qpower
#> 2 0 sample1 neutral ABC 1.0000000     qpower
#> 3 1 sample1 neutral   C 0.2000000     qpower
#> 4 1 sample1 neutral ABC 1.0000000     qpower
#> 5 2 sample1 neutral   C 0.1052632     qpower
#> 6 2 sample1 neutral ABC 1.0000000     qpower
```

For neutral diversity, the weight of taxon `i` is its relative abundance
raised to `q`, divided by the sum of those weights. At `q = 0`, each
present taxon receives one third. At `q = 1`, C receives its abundance
share, 0.2. At `q = 2`, its share is `0.2^2 / (0.5^2 + 0.3^2 + 0.2^2)`,
about 0.105. ABC contains every taxon, so its share is always one.

## Shared phylogenetic branches

Phylogenetic shares also account for the branch length connecting each
taxon to the tree. A branch of length $`L_b`$ with descendant abundance
$`a_b`$ has weight $`L_b a_b^q`$.
[`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md)
divides this weight among the present descendants in proportion to their
individual abundance raised to `q`, then normalises across all branches.

``` r

tree <- ape::read.tree(text = "((A:1,B:1):1,C:2);")
hillshare(counts, q = c(0, 1, 2), tree = tree,
          type = "phylogenetic", sets = sets)
#>   q  sample         type set     share allocation
#> 1 0 sample1 phylogenetic   C 0.4000000     qpower
#> 2 0 sample1 phylogenetic ABC 1.0000000     qpower
#> 3 1 sample1 phylogenetic   C 0.2000000     qpower
#> 4 1 sample1 phylogenetic ABC 1.0000000     qpower
#> 5 2 sample1 phylogenetic   C 0.0754717     qpower
#> 6 2 sample1 phylogenetic ABC 1.0000000     qpower
```

At `q = 0`, C has its own branch of length 2 out of a total branch
length of 5, so its share is 0.4. The branch shared by A and B is
divided equally between them. At `q = 1` on this equal-depth tree, the
shares equal the taxa’s relative abundances.

## Functional similarity

Functional shares use a distance matrix over the taxa and the same
similarity calculation as
[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md).
With threshold `tau` (default: the largest distance), similarity is
$`s_{ij}=1-\min(d_{ij},\tau)/\tau`$. Taxon `i` has ordinariness
$`a_i=\sum_j s_{ij}p_j`$, where $`p_j`$ is relative abundance. Its
power-sum weight is $`p_i a_i^{q-1}`$, which
[`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md)
normalises over taxa.

``` r

fdist <- matrix(c(0, 0.5, 1,
                  0.5, 0, 0.5,
                  1, 0.5, 0), 3,
                dimnames = list(names(counts), names(counts)))
hillshare(counts, q = c(0, 1, 2), dist = fdist,
          type = "functional", sets = sets)
#>   q  sample       type set     share allocation
#> 1 0 sample1 functional   C 0.3170732     qpower
#> 2 0 sample1 functional ABC 1.0000000     qpower
#> 3 1 sample1 functional   C 0.2000000     qpower
#> 4 1 sample1 functional ABC 1.0000000     qpower
#> 5 2 sample1 functional   C 0.1186441     qpower
#> 6 2 sample1 functional ABC 1.0000000     qpower
```

Here the ordinariness values are 0.65 for A and B and 0.35 for C. At
`q = 2`, C has weight $`0.2 \times 0.35=0.07`$; the sum of all three
weights is 0.59, so C receives about 0.119. At `q = 1`, functional
shares equal relative abundance, even though the *effective functional
diversity* still depends on the distances. If all taxa are functionally
distinct, the functional rule reduces to the neutral one.

## Several samples and sets

By default, `sets = NULL` returns one row per taxon; their shares sum to
one within every nonempty sample, type and `q`. For a selected group,
pass its taxon names as one element of `sets`. Sets may overlap, so
shares of overlapping sets should not be added together. A taxon absent
from a sample has share zero; an all-zero sample has undefined (`NA`)
shares.

``` r

samples <- cbind(s1 = counts, s2 = c(A = 2, B = 3, C = 5))
hillshare(samples, q = 2, tree = tree, dist = fdist,
          sets = list(C = "C", AB = c("A", "B")))
#>    q sample         type set     share allocation
#> 1  2     s1      neutral   C 0.1052632     qpower
#> 2  2     s1      neutral  AB 0.8947368     qpower
#> 3  2     s2      neutral   C 0.6578947     qpower
#> 4  2     s2      neutral  AB 0.3421053     qpower
#> 5  2     s1 phylogenetic   C 0.0754717     qpower
#> 6  2     s1 phylogenetic  AB 0.9245283     qpower
#> 7  2     s2 phylogenetic   C 0.5681818     qpower
#> 8  2     s2 phylogenetic  AB 0.4318182     qpower
#> 9  2     s1   functional   C 0.1186441     qpower
#> 10 2     s1   functional  AB 0.8813559     qpower
#> 11 2     s2   functional   C 0.5508475     qpower
#> 12 2     s2   functional  AB 0.4491525     qpower
```

Each sample is normalised separately. If samples belong to treatment
groups, you can summarise their shares by group after the call. Such a
summary describes how the power-sum allocation differs between groups;
uncertainty and statistical significance require a separate analysis.
