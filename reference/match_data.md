# Match and align a count table to a tree or distance matrix

Subsets and reorders a count table so that its taxa match those of a
phylogenetic tree or a functional distance matrix, dropping taxa absent
from the reference. This realises the `match_data()` helper that
hilldiv2's documentation referred to but never provided.

## Usage

``` r
match_data(data, tree = NULL, dist = NULL)
```

## Arguments

- data:

  A count matrix/data.frame (taxa x samples) with row names.

- tree:

  A `phylo` tree (optional).

- dist:

  A distance matrix (optional).

## Value

The count matrix restricted to and ordered by the shared taxa.

## Examples

``` r
counts <- matrix(1:6, nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
tree <- ape::read.tree(text = "((t1:1,t2:1):1,t4:2);")
match_data(counts, tree = tree)
#> Dropped 1 taxon from `data` not in the tree tips.
#> 1 taxon in the tree tips has no counts; prune it before downstream analysis.
#>    s1 s2
#> t1  1  4
#> t2  2  5
```
