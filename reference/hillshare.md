# Allocate Hill power-sum weight to taxa or sets

Calculate the share of a community's neutral, phylogenetic, or
functional Hill power sum assigned to each taxon or a specified set of
taxa. Shares range from zero to one within a nonempty sample. They
describe an additive allocation of the power sum, not a percentage of
the effective Hill number.

## Usage

``` r
hillshare(
  data,
  q = c(0, 1, 2),
  tree = NULL,
  dist = NULL,
  tau = NULL,
  type = c("auto", "neutral", "phylogenetic", "functional"),
  sets = NULL,
  allocation = "qpower"
)
```

## Arguments

- data:

  Counts: a numeric vector or a taxa-by-samples matrix/data frame. The
  input adapters accepted by
  [`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)
  are also supported.

- q:

  Nonnegative diversity orders. Defaults to c(0, 1, 2).

- tree:

  Optional rooted phylogenetic tree of class phylo.

- dist:

  Optional functional distance matrix over the taxa.

- tau:

  Functional distance threshold. Defaults to max(dist).

- type:

  Diversity type(s). "auto" (default) computes neutral shares and also
  phylogenetic and functional shares when their inputs are supplied.
  Explicit types restrict the output as in
  [`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md).

- sets:

  A named list whose elements are character vectors of taxon names. Each
  element can contain one taxon or several. By default, every taxon is
  returned as a singleton set. Sets may overlap, and each set's share is
  the sum of its members' allocations. An unnamed character vector is
  accepted as a request for singleton taxa.

- allocation:

  Allocation rule. The first iteration supports "qpower".

## Value

A data frame with columns q, sample, type, set, share, and allocation.
Shares for distinct singleton taxa sum to one in each nonempty sample,
diversity type, and order.

## Details

Let p_i be the relative abundance of taxon i. Neutral weight is p_i^q,
normalised by its sum over present taxa.

For a phylogenetic branch b with length L_b and descendant abundance
a_b, branch weight is L_b a_b^q. That weight is divided among its
present descendants in proportion to p_i^q. At q = 0, present
descendants share the branch equally. Branch weights and taxon
allocations are normalised by the sum of all branch weights.

For functional diversity, similarity is s_ij = 1 - min(d_ij, tau) / tau
and a_i = sum_j s_ij p_j. The existing functional Hill power sum in
hilldiv3 is sum_i p_i a_i^(q-1), so taxon i receives p_i a_i^(q-1) /
sum_j p_j a_j^(q-1). With identity similarity this reduces to the
neutral rule. At q = 1, functional shares equal p_i, although the
effective functional Hill number still depends on similarities.

An empty sample has no power sum to allocate, so its shares are NA.

## References

Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
based on Hill numbers. Phil. Trans. R. Soc. B, 365, 3599-3609.

Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
of species similarity. Ecology, 93, 477-489.

## See also

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md),
[`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md)

## Examples

``` r
counts <- c(A = 394, B = 202, C = 354, D = 50)
tree <- ape::read.tree(
  text = "(((A:0.300,B:0.300):0.185,C:0.485):0.009,D:0.494);"
)
hillshare(counts, q = 3, tree = tree,
          sets = list(C = "C", ABC = c("A", "B", "C")))
#>   q  sample         type set     share allocation
#> 1 3 sample1      neutral   C 0.3895073     qpower
#> 2 3 sample1      neutral ABC 0.9989025     qpower
#> 3 3 sample1 phylogenetic   C 0.2746866     qpower
#> 4 3 sample1 phylogenetic ABC 0.9993084     qpower
d <- as.matrix(stats::dist(c(A = 0, B = 0.5, C = 1, D = 2)))
hillshare(counts, q = c(0, 1, 2), dist = d, type = "functional",
          sets = list(C = "C", ABC = c("A", "B", "C")))
#>   q  sample       type set     share allocation
#> 1 0 sample1 functional   C 0.3304382     qpower
#> 2 0 sample1 functional ABC 0.8776435     qpower
#> 3 1 sample1 functional   C 0.3540000     qpower
#> 4 1 sample1 functional ABC 0.9500000     qpower
#> 5 2 sample1 functional   C 0.3613380     qpower
#> 6 2 sample1 functional ABC 0.9805325     qpower
```
