# Preparing your data

``` r

library(hilldiv3)
```

Every `hilldiv3` function accepts the same kinds of input and shares one
validation/alignment layer, so the data-wrangling you do once carries
through the whole workflow. This article covers the accepted input
formats, how counts are matched to a tree or distance matrix, and the
helpers
[`tss()`](https://alberdilab.github.io/hilldiv3/reference/tss.md),
[`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md)
and
[`match_data()`](https://alberdilab.github.io/hilldiv3/reference/match_data.md).

## The count table

The core input is a **count table** with taxa (OTUs/ASVs/MAGs) in rows
and samples in columns. Row names are taxa, column names are samples.

``` r

counts <- matrix(
  c(10, 0, 5,
     2, 8, 1,
     3, 4, 0,
     6, 2, 7),
  nrow = 3,
  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3", "s4"))
)
counts
#>    s1 s2 s3 s4
#> t1 10  2  3  6
#> t2  0  8  4  2
#> t3  5  1  0  7
```

## Accepted input objects

You can pass any of the following — they are all coerced internally to
the same representation, so you call the functions the same way:

- a **numeric matrix** (taxa × samples);
- a **data frame** or **tibble** of counts with taxa as row names;
- a **numeric vector** for a single sample (names are the taxa);
- a [`phyloseq`](https://joey711.github.io/phyloseq/) object — its OTU
  table is read, transposed if needed, and an embedded phylogeny is
  picked up automatically;
- a
  [`TreeSummarizedExperiment`](https://bioconductor.org/packages/TreeSummarizedExperiment/)
  — its assay is used as counts and its `rowTree` as the phylogeny.

``` r

# A single sample as a named vector.
hilldiv(c(t1 = 10, t2 = 2, t3 = 3))
#> Computing neutral Hill numbers of "q0", "q1", and
#> "q2".
#> <hilldiv3 result: neutral>
#> 3 rows x 3 cols
#> 
#>   q  sample    value
#> 1 0 sample1 3.000000
#> 2 1 sample1 2.365174
#> 3 2 sample1 1.991150
```

When you pass a `phyloseq` or `TreeSummarizedExperiment` that already
carries a tree, phylogenetic diversity works without supplying `tree`
yourself:

``` r

library(phyloseq)
ps <- phyloseq(
  otu_table(counts, taxa_are_rows = TRUE),
  phy_tree(ape::read.tree(text = "((t1:1,t2:1):1,t3:2);"))
)
hilldiv(ps)        # phylogenetic, tree taken from the object
#> Computing phylogenetic Hill numbers of "q0", "q1", and
#> "q2".
#> <hilldiv3 result: phylogenetic>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 1.500000
#> 2  1     s1 1.500000
#> 3  2     s1 1.500000
#> 4  0     s2 2.000000
#> 5  1     s2 1.598520
#> 6  2     s2 1.406977
#> 7  0     s3 1.500000
#> 8  1     s3 1.406992
#> 9  2     s3 1.324324
#> 10 0     s4 2.000000
#> 11 1     s4 1.677678
#> 12 2     s4 1.500000
```

## Normalisation with `tss()`

Hill numbers are computed from relative abundances, and the `hill*`
functions normalise internally — you do **not** need to pre-normalise
for them.
[`tss()`](https://alberdilab.github.io/hilldiv3/reference/tss.md)
(total-sum scaling) is exposed for when you want normalised values
yourself, or to inspect compositions:

``` r

tss(counts)             # each column sums to 1
#>           s1         s2        s3        s4
#> t1 0.6666667 0.18181818 0.4285714 0.4000000
#> t2 0.0000000 0.72727273 0.5714286 0.1333333
#> t3 0.3333333 0.09090909 0.0000000 0.4666667
```

It accepts a vector or a matrix and safely maps all-zero columns to zero
rather than `NaN`.

## Functional distances with `traits2dist()`

Functional diversity needs a pairwise **distance matrix** over taxa.
Build one from a trait table with
[`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md),
which defaults to Gower distance so mixed continuous and categorical
traits are handled correctly. Constant traits (which carry no
information) are dropped automatically.

``` r

traits <- data.frame(
  body_mass = c(1.0, 0.2, 0.9),
  diet      = factor(c("herb", "carn", "carn")),
  row.names = c("t1", "t2", "t3")
)
fdist <- traits2dist(traits)
round(fdist, 3)
#>       t1    t2    t3
#> t1 0.000 1.000 0.562
#> t2 1.000 0.000 0.437
#> t3 0.562 0.437 0.000

hilldiv(counts, dist = fdist)
#> Computing functional Hill numbers of "q0", "q1", and
#> "q2".
#> <hilldiv3 result: functional>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 1.353846
#> 2  1     s1 1.343253
#> 3  2     s1 1.333333
#> 4  0     s2 1.911682
#> 5  1     s2 1.658248
#> 6  2     s2 1.517241
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 1.650074
#> 11 1     s4 1.617041
#> 12 2     s4 1.590106
```

## Matching counts to a tree or distances

Phylogenetic and functional analyses require the taxa in the count table
to line up with the tree tips or distance-matrix names. `hilldiv3`
enforces this in one place:

- it checks the two name sets with
  [`setequal()`](https://rdrr.io/r/base/sets.html) (order does not
  matter);
- it **reorders** the counts to match the reference — eliminating the
  silent misalignment that plagued earlier tools;
- if the sets disagree, it stops with a message listing the offending
  taxa.

``` r

bad_tree <- ape::read.tree(text = "((t1:1,t2:1):1,t9:2);")
hilldiv(counts, tree = bad_tree)
#> Error in `.abort_mismatch()`:
#> ! Taxa names in the count data and the tree tips do not match.
#> ℹ Only in counts: "t3"
#> ℹ Only in tree tips: "t9"
```

### Intentionally dropping taxa with `match_data()`

Real datasets rarely match perfectly — a reference tree may include taxa
you never observed, or your table may carry taxa absent from the tree.
[`match_data()`](https://alberdilab.github.io/hilldiv3/reference/match_data.md)
restricts and reorders the count table to the shared taxa, reporting
what it dropped, so you can curate alignment *before* analysis:

``` r

tree <- ape::read.tree(text = "((t1:1,t2:1):1,t4:2);")   # has t4, lacks t3
matched <- match_data(counts, tree = tree)
#> Dropped 1 taxon from `data` not in the tree tips.
#> 1 taxon in the tree tips has no counts; prune it before downstream analysis.
matched
#>    s1 s2 s3 s4
#> t1 10  2  3  6
#> t2  0  8  4  2
```

Feed the aligned table back in — now the name sets agree and the
analysis runs:

``` r

tree3 <- ape::read.tree(text = "(t1:1,t2:1);")
hilldiv(match_data(counts, tree = tree3), tree = tree3)
#> Dropped 1 taxon from `data` not in the tree tips.
#> Computing phylogenetic Hill numbers of "q0", "q1", and "q2".
#> <hilldiv3 result: phylogenetic>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 1.000000
#> 2  1     s1 1.000000
#> 3  2     s1 1.000000
#> 4  0     s2 2.000000
#> 5  1     s2 1.649385
#> 6  2     s2 1.470588
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 2.000000
#> 11 1     s4 1.754765
#> 12 2     s4 1.600000
```

[`match_data()`](https://alberdilab.github.io/hilldiv3/reference/match_data.md)
works the same way with a `dist` argument for functional data.

## A typical workflow

Putting it together:

``` r

# 1. Counts + traits in hand.
traits <- data.frame(
  body_mass = c(1.0, 0.2, 0.9),
  diet      = factor(c("herb", "carn", "carn")),
  row.names = c("t1", "t2", "t3")
)

# 2. Derive a functional distance matrix.
fdist <- traits2dist(traits)

# 3. Align counts to it (here they already match).
counts2 <- match_data(counts, dist = fdist)

# 4. Measure and compare.
hilldiv(counts2, dist = fdist)
#> Computing functional Hill numbers of "q0", "q1", and
#> "q2".
#> <hilldiv3 result: functional>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 1.353846
#> 2  1     s1 1.343253
#> 3  2     s1 1.333333
#> 4  0     s2 1.911682
#> 5  1     s2 1.658248
#> 6  2     s2 1.517241
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 1.650074
#> 11 1     s4 1.617041
#> 12 2     s4 1.590106
hilldiss(counts2, dist = fdist, q = 1)
#> dissimilarity from functional Hill numbers of "q1".
#> <hilldiv3 result: functional>
#> 4 rows x 3 cols
#> 
#>   q metric      value
#> 1 1      S 0.16295408
#> 2 1      C 0.09403070
#> 3 1      U 0.09403070
#> 4 1      V 0.04641062
```
