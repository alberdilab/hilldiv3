# Getting started with hilldiv3

``` r

library(hilldiv3)
```

## What hilldiv3 does

`hilldiv3` measures and compares the diversity of biological communities
(OTU/ASV/MAG tables) using **Hill numbers**, a single family of metrics
that unifies richness, Shannon and Simpson diversity through one
parameter, the *diversity order* `q`. From the same framework it derives
diversity **partitioning**, **(dis)similarity**, **profiles**,
**evenness** and **redundancy**, for three flavours of diversity:

| Flavour          | What it accounts for     | How you ask for it |
|------------------|--------------------------|--------------------|
| **Neutral**      | abundances only          | counts             |
| **Phylogenetic** | evolutionary relatedness | counts + `tree`    |
| **Functional**   | trait dissimilarity      | counts + `dist`    |

The diversity *type* is inferred from the inputs you pass — you call the
same functions either way.

## The data

Every function takes a **count table** with taxa (OTUs/ASVs/MAGs) in
rows and samples in columns. A matrix is the simplest form:

``` r

counts <- matrix(
  c(10, 0, 5,
     2, 8, 1,
     3, 4, 0,
     6, 2, 7),
  nrow = 3, byrow = FALSE,
  dimnames = list(c("t1", "t2", "t3"),
                  c("s1", "s2", "s3", "s4"))
)
counts
#>    s1 s2 s3 s4
#> t1 10  2  3  6
#> t2  0  8  4  2
#> t3  5  1  0  7
```

Data frames, tibbles, `phyloseq` objects and `TreeSummarizedExperiment`
objects work too — see `vignette("preparing-data")` *(article on the
website)*.

## Alpha diversity

[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md)
returns Hill numbers per sample. By default it computes orders `q = 0`
(richness), `q = 1` (Shannon diversity) and `q = 2` (Simpson diversity):

``` r

hilldiv(counts)
#> Computing neutral Hill numbers of "q0", "q1", and
#> "q2".
#>          s1       s2       s3       s4
#> q0 2.000000 3.000000 2.000000 3.000000
#> q1 1.889882 2.137309 1.979626 2.693484
#> q2 1.800000 1.753623 1.960000 2.528090
```

Higher `q` down-weights rare taxa, so `qD` decreases as `q` grows unless
the sample is perfectly even. Add a tree or a distance matrix to switch
flavour:

``` r

tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:2);")
hilldiv(counts, tree = tree)              # phylogenetic
#> Computing phylogenetic Hill numbers of "q0", "q1", and
#> "q2".
#>     s1       s2       s3       s4
#> q0 1.5 2.000000 1.500000 2.000000
#> q1 1.5 1.598520 1.406992 1.677678
#> q2 1.5 1.406977 1.324324 1.500000
```

## Partitioning and dissimilarity

[`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md)
splits diversity across samples into **alpha** (within-sample),
**gamma** (pooled) and **beta** (`gamma / alpha`, the number of
effectively distinct communities):

``` r

hillpart(counts)
#> Partitioning neutral Hill numbers of "q0", "q1", and
#> "q2".
#>       alpha    gamma     beta
#> q0 2.500000 3.000000 1.200000
#> q1 2.154269 2.905735 1.348827
#> q2 1.968927 2.828374 1.436505
```

[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md)
and
[`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md)
turn beta into bounded dissimilarity / similarity metrics (Sorensen-,
Jaccard-, and Unifrac-type), and
[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md)
returns a `dist` object of pairwise dissimilarities ready for
ordination:

``` r

hilldiss(counts, q = 1)
#> dissimilarity from neutral Hill numbers of "q1".
#>            S         C         U         V
#> q1 0.3448199 0.2158525 0.2158525 0.1162756
```

## Where to next

The website carries in-depth articles:

- **Diversity types** — neutral, phylogenetic and functional
  measurement.
- **Partitioning and (dis)similarity** — alpha/beta/gamma, the S/C/U/V
  metrics, and pairwise dissimilarity for ordination.
- **Profiles, evenness and redundancy** —
  [`hillprof()`](https://alberdilab.github.io/hilldiv3/reference/hillprof.md),
  [`hilleven()`](https://alberdilab.github.io/hilldiv3/reference/hilleven.md),
  [`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md).
- **Preparing your data** — input formats,
  `phyloseq`/`TreeSummarizedExperiment`,
  [`tss()`](https://alberdilab.github.io/hilldiv3/reference/tss.md),
  [`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md)
  and
  [`match_data()`](https://alberdilab.github.io/hilldiv3/reference/match_data.md).

## References

- Hill, M.O. (1973). Diversity and evenness. *Ecology*, 54, 427–432.
- Jost, L. (2007). Partitioning diversity into independent alpha and
  beta components. *Ecology*, 88, 2427–2439.
- Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of
  Hill numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*,
  19, 804–817.
