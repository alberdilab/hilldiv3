# Diversity types: neutral, phylogenetic and functional

``` r

library(hilldiv3)
```

`hilldiv3` measures three flavours of diversity from the same count
table. You never switch functions — you switch *inputs*, and the
diversity type is inferred:

| You pass        | Diversity type   |
|-----------------|------------------|
| counts only     | **neutral**      |
| counts + `tree` | **phylogenetic** |
| counts + `dist` | **functional**   |

This article walks through all three with
[`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md),
the alpha-diversity workhorse. The same `tree` / `dist` logic applies to
every other `hill*` function.

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

## The diversity order `q`

Hill numbers are a one-parameter family. The order `q` sets how much
weight rare taxa carry:

- `q = 0` — **richness**: every taxon counts equally, abundance ignored.
- `q = 1` — **Shannon diversity**: taxa weighted by their frequency (the
  exponential of Shannon entropy).
- `q = 2` — **Simpson diversity**: dominated by common taxa (inverse
  Simpson concentration).

``` r

hilldiv(counts, q = c(0, 1, 2))
#> Computing neutral Hill numbers of "q0", "q1", and
#> "q2".
#>          s1       s2       s3       s4
#> q0 2.000000 3.000000 2.000000 3.000000
#> q1 1.889882 2.137309 1.979626 2.693484
#> q2 1.800000 1.753623 1.960000 2.528090
```

All Hill numbers are in **effective number of taxa** (“how many
equally-abundant taxa would give this diversity”), so values across `q`
are directly comparable. A diversity *profile* sweeps `q` continuously —
see the profiles article.

## Neutral diversity

With counts alone you get neutral diversity, which treats every taxon as
equally distinct:

``` r

hilldiv(counts)
#> Computing neutral Hill numbers of "q0", "q1", and
#> "q2".
#>          s1       s2       s3       s4
#> q0 2.000000 3.000000 2.000000 3.000000
#> q1 1.889882 2.137309 1.979626 2.693484
#> q2 1.800000 1.753623 1.960000 2.528090
```

## Phylogenetic diversity

Supply a phylogenetic `tree` (class `phylo`, e.g. from
[`ape::read.tree()`](https://rdrr.io/pkg/ape/man/read.tree.html)) whose
tip labels match the taxa. Diversity is then measured in units of branch
length, crediting samples that span deeper, more divergent lineages
(Chao et al. 2010).

``` r

tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:2);")
hilldiv(counts, tree = tree)
#> Computing phylogenetic Hill numbers of "q0", "q1", and
#> "q2".
#>     s1       s2       s3       s4
#> q0 1.5 2.000000 1.500000 2.000000
#> q1 1.5 1.598520 1.406992 1.677678
#> q2 1.5 1.406977 1.324324 1.500000
```

The tip labels must cover the same taxa as the count table. `hilldiv3`
**reorders** the counts to match the tree internally — a frequent source
of silent error in older tools — and raises an informative error if the
name sets disagree. To intentionally drop taxa that are missing from one
side, align them first with
[`match_data()`](https://alberdilab.github.io/hilldiv3/reference/match_data.md)
(see the data-preparation article).

At `q = 0`, phylogenetic Hill numbers relate directly to Faith’s PD.

## Functional diversity

Supply a `dist` matrix of pairwise functional distances between taxa.
Diversity then accounts for how trait-different the taxa are: a sample
of three very similar taxa is less functionally diverse than three
contrasting ones (Chiu & Chao 2014).

``` r

# A toy functional distance matrix over the three taxa.
fdist <- as.matrix(dist(
  data.frame(body = c(1, 0.2, 0.9), diet = c(0, 1, 1),
             row.names = c("t1", "t2", "t3"))
))
hilldiv(counts, dist = fdist)
#> Computing functional Hill numbers of "q0", "q1", and
#> "q2".
#>          s1       s2       s3       s4
#> q0 1.601908 2.046926 2.000000 1.946876
#> q1 1.566804 1.739361 1.979626 1.909900
#> q2 1.535588 1.569081 1.960000 1.878525
```

### Building distances from traits

Usually you start from a **trait table** (taxa in rows, traits in
columns) and convert it with
[`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md),
which uses Gower distance by default so that mixed
continuous/categorical traits are handled sensibly:

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
```

### The distance threshold `tau`

Functional Hill numbers depend on a threshold `tau`: taxa farther apart
than `tau` are treated as fully distinct. By default `tau = max(dist)`,
which makes the measure use the full range of distances. Lowering `tau`
makes more taxa “functionally equivalent” and reduces functional
diversity:

``` r

hilldiv(counts, dist = fdist, tau = max(fdist))   # default
#> Computing functional Hill numbers of "q0", "q1", and
#> "q2".
#>          s1       s2       s3       s4
#> q0 1.353846 1.911682 2.000000 1.650074
#> q1 1.343253 1.658248 1.979626 1.617041
#> q2 1.333333 1.517241 1.960000 1.590106
hilldiv(counts, dist = fdist, tau = max(fdist) / 2)
#> Computing functional Hill numbers of "q0", "q1", and
#> "q2".
#>          s1       s2       s3       s4
#> q0 2.000000 2.484615 2.000000 2.661169
#> q1 1.889882 1.984284 1.979626 2.524573
#> q2 1.800000 1.704225 1.960000 2.432432
```

## One taxon, many flavours

Because all three share the Hill-number scale, you can compare what each
lens sees in the *same* samples:

``` r

list(
  neutral       = hilldiv(counts),
  phylogenetic  = hilldiv(counts, tree = tree),
  functional    = hilldiv(counts, dist = fdist)
)
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#> Computing phylogenetic Hill numbers of "q0", "q1", and "q2".
#> Computing functional Hill numbers of "q0", "q1", and "q2".
#> $neutral
#>          s1       s2       s3       s4
#> q0 2.000000 3.000000 2.000000 3.000000
#> q1 1.889882 2.137309 1.979626 2.693484
#> q2 1.800000 1.753623 1.960000 2.528090
#> 
#> $phylogenetic
#>     s1       s2       s3       s4
#> q0 1.5 2.000000 1.500000 2.000000
#> q1 1.5 1.598520 1.406992 1.677678
#> q2 1.5 1.406977 1.324324 1.500000
#> 
#> $functional
#>          s1       s2       s3       s4
#> q0 1.353846 1.911682 2.000000 1.650074
#> q1 1.343253 1.658248 1.979626 1.617041
#> q2 1.333333 1.517241 1.960000 1.590106
```

A sample can be neutrally diverse yet phylogenetically or functionally
redundant — quantifying exactly that gap is what
[`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md)
does (see the profiles/evenness/redundancy article).

## References

- Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity
  measures based on Hill numbers. *Phil. Trans. R. Soc. B*, 365,
  3599–3609.
- Chiu, C.-H. & Chao, A. (2014). Distance-based functional diversity
  measures and their decomposition. *PLoS ONE*, 9, e100014.
- Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of
  Hill numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*,
  19, 804–817.
