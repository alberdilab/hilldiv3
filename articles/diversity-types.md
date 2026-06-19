# Diversity types: neutral, phylogenetic and functional

``` r

library(hilldiv3)
```

“Diversity” can mean three different things, and which one you want
depends on the question you are asking:

- **Neutral** diversity counts taxa and weighs them by abundance,
  treating every taxon as equally distinct. Two communities with the
  same number of equally common taxa are equally diverse, no matter
  *which* taxa they are.
- **Phylogenetic** diversity also asks how *related* the taxa are. A
  community of three close cousins is less diverse than three taxa from
  distant branches of the tree of life, even if both have three taxa.
- **Functional** diversity asks how *different* the taxa are in their
  traits (size, diet, metabolism…). Three taxa that do very similar
  things are less diverse than three that do contrasting things.

`hilldiv3` measures all three from the same count table. You never
switch functions — you switch *inputs*, and the diversity type is
inferred from what you supply:

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
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> <hilldiv3 result: neutral>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 2.000000
#> 2  1     s1 1.889882
#> 3  2     s1 1.800000
#> 4  0     s2 3.000000
#> 5  1     s2 2.137309
#> 6  2     s2 1.753623
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 3.000000
#> 11 1     s4 2.693484
#> 12 2     s4 2.528090
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
#> Computing neutral Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> <hilldiv3 result: neutral>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 2.000000
#> 2  1     s1 1.889882
#> 3  2     s1 1.800000
#> 4  0     s2 3.000000
#> 5  1     s2 2.137309
#> 6  2     s2 1.753623
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 3.000000
#> 11 1     s4 2.693484
#> 12 2     s4 2.528090
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
#> Computing phylogenetic Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> <hilldiv3 result: phylogenetic>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 2.000000
#> 2  1     s1 1.889882
#> 3  2     s1 1.800000
#> 4  0     s2 2.500000
#> 5  1     s2 1.702490
#> 6  2     s2 1.423529
#> 7  0     s3 1.500000
#> 8  1     s3 1.406992
#> 9  2     s3 1.324324
#> 10 0     s4 2.500000
#> 11 1     s4 2.318405
#> 12 2     s4 2.227723
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
#> Computing functional Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> <hilldiv3 result: functional>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 1.601908
#> 2  1     s1 1.566804
#> 3  2     s1 1.535588
#> 4  0     s2 2.046926
#> 5  1     s2 1.739361
#> 6  2     s2 1.569081
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 1.946876
#> 11 1     s4 1.909900
#> 12 2     s4 1.878525
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
#> Computing functional Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
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
hilldiv(counts, dist = fdist, tau = max(fdist) / 2)
#> Computing functional Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> <hilldiv3 result: functional>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 2.000000
#> 2  1     s1 1.889882
#> 3  2     s1 1.800000
#> 4  0     s2 2.484615
#> 5  1     s2 1.984284
#> 6  2     s2 1.704225
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 2.661169
#> 11 1     s4 2.524573
#> 12 2     s4 2.432432
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
#> ℹ 3 taxa across 4 samples.
#> Computing phylogenetic Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> Computing functional Hill numbers of "q0", "q1", and "q2".
#> ℹ 3 taxa across 4 samples.
#> $neutral
#> <hilldiv3 result: neutral>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 2.000000
#> 2  1     s1 1.889882
#> 3  2     s1 1.800000
#> 4  0     s2 3.000000
#> 5  1     s2 2.137309
#> 6  2     s2 1.753623
#> 7  0     s3 2.000000
#> 8  1     s3 1.979626
#> 9  2     s3 1.960000
#> 10 0     s4 3.000000
#> 11 1     s4 2.693484
#> 12 2     s4 2.528090
#> 
#> $phylogenetic
#> <hilldiv3 result: phylogenetic>
#> 12 rows x 3 cols
#> 
#>    q sample    value
#> 1  0     s1 2.000000
#> 2  1     s1 1.889882
#> 3  2     s1 1.800000
#> 4  0     s2 2.500000
#> 5  1     s2 1.702490
#> 6  2     s2 1.423529
#> 7  0     s3 1.500000
#> 8  1     s3 1.406992
#> 9  2     s3 1.324324
#> 10 0     s4 2.500000
#> 11 1     s4 2.318405
#> 12 2     s4 2.227723
#> 
#> $functional
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

A sample can be neutrally diverse yet phylogenetically or functionally
redundant — quantifying exactly that gap is what
[`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md)
does (see the profiles/evenness/redundancy article).

## Contrasts with other packages

> **Advanced — skip on a first read.** This section is for users coming
> from other Hill-number tools who want to know why the numbers can
> differ. If you are just getting started, you can move on to the next
> article.

Other R packages compute Hill numbers, and their **phylogenetic**
results can disagree with `hilldiv3` even on identical data. These are
not bugs: they are different conventions for the *reference depth* `T` —
the tree depth at which an assemblage is read when branch lengths are
converted into effective numbers of lineages. Knowing them makes
migration painless. (Neutral and functional results match across
packages; only the phylogenetic depth convention varies.)

On an **ultrametric** tree every sample shares the same depth, so all
conventions coincide. The differences below surface only on
**non-ultrametric** trees (e.g. most genome or gene trees), and grow
with how uneven the abundances are.

### `hillR::hill_phylo`

Two adjustments line it up with `hilldiv3`:

- By default `hill_phylo()` returns the raw phylogenetic diversity `PD`
  (branch-length units; Faith’s PD at `q = 0`), **not** an effective
  number of taxa. Pass `return_dt = TRUE` and read the `D_t` row to get
  the Hill number `PD / T`.
- It reads each sample at its *own* depth `T_j = sum(L_i a_i)`, so it
  corresponds to `reference = "sample"`, not `reference = "pool"`.

``` r

# hillR expects sites x species; hilldiv3 expects taxa x samples.
hillR::hill_phylo(t(counts), tree, q = 1, return_dt = TRUE)["D_t", ]
hilldiv(counts, tree = tree, q = 1, reference = "sample")
```

With those two changes the packages agree **exactly for `q >= 1`**. They
differ only at `q = 0`: `hillR` rebuilds `T` from a presence/absence
(uniform) community for the richness case, whereas `hilldiv3` keeps a
single abundance-based `T` across the whole `q` profile, following Chao
et al. (2010). The Faith’s PD numerator is identical; only the `T` in
the denominator differs.

### `hilldiv2::hilldiv`

`reference = "sample"` reproduces `hilldiv2`’s **single-vector**
phylogenetic formula exactly. `hilldiv2`’s **matrix** path is different:
it folds a shared even-pool depth *into* the diversity transform. That
is fine on ultrametric trees but breaks down on non-ultrametric ones —
the branch weights no longer sum to one, and the `q`-profile becomes
discontinuous around `q = 1`. `hilldiv3`’s `reference = "pool"` delivers
comparable, common-depth values without that pathology, because it
applies the reference depth as a final rescaling rather than inside the
exponentiated sum.

| Want to match | Use in `hilldiv3` |
|----|----|
| `hillR::hill_phylo(..., return_dt = TRUE)` | `reference = "sample"` (exact for `q >= 1`) |
| `hilldiv2::hilldiv(vector, tree)` | `reference = "sample"` (exact) |
| `hilldiv2::hilldiv(matrix, tree)` | no exact match — `reference = "pool"` is the sound analogue |

In short: when results diverge, check (1) whether the other package
reports `PD` or an effective number, and (2) which reference depth it
uses. `hilldiv3` makes that choice explicit through `reference`.

## References

- Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity
  measures based on Hill numbers. *Phil. Trans. R. Soc. B*, 365,
  3599–3609.
- Chiu, C.-H. & Chao, A. (2014). Distance-based functional diversity
  measures and their decomposition. *PLoS ONE*, 9, e100014.
- Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of
  Hill numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*,
  19, 804–817.
