# Partitioning and (dis)similarity

``` r

library(hilldiv3)
```

Measuring *how different* samples are from one another is a separate
question from how diverse each one is. `hilldiv3` answers it by
**partitioning** the total (gamma) diversity of a set of samples into a
within-sample (alpha) and a between-sample (beta) component, then
converting beta into bounded (dis)similarity metrics. All of this works
for neutral, phylogenetic and functional diversity — just add a `tree`
or `dist` as in the [diversity-types
article](https://alberdilab.github.io/hilldiv3/articles/diversity-types.md).

``` r

counts <- matrix(
  c(10, 0, 5,
     2, 8, 1,
     3, 4, 0,
     6, 2, 7),
  nrow = 3,
  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3", "s4"))
)
```

## Alpha, gamma and beta

[`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md)
returns the three components per diversity order:

``` r

hillpart(counts)
#> Partitioning neutral Hill numbers of "q0", "q1", and
#> "q2".
#>       alpha    gamma     beta
#> q0 2.500000 3.000000 1.200000
#> q1 2.154269 2.905735 1.348827
#> q2 1.968927 2.828374 1.436505
```

- **alpha** — the effective diversity of a typical single sample.
- **gamma** — the effective diversity of all samples pooled.
- **beta** — `gamma / alpha`, the effective number of *distinct*
  communities. It ranges from `1` (all samples identical) to `N` (all
  samples completely distinct), where `N` is the number of samples.

Beta is the multiplicative decomposition of Jost (2007): it is
independent of alpha, which is what makes it comparable across datasets
with different richness.

## From beta to (dis)similarity

Beta lives on `[1, N]`, which is awkward to compare across studies with
different sample counts.
[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md)
and
[`hillsim()`](https://alberdilab.github.io/hilldiv3/reference/hillsim.md)
rescale it to `[0, 1]` using the four overlap measures of Chiu et
al. (2014):

| Code | Measure | Type |
|----|----|----|
| `S` | Sorensen-type | local overlap, occurrence-flavoured |
| `C` | Jaccard-type (Morisita-Horn at `q = 2`) | local overlap, abundance-flavoured |
| `U` | one-minus regional overlap | regional |
| `V` | homogeneity / Sorensen complement | regional |

``` r

hilldiss(counts)            # dissimilarities, one column per metric
#> dissimilarity from neutral Hill numbers of "q0", "q1",
#> and "q2".
#>            S          C         U          V
#> q0 0.2222222 0.06666667 0.2222222 0.06666667
#> q1 0.3448199 0.21585249 0.2158525 0.11627556
#> q2 0.4051548 0.40515478 0.1455017 0.14550174
hillsim(counts)             # similarities = 1 - dissimilarity
#> similarity from neutral Hill numbers of "q0", "q1", and
#> "q2".
#>            S         C         U         V
#> q0 0.7777778 0.9333333 0.7777778 0.9333333
#> q1 0.6551801 0.7841475 0.7841475 0.8837244
#> q2 0.5948452 0.5948452 0.8544983 0.8544983
```

Request a single metric and order to get a plain matrix slice:

``` r

hilldiss(counts, q = 1, metric = "C")
#> dissimilarity from neutral Hill numbers of "q1".
#> [1] 0.2158525
```

The choice of `q` matters as much as the metric: `q = 0` weights
presence/ absence (turnover of rare taxa), while `q = 2` reflects shifts
among dominant taxa.

## Pairwise dissimilarity for ordination

[`hilldiss()`](https://alberdilab.github.io/hilldiv3/reference/hilldiss.md)
summarises *all* samples at once. To compare samples two-by-two — the
input ordination methods expect — use
[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md),
which returns a `dist` object (or a long tibble):

``` r

d <- hillpair(counts, q = 1, metric = "C")
#> Computing neutral pairwise dissimilarity for 6 sample pairs.
d
#>            s1         s2         s3
#> s2 0.52298484                      
#> s3 0.47119926 0.08924381           
#> s4 0.09902103 0.29281463 0.33948330
```

That `dist` drops straight into ordination and clustering:

``` r

pcoa <- cmdscale(d, k = 2)
plot(pcoa, type = "n", xlab = "PCoA 1", ylab = "PCoA 2")
text(pcoa, labels = rownames(pcoa))
```

![PCoA of pairwise Hill-number
dissimilarities](partitioning-and-dissimilarity_files/figure-html/ordination-1.png)

Ask for several orders/metrics at once and you get a **named list** of
`dist` objects (e.g. `q0S`, `q1C`):

``` r

ds <- hillpair(counts, q = c(0, 2), metric = c("C", "U"))
#> Computing neutral pairwise dissimilarity for 6 sample pairs.
names(ds)
#> [1] "q0C" "q2C" "q0U" "q2U"
```

Or a tidy long table, convenient for `ggplot2` or `dplyr`:

``` r

head(hillpair(counts, q = 1, metric = "C", out = "tibble"))
#> Computing neutral pairwise dissimilarity for 6 sample pairs.
#>   first second q metric      value
#> 1    s1     s2 1      C 0.52298484
#> 2    s1     s3 1      C 0.47119926
#> 3    s1     s4 1      C 0.09902103
#> 4    s2     s3 1      C 0.08924381
#> 5    s2     s4 1      C 0.29281463
#> 6    s3     s4 1      C 0.33948330
```

### Larger datasets

[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md)
computes every pair through the same partitioning engine, so cost grows
with the number of pairs. For many samples, run pairs in parallel by
setting `parallel = TRUE` with a `future` plan and the `furrr` package
installed:

``` r

library(future)
plan(multisession)
hillpair(counts, q = 1, metric = "C", parallel = TRUE)
```

## References

- Jost, L. (2007). Partitioning diversity into independent alpha and
  beta components. *Ecology*, 88, 2427–2439.
- Chiu, C.-H., Jost, L. & Chao, A. (2014). Phylogenetic beta diversity,
  similarity, and differentiation measures based on Hill numbers.
  *Ecological Monographs*, 84, 21–44.
