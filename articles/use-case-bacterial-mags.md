# Case 2: Three faces of diversity in a gut microbiome

``` r

library(hilldiv3)
library(ggplot2)

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "top"))

# hilldiv3 results are long-format data.frame subclasses; drop the class so
# base `$<-` / `rbind` (and ggplot2) work without coercion.
as_df <- function(x) { class(x) <- "data.frame"; x }
pal <- c(control = "#1F77B4", treatment = "#D62728")
```

> **New to Hill numbers?** This is a complete worked example. If terms
> like diversity order `q`, alpha/beta diversity or
> phylogenetic/functional diversity are unfamiliar, start with the
> [Getting
> started](https://alberdilab.github.io/hilldiv3/articles/hilldiv3.md)
> guide and the [Diversity
> types](https://alberdilab.github.io/hilldiv3/articles/diversity-types.md)
> article, then come back here.

Genome-resolved metagenomics summarises a microbial community as a table
of metagenome-assembled genome (MAG) abundances. Each MAG carries two
extra layers of information that a raw count ignores: its **position on
the bacterial phylogeny** and its **functional attributes** (genome
size, GC content, oxygen tolerance, encoded capabilities). These layers
can move independently: an intervention may reshuffle *which* organisms
— and *which lineages* — dominate a community while leaving *what the
community does* untouched, if the incoming genomes are functionally
equivalent to the ones they replace. `hilldiv3` measures all three
flavours from the same count table and the same calls, so taxonomic,
phylogenetic and functional change can be compared on a common,
interpretable scale (effective numbers of lineages) — and, crucially,
told apart.

## The data

We use the bundled simulated data set: **24 MAGs across 12 host gut
samples**, split into a `control` and a `treatment` group of six. The
MAGs fall into two deep bacterial clades, and the intervention **swaps
which clade dominates**: control guts are dominated by clade A,
treatment guts by clade B, with taxon richness left untouched. The two
clades are *functional mirrors* of one another — for every genome in one
there is a genome in the other with the same trait profile — so a swap
that is dramatic phylogenetically is invisible functionally. `gut_tree`
is the genome phylogeny and `gut_traits` a trait table mixing continuous
(`genome_size`, `gc_content`), categorical (`oxygen`) and binary
(`motility`) attributes.

``` r

group <- rep(c("control", "treatment"), each = 6)
names(group) <- colnames(gut_counts)
metadata <- data.frame(group = group, row.names = colnames(gut_counts))
grp <- function(s) group[s]

dim(gut_counts)                 # 24 MAGs x 12 samples
#> [1] 24 12
range(colSums(gut_counts > 0))  # per-sample richness (near-constant)
#> [1] 23 24
head(gut_traits)
#>       genome_size gc_content         oxygen motility
#> mag01        1.64      0.343       anaerobe        0
#> mag02        1.97      0.354 microaerophile        1
#> mag03        2.34      0.400    facultative        0
#> mag04        2.75      0.420         aerobe        1
#> mag05        3.20      0.447       anaerobe        0
#> mag06        3.66      0.496 microaerophile        1
```

[`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md)
turns that mixed trait table into a functional distance with Gower’s
coefficient, which handles the different variable types automatically:

``` r

fdist <- traits2dist(gut_traits)
range(fdist)
#> [1] 0.0000000 0.9988839
```

## Neutral, phylogenetic and functional diversity from one call

The diversity *type* follows from what you supply, cumulatively: counts
give neutral, `tree =` adds phylogenetic, `dist =` adds functional.
Supplying both a tree and a distance matrix returns all three at once in
a single tibble, with a `type` column telling them apart (`type =`
restricts the output). The interface never changes.

``` r

# One call with both `tree =` and `dist =` returns neutral, phylogenetic and
# functional diversity together, tagged by a `type` column.
labs <- c(neutral = "Neutral (qD)", phylogenetic = "Phylogenetic (qPD)",
          functional = "Functional (qFD)")
alpha <- as_df(hilldiv(gut_counts, q = c(0, 1, 2),
                       tree = gut_tree, dist = fdist))
alpha$flavour <- factor(labs[alpha$type], labs)
alpha$group   <- grp(alpha$sample)

aggregate(value ~ flavour + q + group, alpha, function(x) round(mean(x), 2))
#>               flavour q     group value
#> 1        Neutral (qD) 0   control 23.83
#> 2  Phylogenetic (qPD) 0   control  3.76
#> 3    Functional (qFD) 0   control  1.93
#> 4        Neutral (qD) 1   control 12.25
#> 5  Phylogenetic (qPD) 1   control  2.10
#> 6    Functional (qFD) 1   control  1.90
#> 7        Neutral (qD) 2   control  8.11
#> 8  Phylogenetic (qPD) 2   control  1.64
#> 9    Functional (qFD) 2   control  1.87
#> 10       Neutral (qD) 0 treatment 23.50
#> 11 Phylogenetic (qPD) 0 treatment  3.75
#> 12   Functional (qFD) 0 treatment  1.87
#> 13       Neutral (qD) 1 treatment 11.26
#> 14 Phylogenetic (qPD) 1 treatment  1.96
#> 15   Functional (qFD) 1 treatment  1.85
#> 16       Neutral (qD) 2 treatment  7.78
#> 17 Phylogenetic (qPD) 2 treatment  1.54
#> 18   Functional (qFD) 2 treatment  1.82
```

``` r

ggplot(alpha, aes(factor(q), value, fill = group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  facet_wrap(~ flavour, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "Diversity order (q)", y = "Effective number of lineages")
```

![Neutral, phylogenetic and functional alpha diversity by group across
q](use-case-bacterial-mags_files/figure-html/alpha-plot-1.png)

At the **alpha** (within-sample) level the two groups look almost
identical in every flavour: a single host carries about the same number
of effective MAGs, spread over about the same amount of the phylogeny
and the same functional space, whether or not it received the treatment.
Read at the alpha level alone, the intervention looks like it did
nothing. The effect is not in *how much* diversity each gut holds, but
in *which* lineages hold it — a compositional change that only
between-sample analysis can see.

## Which MAGs carry the power sum?

[`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md)
shows how much of the Hill power sum is assigned to one MAG or a set of
MAGs in each sample. Here we inspect `mag01` and the first phylogenetic
clade in one control sample. Supplying both the tree and functional
distances returns neutral, phylogenetic and functional shares together.

``` r

clade_A <- paste0("mag", sprintf("%02d", 1:12))
shares <- hillshare(gut_counts[, "ctrl01"], q = c(0, 1, 2),
                    tree = gut_tree, dist = fdist,
                    sets = list(mag01 = "mag01", clade_A = clade_A))
shown <- shares
shown$share <- round(shown$share, 3)
shown
#>    q  sample         type     set share allocation
#> 1  0 sample1      neutral   mag01 0.042     qpower
#> 2  0 sample1      neutral clade_A 0.500     qpower
#> 3  1 sample1      neutral   mag01 0.075     qpower
#> 4  1 sample1      neutral clade_A 0.846     qpower
#> 5  2 sample1      neutral   mag01 0.045     qpower
#> 6  2 sample1      neutral clade_A 0.967     qpower
#> 7  0 sample1 phylogenetic   mag01 0.079     qpower
#> 8  0 sample1 phylogenetic clade_A 0.500     qpower
#> 9  1 sample1 phylogenetic   mag01 0.075     qpower
#> 10 1 sample1 phylogenetic clade_A 0.846     qpower
#> 11 2 sample1 phylogenetic   mag01 0.037     qpower
#> 12 2 sample1 phylogenetic clade_A 0.968     qpower
#> 13 0 sample1   functional   mag01 0.070     qpower
#> 14 0 sample1   functional clade_A 0.844     qpower
#> 15 1 sample1   functional   mag01 0.075     qpower
#> 16 1 sample1   functional clade_A 0.846     qpower
#> 17 2 sample1   functional   mag01 0.078     qpower
#> 18 2 sample1   functional clade_A 0.851     qpower
```

``` r

plot_shares <- shares
plot_shares$set <- factor(plot_shares$set, levels = c("mag01", "clade_A"))
plot_shares$type <- factor(plot_shares$type,
                           levels = c("neutral", "phylogenetic", "functional"))
ggplot(plot_shares, aes(factor(q), share, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  facet_wrap(~ set, scales = "free_y",
             labeller = labeller(set = c(mag01 = "MAG 01",
                                         clade_A = "Clade A (12 MAGs)"))) +
  scale_fill_manual(values = c(neutral = "#4E79A7",
                               phylogenetic = "#59A14F",
                               functional = "#E15759"),
                    breaks = c("neutral", "phylogenetic", "functional"),
                    labels = c("Neutral", "Phylogenetic", "Functional"),
                    name = NULL) +
  labs(x = "Diversity order (q)", y = "Share of Hill power sum",
       caption = "Vertical scales differ between panels to show the MAG trend.")
```

![Grouped bars show the shares of mag01 and clade A in control sample
ctrl01 at q 0, 1 and 2, for neutral, phylogenetic and functional
diversity. The two panels use different vertical
scales.](use-case-bacterial-mags_files/figure-html/shares-plot-1.png)

Each share is between zero and one within its sample, type and order.
The `clade_A` share adds the allocations of its 12 MAGs, including
`mag01`; the two rows therefore overlap and should not be added
together. These are shares of the underlying power sum, **not**
percentages of the effective numbers shown above. The functional
allocation accounts for trait similarity among MAGs. In this sample,
neutral and phylogenetic shares of clade A rise with `q`, while its
functional share changes little.

## Where does the community turn over — and where doesn’t it?

[`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md)
partitions diversity between the two groups and reports beta, the
effective number of distinct communities. Running it for each flavour
shows *which* axis the treatment shifts:

``` r

part <- function(type, ...) {
  m <- as_df(hillpart(gut_counts, q = c(0, 1, 2), hierarchy = ~ group,
                      metadata = metadata, ...))
  m <- m[m$scale == "total", c("q", "beta")]
  m$flavour <- type
  m
}
beta <- rbind(part("Neutral"),
              part("Phylogenetic", tree = gut_tree),
              part("Functional",   dist = fdist))
beta$flavour <- factor(beta$flavour, c("Neutral", "Phylogenetic", "Functional"))

ggplot(beta, aes(factor(q), beta, fill = flavour)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = "Diversity order (q)",
       y = expression("Between-group turnover (" * beta * ")"))
```

![Between-group turnover by flavour across
q](use-case-bacterial-mags_files/figure-html/partition-1.png)

Between-group beta grows strongly with `q` for the **neutral and
phylogenetic** decompositions but stays essentially flat (β ≈ 1) for the
**functional** one. The clade swap replaces the dominant genomes with
phylogenetically distant ones, so taxonomic *and* evolutionary
composition turn over almost in lock-step — yet because the two clades
are functional mirrors, the *functional* make-up of the community barely
moves. Distinguishing taxonomic, phylogenetic and functional turnover on
one common beta scale is a core strength of the Hill-number framework as
implemented here, and here it isolates a change that is phylogenetic but
not functional.

## Which MAGs account for the dissimilarity?

The beta results tell us that samples differ;
[`hillcontrib()`](https://alberdilab.github.io/hilldiv3/reference/hillcontrib.md)
shows which MAGs account for a particular dissimilarity. Here we compare
one control and one treatment sample at `q = 2`, which gives abundant
MAGs more weight. The rows are ordered from largest to smallest
contribution.

``` r

pair_counts <- gut_counts[, c("ctrl01", "trt01")]
pair_contrib <- hillcontrib(pair_counts, q = 2, metric = "C")
head(pair_contrib[, c("set", "contribution", "share", "rank")], 5)
#>     set contribution      share rank
#> 1 mag05   0.29972486 0.38273484    1
#> 2 mag14   0.18242450 0.23294768    2
#> 3 mag23   0.09216861 0.11769506    3
#> 4 mag07   0.04400636 0.05619409    4
#> 5 mag18   0.04384198 0.05598419    5

# With one row per MAG, the contributions add up to the pair's dissimilarity.
c(attributed = sum(pair_contrib$contribution),
  pairwise = as.numeric(hillpair(pair_counts, q = 2, metric = "C")))
#> attributed   pairwise 
#>  0.7831136  0.7831136
```

The figure shows the eight largest shares. Colour indicates which sample
has the higher *relative abundance* of each MAG.

``` r

top_contrib <- head(pair_contrib, 8)
pair_relative <- sweep(pair_counts, 2, colSums(pair_counts), "/")
top_contrib$higher_in <- ifelse(
  pair_relative[top_contrib$set, "ctrl01"] >
    pair_relative[top_contrib$set, "trt01"], "ctrl01", "trt01")

ggplot(top_contrib, aes(reorder(set, share), share, fill = higher_in)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c(ctrl01 = pal[["control"]],
                               trt01 = pal[["treatment"]]),
                    labels = c(ctrl01 = "Control (ctrl01)",
                               trt01 = "Treatment (trt01)"),
                    name = "Higher relative abundance") +
  labs(x = NULL, y = "Share of pairwise dissimilarity (q = 2, C)")
```

![Horizontal bars show the eight MAGs contributing most to dissimilarity
between ctrl01 and trt01. mag05 contributes 38 percent and mag14
contributes 23 percent. Bar colours indicate the sample with higher
relative
abundance.](use-case-bacterial-mags_files/figure-html/contrib-plot-1.png)

`mag05` and `mag14` rank highest: `mag05` has 725 counts in `ctrl01` but
35 in `trt01`, while `mag14` has 29 and 690, respectively. Their shares
are about 38% and 23% of this pair’s total dissimilarity; the eight
shown account for about 91% together.

We can also examine **collective turnover within each group**. This
compares the six control samples with one another and, separately, the
six treatment samples with one another; it does not compare the two
groups against each other.

``` r

sample_groups <- split(names(group), group)
group_contrib <- hillcontrib(gut_counts, q = 2, metric = "C",
                            by = "collective", groups = sample_groups)
top_three <- do.call(rbind, lapply(split(group_contrib, group_contrib$group),
                                   function(x) head(x, 3)))
top_three[, c("group", "set", "contribution", "share")]
#>                  group   set contribution     share
#> control.1      control mag12   0.10110789 0.3104164
#> control.2      control mag05   0.09652995 0.2963615
#> control.3      control mag02   0.06068756 0.1863199
#> treatment.25 treatment mag19   0.15466253 0.3405158
#> treatment.26 treatment mag23   0.14308196 0.3150192
#> treatment.27 treatment mag14   0.04910499 0.1081130
```

## Ordination in three spaces

[`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md)
accepts the same `tree =` / `dist =` switch, so one workflow produces
ordinations in *neutral*, *phylogenetic* and *functional* space:

``` r

ord_one <- function(d, lab) {
  pc <- cmdscale(d, k = 2)
  o  <- data.frame(pc, mag = rownames(pc))
  names(o)[1:2] <- c("PCoA1", "PCoA2")
  o$group <- grp(o$mag); o$space <- lab
  o
}
ord <- rbind(
  ord_one(hillpair(gut_counts, q = 1, metric = "C"),                 "Neutral (q = 1)"),
  ord_one(hillpair(gut_counts, q = 1, metric = "C", tree = gut_tree), "Phylogenetic (q = 1)"),
  ord_one(hillpair(gut_counts, q = 1, metric = "C", dist = fdist),   "Functional (q = 1)")
)
ord$space <- factor(ord$space,
                    c("Neutral (q = 1)", "Phylogenetic (q = 1)", "Functional (q = 1)"))

ggplot(ord, aes(PCoA1, PCoA2, colour = group)) +
  geom_point(size = 2.4) +
  stat_ellipse(level = 0.68) +
  facet_wrap(~ space, scales = "free") +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "PCoA 1", y = "PCoA 2")
```

![PCoA in neutral, phylogenetic and functional space, coloured by
group](use-case-bacterial-mags_files/figure-html/ordination-1.png)

The same samples separate cleanly into control and treatment in neutral
space and in phylogenetic space, but collapse onto a single overlapping
cloud in functional space. The picture the three panels paint together —
*who* and *which lineage* differ, *what they do* does not — is one no
single distance could give.

## Functional redundancy

[`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md)
quantifies the **functional redundancy** that underlies this pattern: it
fits the saturating relationship between neutral diversity (number of
genomes) and functional diversity (number of distinct trait profiles)
across samples. A curve that plateaus well below the spread of neutral
diversity means many genomes share functions — so the community can
replace genomes without losing function.

``` r

red <- hillred(gut_counts, q = c(1, 2), dist = fdist)
as_df(red)[, c("q", "redundancy")]
#>   q redundancy
#> 1 1  0.6957595
#> 2 2  0.6024721
plot(red)
```

![Per-sample neutral versus functional diversity with fitted saturating
curves](use-case-bacterial-mags_files/figure-html/redundancy-1.png)

The fitted redundancy is high (around 0.6–0.7): functional diversity
saturates well below the spread of neutral diversity, so adding genomes
contributes far less than proportional new function. This is exactly
*why* the treatment could overturn the community’s taxonomic and
phylogenetic composition without touching its function — every
functional role lost with clade A is recovered from its mirror in clade
B. The redundancy that
[`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md)
measures is the mechanism the beta and ordination analyses revealed.

## Summary

From one MAG table plus a genome tree and a trait table, `hilldiv3`
produced a complete three-flavour analysis through one interface:
trait-to-distance conversion
([`traits2dist()`](https://alberdilab.github.io/hilldiv3/reference/traits2dist.md));
neutral, phylogenetic and functional alpha diversity
([`hilldiv()`](https://alberdilab.github.io/hilldiv3/reference/hilldiv.md));
MAG and clade power-sum shares
([`hillshare()`](https://alberdilab.github.io/hilldiv3/reference/hillshare.md));
between-group partitioning for each flavour
([`hillpart()`](https://alberdilab.github.io/hilldiv3/reference/hillpart.md));
MAG contributions to pairwise and within-group dissimilarity
([`hillcontrib()`](https://alberdilab.github.io/hilldiv3/reference/hillcontrib.md));
ordinations in neutral, phylogenetic and functional space
([`hillpair()`](https://alberdilab.github.io/hilldiv3/reference/hillpair.md));
and functional redundancy
([`hillred()`](https://alberdilab.github.io/hilldiv3/reference/hillred.md)).
The unified, type-switching interface is what lets a single study
separate a change that is taxonomic and phylogenetic from one that is
functional — here revealing a community reorganised in identity but
conserved in function.

## References

- Jost, L. (2007). Partitioning diversity into independent alpha and
  beta components. *Ecology*, 88, 2427–2439.
- Chiu, C.-H., Jost, L. & Chao, A. (2014). Phylogenetic beta diversity,
  similarity, and differentiation measures based on Hill numbers.
  *Ecological Monographs*, 84, 21–44.
- Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of
  Hill numbers to DNA-based diversity analyses. *Mol. Ecol. Resour.*,
  19, 804–817. \`\`\`
