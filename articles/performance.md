# Performance benchmarks

The benchmark suite compares `hilldiv3` with the legacy `hilldiv2` code
and the external packages `hillR`, `entropart`, `vegan` and `BAT`.

The goal is not to force every package into the same table when the
operation is not equivalent. Each package-operation row is labelled as:

| Status | Meaning |
|----|----|
| `yes` | The package has a directly comparable operation. |
| `partial` | The package has a related operation, but the output or statistical definition differs. |
| `no` | The package does not expose that operation. |
| `unavailable` | The package was not installed when the benchmark was run. |
| `failed` | The package was installed, but the benchmark call errored. |

## Reproducing the benchmark

Install the comparison packages, then run the benchmark script from the
package root:

``` r

install.packages(c("bench", "hillR", "entropart", "vegan", "BAT"))
remotes::install_github("anttonalberdi/hilldiv2")

Sys.setenv(
  BENCH_ITERATIONS = 5,
  BENCH_N_TAXA = 200,
  BENCH_N_SAMPLES = 40
)
source("inst/benchmarks/run-benchmarks.R")
```

The script writes:

| File | Contents |
|----|----|
| `inst/benchmarks/results/performance.csv` | Support status, timing and memory table. |
| `inst/benchmarks/results/session-info.txt` | R version, platform and package versions. |

Increase `BENCH_ITERATIONS`, `BENCH_N_TAXA` or `BENCH_N_SAMPLES` for the
final publication run. Keep `session-info.txt` with the published
results because benchmark times depend on hardware, BLAS, R version and
package versions.

No benchmark results have been generated yet. Run
`inst/benchmarks/run-benchmarks.R` and rebuild the pkgdown site to
populate this page.

## Support matrix

## Timing and memory

Rows marked `partial` are timed because they are useful context, but
they should not be interpreted as exact replacements for the
corresponding `hilldiv3` operation.

## Notes on equivalence

`hillR` is the closest external comparator for Hill-number alpha
diversity, partitioning and pairwise comparisons across taxonomic,
phylogenetic and functional diversity.

`entropart` supports metacommunity alpha, beta and gamma diversity,
including phylogenetic and similarity-based diversity, but some outputs
and assumptions differ from the `hilldiv3` API.

`vegan` provides neutral Hill numbers through Renyi diversity and many
pairwise community dissimilarities, but it does not implement the
phylogenetic, functional or Hill-number S/C/U/V dissimilarity operations
used by `hilldiv3`.

`BAT` includes biodiversity assessment tools for taxonomic, phylogenetic
and functional diversity. Its Hill-number function is a neutral
alpha-diversity comparator; its beta-diversity tools are related but not
the same Hill-number partition/dissimilarity operations.
