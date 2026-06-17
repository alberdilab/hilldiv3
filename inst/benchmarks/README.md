# Performance benchmarks

Run the benchmark suite from the package root:

```r
source("inst/benchmarks/run-benchmarks.R")
```

The script writes `results/performance.csv`, `results/performance-summary.csv`,
`results/performance-times.csv` and `results/session-info.txt`. `performance.csv`
contains one row per benchmark iteration; `performance-summary.csv` contains the
aggregated one-row-per-package-operation view. The pkgdown article
`vignettes/articles/performance.Rmd` reads those files and renders the
benchmark tables and boxplots.

Install the optional comparison packages before the final publication run:

```r
install.packages(c("bench", "hillR", "entropart", "vegan", "BAT"))
remotes::install_github("anttonalberdi/hilldiv2")
```

Use these environment variables to scale the run:

- `BENCH_ITERATIONS`: repetitions per package-operation call.
- `BENCH_N_TAXA`: number of simulated taxa.
- `BENCH_N_SAMPLES`: number of simulated samples.
- `BENCH_OUTPUT_DIR`: output directory for CSV/session files.
- `BENCH_MEMORY_LIMIT_GB`: per-run resident memory cap (default `10`). A run
  exceeding this is killed, marked `out_of_memory` in the summary, and logged as
  `Out of memory ...` in `progress.log`. Enforced by running each call in a
  forked child whose RSS is polled (Unix only; set to `0` to disable). Install
  the optional `ps` package for more accurate memory readings.

Each run is executed in a forked child whose peak resident memory (RSS) is
recorded in the `memory_bytes` column and logged as `peak_memory_bytes` in
`progress.log`, independent of the `BENCH_MEMORY` allocation profiler.

The default benchmark run uses 10 iterations on a synthetic dataset with
200 taxa and 50 samples. The active benchmark set is:

| Benchmark | Operation | Facet |
|---|---|---|
| B1.1 | Alpha diversity at q = 0, 1, 2 | Neutral |
| B1.2 | Alpha diversity at q = 0, 1, 2 | Phylogenetic |
| B1.3 | Alpha diversity at q = 0, 1, 2 | Functional |
| B2.1 | Diversity partitioning at q = 0, 1, 2 | Neutral |
| B2.2 | Diversity partitioning at q = 0, 1, 2 | Phylogenetic |
| B2.3 | Diversity partitioning at q = 0, 1, 2 | Functional |
| B3.1 | Nested diversity partitioning at q = 0, 1, 2 | Neutral |
| B3.2 | Nested diversity partitioning at q = 0, 1, 2 | Phylogenetic |
| B3.3 | Nested diversity partitioning at q = 0, 1, 2 | Functional |
| B4.1 | Pairwise beta diversity matrix at q = 1 | Neutral |
| B4.2 | Pairwise beta diversity matrix at q = 1 | Phylogenetic |
| B4.3 | Pairwise beta diversity matrix at q = 1 | Functional |
