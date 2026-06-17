# Performance benchmarks

Run the benchmark suite from the package root:

```r
source("inst/benchmarks/run-benchmarks.R")
```

The script writes `results/performance.csv` and `results/session-info.txt`.
The pkgdown article `vignettes/articles/performance.Rmd` reads those files and
renders the benchmark tables.

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
