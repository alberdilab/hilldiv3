## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* The check reports a NOTE for "New submission". This is the first
  submission of 'hilldiv3' to CRAN.

## Suggested packages on Bioconductor

* Three suggested packages used by optional input adapters
  ('phyloseq', 'SummarizedExperiment', 'TreeSummarizedExperiment')
  are distributed via Bioconductor. They are listed in 'Suggests' and
  every use is guarded with `rlang::check_installed()` /
  `requireNamespace()`, so the package builds, checks and runs without
  them. The Bioconductor repository is declared in
  'Additional_repositories'.

## Downstream dependencies

* There are currently no downstream dependencies (new package).

## Test environments

* local macOS, R 4.3.3
* win-builder (devel and release)  [to be run before submission]
* R-hub (Windows, macOS, Linux)    [to be run before submission]
