## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* The check reports a NOTE for "New submission". This is the first
  submission of 'hilldiv3' to CRAN.

* The same NOTE flags possibly misspelled words in DESCRIPTION
  ('ASV', 'OTU', 'Rao', 'Sorensen', 'Unifrac', 'Jost', 'Chao', 'Chiu',
  'Alberdi', 'et', 'al'). These are domain terms, author surnames and
  the abbreviation "et al." in the cited references; they are spelled
  correctly.

## Suggested packages on Bioconductor

* Three suggested packages used by optional input adapters
  ('phyloseq', 'SummarizedExperiment', 'TreeSummarizedExperiment')
  are distributed via Bioconductor. They are listed in 'Suggests' and
  every use is guarded with `rlang::check_installed()` /
  `requireNamespace()`, so the package builds, checks and runs without
  them. The Bioconductor repository is declared in
  'Additional_repositories'. The "Additional_repositories" availability
  line in the incoming-feasibility NOTE refers to this.

## Downstream dependencies

* There are currently no downstream dependencies (new package).

## Test environments

* local macOS, R 4.3.3 -- 0 errors | 0 warnings | 1 note
* win-builder, R-devel (2026-06-26 r90195 ucrt) -- 0 errors | 0 warnings |
  1 note (the "New submission" NOTE described above)
