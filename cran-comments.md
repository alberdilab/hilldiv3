## Update (3.2.0)

This is an update to 'hilldiv3', which was published on CRAN as 3.0.0 on
6 October 2026. It adds multi-scale (hierarchical) Hill partitioning to
`hillpart()`, `hilldiss()`, `hillsim()` and `hillpair()`, a new
`hillcontrib()` for taxon contributions to dissimilarity, and `hillshare()`.
It also fixes the README links to repository files that are not part of the
built package, which had caused the incoming checks to report possibly
invalid file URIs for an earlier 3.1.x submission. See NEWS.md for details.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps: "unable to verify current time". This
  comes from the check machine having no way to query the time, not from the
  package.

## Downstream dependencies

* There are currently no downstream dependencies.

## Test environments

* local macOS (aarch64), R 4.3.3, `R CMD check --as-cran --no-manual`
