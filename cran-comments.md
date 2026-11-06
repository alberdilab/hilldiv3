## Update (3.2.0)

This is an update to 'hilldiv3', which is on CRAN as 3.1.1. It adds
multi-scale (hierarchical) Hill partitioning to `hillpart()`, `hilldiss()`,
`hillsim()` and `hillpair()`, and a new `hillcontrib()` for taxon
contributions to dissimilarity. See NEWS.md for details.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps: "unable to verify current time". This
  comes from the check machine having no way to query the time, not from the
  package.

## Downstream dependencies

* There are currently no downstream dependencies.

## Test environments

* local macOS (aarch64), R 4.3.3, `R CMD check --as-cran --no-manual`
