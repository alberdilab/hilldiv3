# Contributing to hilldiv3

Thank you for your interest in `hilldiv3`. Bug reports, questions,
feature ideas, documentation fixes and code contributions are all
welcome.

## Getting help

- **Questions about using the package** — first check the [documentation
  website](https://alberdilab.github.io/hilldiv3/), in particular the
  *Get started* vignette and the articles. If your question is not
  answered there, open an
  [issue](https://github.com/alberdilab/hilldiv3/issues) and tag it as a
  question.
- **Questions about Hill numbers themselves** (which `q` to use, how to
  interpret beta, etc.) — Alberdi & Gilbert (2019),
  <doi:10.1111/1755-0998.13014>, is a practical guide; issues are also
  fine for these.

## Reporting bugs

Open an
[issue](https://github.com/alberdilab/hilldiv3/issues/new/choose) and
include:

1.  A **minimal reproducible example** — ideally using the bundled data
    (`gut_counts`, `gut_tree`, `gut_traits`) or a small simulated table.
    The [reprex](https://reprex.tidyverse.org/) package makes this easy.
2.  What you expected to happen and what happened instead (full error
    message).
3.  The output of
    [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) (or
    [`sessioninfo::session_info()`](https://sessioninfo.r-lib.org/reference/session_info.html)).

If you suspect a numerical error, please say which reference value you
compared against (a hand calculation, another package, or a published
formula), so the discrepancy can be traced to either a bug or a
convention difference. Known convention differences with other packages
are documented in the *Diversity types* and *Partitioning and
(dis)similarity* articles.

## Suggesting features

Open an issue describing the analysis you want to do, why the current
functions do not cover it, and — where relevant — the methodological
reference that defines the measure. Please discuss larger features in an
issue before writing code, so we can agree on the design first.

## Contributing code or documentation

1.  Fork the repository and create a branch from `main`.
2.  Install the development dependencies:
    `devtools::install_dev_deps()`.
3.  Make your change.
    - Diversity maths belongs in the internal engine (`R/engine-*.R`);
      the user-facing `hill*()` functions validate inputs and call the
      engine.
    - Add or update tests in `tests/testthat/`. New measures should have
      a test against an independent reference value (see
      `test-golden.R`).
    - Document exported functions with roxygen2 and run
      `devtools::document()`.
    - Add a bullet to the top section of `NEWS.md` describing the
      change.
4.  Check that everything passes locally: `devtools::check()` and
    [`lintr::lint_package()`](https://lintr.r-lib.org/reference/lint.html).
5.  Open a pull request against `main` and describe what changed and
    why. Continuous integration runs `R CMD check` on Linux, macOS and
    Windows, the linter and test coverage.

Small fixes (typos, broken links, clarifications) can be submitted
directly as a pull request without opening an issue first.

## Support expectations

`hilldiv3` is maintained by the Alberdi Lab at the University of
Copenhagen. We aim to respond to new issues and pull requests within two
weeks. Bug reports affecting the correctness of results are given
priority. The package follows [semantic
versioning](https://semver.org/); breaking changes are announced in
`NEWS.md`.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://alberdilab.github.io/hilldiv3/CODE_OF_CONDUCT.md). By
participating in this project you agree to abide by its terms.
