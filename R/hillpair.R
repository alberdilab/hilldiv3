#' Pairwise Hill numbers-based dissimilarity
#'
#' Compute dissimilarity metrics for every pair of samples, returning distance
#' objects suitable for ordination (e.g. NMDS, PCoA).
#'
#' Unlike hilldiv2, the planned implementation reuses precomputed branch /
#' distance structure across pairs and can run in parallel via the `future`
#' framework, instead of re-running a full partition for each pair.
#'
#' @inheritParams hilldiss
#' @param out Output type: `"dist"` (default) returns a `dist` object per
#'   requested metric/order; `"tibble"` returns a long-format table.
#' @param parallel Logical; if `TRUE` and the `future`/`furrr` packages are
#'   installed, compute pairs in parallel.
#'
#' @return A `dist` object, a named list of `dist` objects, or a long-format
#'   table (see `out`).
#'
#' @seealso [hilldiss()], [hilldiv()]
#' @export
hillpair <- function(data, q = c(0, 1, 2), metric = c("S", "C", "U", "V"),
                     tree = NULL, dist = NULL, tau = NULL,
                     out = c("dist", "tibble"), parallel = FALSE) {
  out <- match.arg(out)
  # TODO(engine): vectorised pairwise computation reusing branch_abundance()
  # / distance structure; optional furrr backend; assemble dist/tibble output.
  .NotYetImplemented()
}
