#' Diversity profile across a range of orders
#'
#' Compute a diversity profile: Hill numbers evaluated over a fine sweep of
#' diversity orders `q`. Profiles are the standard diagnostic for comparing the
#' diversity of assemblages, since the ranking of samples can change with `q`.
#'
#' @inheritParams hilldiv
#' @param q Numeric vector of diversity orders to evaluate. Defaults to a fine
#'   sweep from 0 to 3.
#' @param out Output type: `"tibble"` (default, long format ready for ggplot2)
#'   or `"matrix"`.
#'
#' @return A long-format table (or matrix) of diversity by `q` and sample.
#'   Objects of class `hill_profile` gain a `plot()` method.
#'
#' @seealso [hilldiv()]
#' @export
hillprof <- function(data, q = seq(0, 3, by = 0.1), tree = NULL, dist = NULL,
                     tau = NULL, out = c("tibble", "matrix")) {
  out <- match.arg(out)
  # TODO(feature): call hill_alpha() across `q`, reshape to long, wrap in a
  # `hill_profile` S3 object, and add plot()/autoplot() methods.
  .NotYetImplemented()
}

#' Hill-number evenness
#'
#' Evenness expressed through Hill numbers, e.g. the ratio of diversity of order
#' `q` to richness (`qD / 0D`), which ranges from 0 to 1.
#'
#' @inheritParams hilldiv
#' @param q Numeric vector of diversity orders (> 0 are meaningful for
#'   evenness). Defaults to `c(1, 2)`.
#'
#' @return A matrix of evenness values (orders in rows, samples in columns).
#' @seealso [hilldiv()]
#' @export
hilleven <- function(data, q = c(1, 2), tree = NULL, dist = NULL, tau = NULL) {
  # TODO(feature): compute hill_alpha() at q and at q = 0, return qD / 0D
  # (and optionally other evenness normalisations).
  .NotYetImplemented()
}
