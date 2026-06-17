#' Diversity profile across a range of orders
#'
#' Compute a diversity profile: Hill numbers evaluated over a fine sweep of
#' diversity orders `q`. Profiles are the standard diagnostic for comparing the
#' diversity of assemblages, since the ranking of samples can change with `q`.
#'
#' @inheritParams hilldiv
#' @param q Numeric vector of diversity orders to evaluate. Defaults to a fine
#'   sweep from 0 to 3.
#' @param out Output type: `"tibble"` (default, long format ready for plotting)
#'   or `"matrix"`.
#'
#' @return A long-format `data.frame` of class `hill_profile` (columns `q`,
#'   `sample`, `value`) with a [plot()][plot.hill_profile] method, or a matrix
#'   (orders in rows, samples in columns) when `out = "matrix"`.
#'
#' @seealso [hilldiv()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' prof <- hillprof(counts)
#' plot(prof)
#' @export
hillprof <- function(data, q = seq(0, 3, by = 0.1), tree = NULL, dist = NULL,
                     tau = NULL, out = c("tibble", "matrix")) {
  out <- match.arg(out)
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q)
  type <- attr(x, "type")
  cli::cli_inform("Computing {type} diversity profile over
                   {length(q)} order{?s}.")
  mat <- hill_alpha(x$counts, q = q, type = type,
                    tree = x$tree, dist = x$dist, tau = tau)
  if (out == "matrix") return(mat)

  samples <- colnames(mat)
  df <- data.frame(
    q = rep(q, times = length(samples)),
    sample = rep(samples, each = length(q)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  structure(df, class = c("hill_profile", "data.frame"))
}

#' Plot a diversity profile
#'
#' Base-graphics plot of a [hillprof()] result: one line per sample showing the
#' Hill number against the diversity order `q`.
#'
#' @param x A `hill_profile` object from [hillprof()].
#' @param ... Further arguments passed to [plot()].
#' @return The `hill_profile` object, invisibly.
#' @exportS3Method graphics::plot hill_profile
plot.hill_profile <- function(x, ...) {
  samples <- unique(x$sample)
  cols <- grDevices::hcl.colors(max(length(samples), 2), "Dark 3")
  plot(NA, xlim = range(x$q), ylim = range(x$value),
       xlab = "Diversity order (q)", ylab = "Hill number (qD)", ...)
  for (i in seq_along(samples)) {
    s <- x[x$sample == samples[i], ]
    s <- s[order(s$q), ]
    graphics::lines(s$q, s$value, col = cols[i], lwd = 2)
  }
  graphics::legend("topright", legend = samples, col = cols[seq_along(samples)],
                   lwd = 2, bty = "n")
  invisible(x)
}

#' Hill-number evenness
#'
#' Evenness expressed through Hill numbers as the ratio of diversity of order
#' `q` to richness (`qD / 0D`), which ranges from 0 to 1.
#'
#' @inheritParams hilldiv
#' @param q Numeric vector of diversity orders (> 0 are meaningful for
#'   evenness). Defaults to `c(1, 2)`.
#'
#' @return A matrix of evenness values (orders in rows, samples in columns).
#' @seealso [hilldiv()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hilleven(counts)
#' @export
hilleven <- function(data, q = c(1, 2), tree = NULL, dist = NULL, tau = NULL) {
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q)
  type <- attr(x, "type")
  cli::cli_inform("Computing {type} evenness of {.val {paste0('q', q)}}.")

  qd <- hill_alpha(x$counts, q = q, type = type,
                   tree = x$tree, dist = x$dist, tau = tau)
  q0 <- hill_alpha(x$counts, q = 0, type = type,
                   tree = x$tree, dist = x$dist, tau = tau)
  even <- sweep(qd, 2, q0[1, ], "/")
  rownames(even) <- paste0("q", q)
  even
}
