#' Hill numbers redundancy
#'
#' Estimate phylogenetic or functional redundancy by fitting the saturating
#' relationship between neutral diversity and phylogenetic/functional diversity
#' across samples: `y = -a * 2^(-x / b) + c`. Redundancy is summarised as
#' `1 - b / max(x)`.
#'
#' @inheritParams hilldiv
#' @param data A count table (taxa x samples); requires either `tree` or `dist`.
#' @param out Output shape: `"tibble"` (default) returns a `data.frame` with one
#'   row per `q` and columns `q`, `redundancy`, `a`, `b`, `c`; `"matrix"`
#'   returns the legacy matrix (orders in rows).
#'
#' @return A `data.frame` of class `hill_redundancy` (default) with a
#'   [plot()][plot.hill_redundancy] method, or a matrix with columns
#'   `redundancy`, `a`, `b`, `c` (one row per `q`) when `out = "matrix"`. The
#'   tibble carries the per-sample neutral and phylogenetic/functional diversity
#'   used for the fit as a `"hill_fit"` attribute, which the plot method draws.
#'
#' @seealso [hilldiv()], [plot.hill_redundancy()]
#' @examples
#' d <- traits2dist(gut_traits)
#' red <- hillred(gut_counts, dist = d)
#' red
#' plot(red)
#' @export
hillred <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL,
                    type = c("auto", "phylogenetic", "functional"),
                    reference = c("pool", "sample"),
                    out = c("tibble", "matrix")) {
  out <- match.arg(out)
  type <- match.arg(type)
  reference <- match.arg(reference)
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q, type)
  type <- attr(x, "type")
  if (type == "neutral") {
    cli::cli_abort("Redundancy requires either a {.arg tree} or a {.arg dist}.")
  }

  relationship <- function(x, a, b, c) -a * 2^(-x / b) + c

  x_all <- hill_alpha(x$counts, q = q, type = "neutral")
  y_all <- hill_alpha(x$counts, q = q, type = type,
                      tree = x$tree, dist = x$dist, tau = tau,
                      reference = reference)

  res <- matrix(NA_real_, nrow = length(q), ncol = 4,
                dimnames = list(paste0("q", q), c("redundancy", "a", "b", "c")))
  for (i in seq_along(q)) {
    xi <- x_all[i, ]
    yi <- y_all[i, ]
    fit <- tryCatch(
      stats::nls(yi ~ relationship(xi, a, b, c),
                 start = list(a = max(yi) - min(yi),
                              b = (max(xi) - min(xi)) / 2,
                              c = max(yi)),
                 control = list(maxiter = 1000)),
      error = function(e) {
        cli::cli_warn("Redundancy for {.val {paste0('q', q[i])}} could not be \\
                       estimated: {conditionMessage(e)}")
        NULL
      }
    )
    if (is.null(fit)) next
    co <- stats::coef(fit)
    res[i, ] <- c(1 - co[["b"]] / max(xi), co[["a"]], co[["b"]], co[["c"]])
  }
  if (out == "matrix") return(res)
  df <- data.frame(q = q, redundancy = res[, "redundancy"], a = res[, "a"],
                   b = res[, "b"], c = res[, "c"], row.names = NULL,
                   stringsAsFactors = FALSE)
  df <- new_hill_result(df, "hill_redundancy", type, value_label = "Redundancy")
  # Retain the points the curve was fitted through (neutral x vs
  # phylogenetic/functional y, per sample and order) so plot()/autoplot() can
  # redraw the saturating relationship, not just the redundancy summary.
  samples <- colnames(x_all)
  attr(df, "hill_fit") <- data.frame(
    q = rep(q, times = length(samples)),
    sample = rep(samples, each = length(q)),
    neutral = as.vector(x_all),
    diversity = as.vector(y_all),
    stringsAsFactors = FALSE
  )
  df
}

#' Plot a redundancy fit
#'
#' Base-graphics plot of a [hillred()] result. For each diversity order `q` it
#' shows the per-sample neutral diversity (x) against phylogenetic/functional
#' diversity (y), overlaid with the fitted saturating curve
#' `y = -a * 2^(-x / b) + c`. A curve that bends sharply and plateaus well below
#' the points' spread indicates high redundancy; a near-linear fit indicates
#' low redundancy. This mirrors the profile plot of [hillprof()].
#'
#' @param x A `hill_redundancy` object from [hillred()].
#' @param ... Further arguments passed to [plot()].
#' @return The `hill_redundancy` object, invisibly.
#' @seealso [hillred()], [plot.hill_profile()]
#' @exportS3Method graphics::plot hill_redundancy
plot.hill_redundancy <- function(x, ...) {
  fit <- attr(x, "hill_fit")
  if (is.null(fit)) {
    # No fit points retained: fall back to the redundancy-vs-q summary.
    x <- x[order(x$q), ]
    plot(x$q, x$redundancy, type = "b", pch = 16, lwd = 2,
         xlab = "Diversity order (q)", ylab = "Redundancy", ...)
    return(invisible(x))
  }
  qs <- sort(unique(fit$q))
  cols <- grDevices::hcl.colors(max(length(qs), 2), "Dark 3")
  plot(NA, xlim = range(fit$neutral, na.rm = TRUE),
       ylim = range(fit$diversity, na.rm = TRUE),
       xlab = "Neutral diversity (qD)",
       ylab = .hill_red_ylab(attr(x, "hill_type")), ...)
  for (i in seq_along(qs)) {
    pts <- fit[fit$q == qs[i], ]
    graphics::points(pts$neutral, pts$diversity, col = cols[i], pch = 16)
    co <- x[x$q == qs[i], c("a", "b", "c")]
    if (all(is.finite(unlist(co)))) {
      xx <- seq(min(pts$neutral), max(pts$neutral), length.out = 100)
      graphics::lines(xx, -co$a * 2^(-xx / co$b) + co$c, col = cols[i], lwd = 2)
    }
  }
  graphics::legend("bottomright", legend = paste0("q", qs),
                   col = cols[seq_along(qs)], lwd = 2, pch = 16, bty = "n")
  invisible(x)
}

# y-axis label for redundancy fits, named after the diversity type when known.
.hill_red_ylab <- function(type) {
  label <- if (is.null(type)) "Phylogenetic/functional" else
    paste0(toupper(substring(type, 1, 1)), substring(type, 2))
  sprintf("%s diversity (qD)", label)
}
