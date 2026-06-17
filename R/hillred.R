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
#' @return A `data.frame` of class `hill_redundancy` (default), or a matrix with
#'   columns `redundancy`, `a`, `b`, `c` (one row per `q`) when
#'   `out = "matrix"`.
#'
#' @seealso [hilldiv()]
#' @export
hillred <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL,
                    type = c("auto", "phylogenetic", "functional"),
                    out = c("tibble", "matrix")) {
  out <- match.arg(out)
  type <- match.arg(type)
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q, type)
  type <- attr(x, "type")
  if (type == "neutral") {
    cli::cli_abort("Redundancy requires either a {.arg tree} or a {.arg dist}.")
  }

  relationship <- function(x, a, b, c) -a * 2^(-x / b) + c

  x_all <- hill_alpha(x$counts, q = q, type = "neutral")
  y_all <- hill_alpha(x$counts, q = q, type = type,
                      tree = x$tree, dist = x$dist, tau = tau)

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
  new_hill_result(df, "hill_redundancy", type, value_label = "Redundancy")
}
