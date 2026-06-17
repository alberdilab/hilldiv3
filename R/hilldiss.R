#' Hill numbers-based dissimilarity
#'
#' Compute overall (multi-sample) dissimilarity metrics from the Hill-number
#' beta diversity following Chiu et al. (2014). These are the complements of the
#' similarities returned by [hillsim()].
#'
#' @inheritParams hillpart
#' @param metric Dissimilarity metric(s) to return, any of `"S"`, `"C"`, `"U"`,
#'   `"V"`. Defaults to all four.
#' @param out Output shape: `"tibble"` (default) returns a long-format
#'   `data.frame` with columns `q`, `metric`, `value`; `"matrix"` returns the
#'   legacy matrix (orders in rows, metrics in columns, dropped to a vector for
#'   a single metric).
#'
#' @return A long-format `data.frame` of class `hill_dissimilarity` (default,
#'   with a `plot()` method), or a matrix/vector of dissimilarities when
#'   `out = "matrix"`.
#'
#' @seealso [hillsim()], [hillpair()], [hillpart()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hilldiss(counts)
#' plot(hilldiss(counts))
#' @export
hilldiss <- function(data, q = c(0, 1, 2), metric = c("S", "C", "U", "V"),
                     tree = NULL, dist = NULL, tau = NULL,
                     type = c("auto", "neutral", "phylogenetic", "functional"),
                     out = c("tibble", "matrix")) {
  .hill_overlap(data, q, metric, tree, dist, tau, kind = "dissimilarity",
                type = match.arg(type), out = match.arg(out))
}

#' Hill numbers-based similarity
#'
#' Compute overall similarity metrics from the Hill-number beta diversity
#' (Chiu et al. 2014). These are `1 -` the dissimilarities from [hilldiss()].
#'
#' @inheritParams hilldiss
#' @return A long-format `data.frame` of class `hill_similarity` (default, with
#'   a `plot()` method), or a matrix/vector of similarities when
#'   `out = "matrix"`.
#' @seealso [hilldiss()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hillsim(counts)
#' plot(hillsim(counts))
#' @export
hillsim <- function(data, q = c(0, 1, 2), metric = c("S", "C", "U", "V"),
                    tree = NULL, dist = NULL, tau = NULL,
                    type = c("auto", "neutral", "phylogenetic", "functional"),
                    out = c("tibble", "matrix")) {
  .hill_overlap(data, q, metric, tree, dist, tau, kind = "similarity",
                type = match.arg(type), out = match.arg(out))
}

# Shared implementation for hilldiss()/hillsim().
.hill_overlap <- function(data, q, metric, tree, dist, tau, kind, type, out) {
  metric <- match.arg(metric, c("S", "C", "U", "V"), several.ok = TRUE)
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q, type)
  type <- attr(x, "type")
  N <- ncol(x$counts)
  cli::cli_inform("{kind} from {type} Hill numbers of {.val {paste0('q', q)}}.")

  part <- hill_partition(x$counts, q = q, type = type,
                         tree = x$tree, dist = x$dist, tau = tau)
  betas <- part[, "beta"]

  fun <- if (kind == "similarity") beta_to_sim else beta_to_dissim
  res <- t(vapply(seq_along(q),
                  function(i) fun(betas[i], N, q[i]),
                  c(S = 0, C = 0, U = 0, V = 0)))
  rownames(res) <- paste0("q", q)
  res <- res[, metric, drop = FALSE]
  if (out == "matrix") {
    return(res[, , drop = length(metric) == 1])
  }
  subclass <- if (kind == "similarity") {
    "hill_similarity"
  } else {
    "hill_dissimilarity"
  }
  new_hill_result(.hill_longify(res, q, "metric"), subclass, type)
}
