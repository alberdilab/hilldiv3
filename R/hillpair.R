#' Pairwise Hill numbers-based dissimilarity
#'
#' Compute dissimilarity metrics for every pair of samples, returning distance
#' objects suitable for ordination (e.g. NMDS, PCoA).
#'
#' The type-specific structure (per-sample normalisation, the tree traversal or
#' the functional similarity product) is computed **once** over all samples via
#' the partitioning engine; each pair then only combines its two precomputed
#' columns into beta, which is turned into the requested overlap metrics. The
#' maths are therefore identical to [hilldiss()] on two samples, without
#' re-running the full engine per pair. When `parallel = TRUE` and the `furrr`
#' package is installed, pairs are computed in parallel via the active `future`
#' plan. A `progressr` progress bar is reported when that package is installed
#' and a handler is active.
#'
#' @inheritParams hilldiss
#' @param out Output type: `"dist"` (default) returns a `dist` object per
#'   requested metric/order combination; `"tibble"` returns a long-format table.
#' @param parallel Logical; if `TRUE` and `furrr` is installed, compute pairs in
#'   parallel.
#'
#' @return For `out = "dist"`, a named list of `dist` objects (one per
#'   order/metric, named e.g. `"q0S"`), collapsed to a single `dist` when only
#'   one combination is requested. For `out = "tibble"`, a long-format
#'   `data.frame` with columns `first`, `second`, `q`, `metric`, `value`.
#'
#' @seealso [hilldiss()], [hilldiv()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1, 3, 4, 0, 6, 2, 7), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"),
#'                                  c("s1", "s2", "s3", "s4")))
#' hillpair(counts, q = 1, metric = "C")
#' @export
hillpair <- function(data, q = c(0, 1, 2), metric = c("S", "C", "U", "V"),
                     tree = NULL, dist = NULL, tau = NULL,
                     type = c("auto", "neutral", "phylogenetic", "functional"),
                     out = c("dist", "tibble"), parallel = FALSE) {
  out <- match.arg(out)
  type <- match.arg(type)
  metric <- match.arg(metric, c("S", "C", "U", "V"), several.ok = TRUE)
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q, type)
  type <- attr(x, "type")
  counts <- x$counts
  samples <- colnames(counts)
  N <- ncol(counts)
  if (N < 2) {
    cli::cli_abort("Pairwise dissimilarity needs at least two samples.")
  }

  pairs <- utils::combn(N, 2, simplify = FALSE)
  cli::cli_inform("Computing {type} pairwise dissimilarity for
                   {length(pairs)} sample pair{?s}.")

  # Compute the type-specific structure once over all samples; each pair only
  # combines its two columns below.
  prep <- part_prep(counts, type, tree = x$tree, dist = x$dist, tau = tau)

  # beta -> dissimilarity for one pair; returns a (q x metric) matrix.
  pair_fun <- function(idx) {
    betas <- part_eval(prep, idx, q)[, "beta"]
    d <- t(vapply(seq_along(q),
                  function(i) beta_to_dissim(betas[i], 2L, q[i]),
                  c(S = 0, C = 0, U = 0, V = 0)))
    d[, metric, drop = FALSE]
  }

  results <- .pair_lapply(pairs, pair_fun, parallel)

  if (out == "tibble") {
    return(.pairs_to_tibble(results, pairs, samples, q, metric))
  }
  .pairs_to_dist(results, q, metric, samples)
}

# lapply over pairs, optionally parallelised with furrr/future, reporting a
# progressr bar when that package is available and a handler is active.
.pair_lapply <- function(x, fun, parallel) {
  if (rlang::is_installed("progressr")) {
    p <- progressr::progressor(steps = length(x))
    step_fun <- function(idx) {
      on.exit(p())
      fun(idx)
    }
  } else {
    step_fun <- fun
  }
  if (parallel) {
    if (rlang::is_installed("furrr")) {
      return(furrr::future_map(x, step_fun))
    }
    cli::cli_warn("{.pkg furrr} is not installed; running sequentially.")
  }
  lapply(x, step_fun)
}

# Assemble per-pair (q x metric) matrices into a long data.frame.
.pairs_to_tibble <- function(results, pairs, samples, q, metric) {
  first <- samples[vapply(pairs, `[`, integer(1), 1L)]
  second <- samples[vapply(pairs, `[`, integer(1), 2L)]
  rows <- lapply(seq_along(pairs), function(k) {
    data.frame(
      first = first[k], second = second[k],
      q = rep(q, times = length(metric)),
      metric = rep(metric, each = length(q)),
      value = as.vector(results[[k]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Assemble per-pair matrices into a (named list of) dist object(s). combn(N, 2)
# emits pairs in the same order dist stores its lower triangle, so values map
# straight onto the dist body.
.pairs_to_dist <- function(results, q, metric, samples) {
  combos <- expand.grid(q = q, metric = metric, stringsAsFactors = FALSE)
  dists <- lapply(seq_len(nrow(combos)), function(c) {
    qi <- match(combos$q[c], q)
    vals <- vapply(results, function(m) m[qi, combos$metric[c]], numeric(1))
    .as_dist(vals, samples)
  })
  names(dists) <- paste0("q", combos$q, combos$metric)
  if (length(dists) == 1L) dists[[1]] else dists
}

# Build a dist object from its lower-triangle values (in combn order).
.as_dist <- function(vals, labels) {
  structure(vals, Size = length(labels), Labels = labels,
            Diag = FALSE, Upper = FALSE, class = "dist")
}
