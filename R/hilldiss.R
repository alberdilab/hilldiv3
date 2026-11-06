#' Hill numbers-based dissimilarity
#'
#' Compute overall (multi-sample) dissimilarity metrics from the Hill-number
#' beta diversity following Chiu et al. (2014). These are the complements of the
#' similarities returned by [hillsim()]. With `hierarchy`, each transition's
#' beta is converted using its data-dependent effective ceiling rather than
#' the number of child units, keeping the metrics in `[0, 1]`.
#'
#' @inheritParams hillpart
#' @param hierarchy Optional one-sided nesting formula, coarsest to finest,
#'   e.g. `~ region / site`. Every transition above `sample`, including
#'   `total`, gets its own dissimilarity metrics. Each run uses one shared tree
#'   `depth` or functional `tau`; a functional `tau` grid runs the full
#'   hierarchy once per threshold.
#' @param metric Dissimilarity metric(s) to return, any of `"S"`, `"C"`, `"U"`,
#'   `"V"`. Defaults to all four. `"S"` (Jaccard-type) and `"V"` (Sorensen-type)
#'   are turnover measures, relative to the pooled assemblage and to a single
#'   assemblage respectively; `"U"` (Jaccard-type) and `"C"` (Sorensen-type;
#'   one minus the Morisita-Horn similarity at `q = 2`) are the corresponding
#'   overlap complements (proportions of non-shared units). At `q = 0`,
#'   `"U"` equals `"S"` and `"C"` equals `"V"`; at `q = 1`, `"U"` equals `"C"`.
#' @param out Output shape: `"tibble"` (default) returns a long-format
#'   `data.frame` with columns `q`, `metric`, `value`; `"matrix"` returns the
#'   legacy matrix (orders in rows, metrics in columns, dropped to a vector for
#'   a single metric). With `hierarchy`, the table adds a `scale` column for
#'   each transition above the sample scale. A `tau` grid adds a `tau` column.
#'   Matrix output has one row per `q` and transition, with a named list of
#'   matrices for a `tau` grid.
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
                     out = c("tibble", "matrix"), hierarchy = NULL,
                     metadata = NULL, weights = "equal", depth = NULL) {
  .hill_overlap(data, q, metric, tree, dist, tau, kind = "dissimilarity",
                type = match.arg(type), out = match.arg(out),
                hierarchy = hierarchy, metadata = metadata, weights = weights,
                depth = depth, weights_supplied = !missing(weights),
                depth_supplied = !missing(depth))
}

#' Hill numbers-based similarity
#'
#' Compute overall similarity metrics from the Hill-number beta diversity
#' (Chiu et al. 2014). These are `1 -` the dissimilarities from [hilldiss()].
#'
#' @inheritParams hilldiss
#' @param metric Similarity metric(s) to return, any of `"S"`, `"C"`, `"U"`,
#'   `"V"`. Defaults to all four. Each is one minus the dissimilarity of the
#'   same code in [hilldiss()], so identical assemblages give 1. `"C"` at
#'   `q = 2` is the Morisita-Horn similarity.
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
                    out = c("tibble", "matrix"), hierarchy = NULL,
                    metadata = NULL, weights = "equal", depth = NULL) {
  .hill_overlap(data, q, metric, tree, dist, tau, kind = "similarity",
                type = match.arg(type), out = match.arg(out),
                hierarchy = hierarchy, metadata = metadata, weights = weights,
                depth = depth, weights_supplied = !missing(weights),
                depth_supplied = !missing(depth))
}

# Shared implementation for hilldiss()/hillsim().
.hill_overlap <- function(data, q, metric, tree, dist, tau, kind, type, out,
                          hierarchy, metadata, weights, depth,
                          weights_supplied, depth_supplied) {
  metric <- match.arg(metric, c("S", "C", "U", "V"), several.ok = TRUE)
  if (is.null(hierarchy) && weights_supplied) {
    cli::cli_abort("{.arg weights} requires {.arg hierarchy}.")
  }
  if (is.null(hierarchy) && depth_supplied) {
    cli::cli_abort("{.arg depth} requires {.arg hierarchy}.")
  }
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q, type)
  type <- attr(x, "type")
  if (!is.null(hierarchy)) {
    if (!is.null(depth) && type != "phylogenetic") {
      cli::cli_abort("{.arg depth} requires phylogenetic diversity and a
                      {.arg tree}.")
    }
    if (length(tau) > 1L && type != "functional") {
      cli::cli_abort("A {.arg tau} grid requires functional diversity and a
                      {.arg dist} matrix.")
    }
    sample_weights <- .hier_sample_weights(weights, x$counts)
    if (is.null(metadata)) metadata <- x$metadata
    spec <- .parse_hierarchy(hierarchy, metadata, colnames(x$counts))
    tau_values <- if (length(tau) > 1L) as.list(tau) else list(tau)
    fun <- if (kind == "similarity") beta_to_sim else beta_to_dissim
    runs <- lapply(tau_values, function(tau_value) {
      prep <- hier_part_prep(x$counts, type = type, tree = x$tree,
                             dist = x$dist, tau = tau_value,
                             sample_weights = sample_weights, depth = depth)
      part <- hier_partition(prep, spec$groupings, q, spec$level_names,
                             with_ceiling = TRUE)
      part <- part[!is.na(part$beta), ]
      res <- t(vapply(seq_len(nrow(part)), function(i) {
        if (part$n_eff[i] <= 1 + 1e-12) {
          if (kind == "similarity") {
            return(c(S = 1, C = 1, U = 1, V = 1))
          }
          return(c(S = 0, C = 0, U = 0, V = 0))
        }
        fun(part$beta[i], part$n_eff[i], part$q[i])
      }, c(S = 0, C = 0, U = 0, V = 0)))
      res <- res[, metric, drop = FALSE]
      if (out == "matrix") {
        rownames(res) <- paste0("q", part$q, "_", part$scale)
        return(res)
      }
      df <- data.frame(q = rep(part$q, times = length(metric)),
                       scale = rep(part$scale, times = length(metric)),
                       metric = rep(metric, each = nrow(part)),
                       value = as.vector(res), stringsAsFactors = FALSE)
      if (length(tau) > 1L) df$tau <- tau_value
      df
    })
    if (out == "matrix") {
      if (length(runs) == 1L) return(runs[[1]])
      names(runs) <- make.unique(paste0("tau_", tau))
      return(runs)
    }
    df <- do.call(rbind, runs)
    rownames(df) <- NULL
    if ("tau" %in% names(df)) {
      df <- df[, c("q", "scale", "tau", "metric", "value")]
    }
    subclass <- if (kind == "similarity") "hill_similarity" else
      "hill_dissimilarity"
    return(new_hill_result(df, subclass, type))
  }
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
