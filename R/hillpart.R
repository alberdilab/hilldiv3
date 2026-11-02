#' Hill numbers diversity partitioning
#'
#' Partition neutral, phylogenetic or functional Hill-number diversity into
#' alpha, gamma and beta components across a set of samples. With a `hierarchy`
#' formula it instead performs *multi-scale* (nested) partitioning, returning
#' one beta per hierarchical level.
#'
#' @inheritParams hilldiv
#' @param type Diversity type: `"auto"` (default) infers it from the inputs
#'   (counts only -> neutral, `+tree` -> phylogenetic, `+dist` -> functional);
#'   an explicit `"neutral"`, `"phylogenetic"` or `"functional"` asserts the
#'   type and is validated against the inputs (e.g. `"phylogenetic"` requires a
#'   `tree`; `"neutral"` ignores any tree/dist carried by the object).
#' @param data A count table (taxa x samples) or a supported object; a single
#'   sample is not meaningful for partitioning.
#' @param hierarchy Optional one-sided nesting formula, coarsest to finest, e.g.
#'   `~ region / site`, requesting multi-scale (nested) partitioning instead of
#'   the default single-level partition. One beta is returned per hierarchical
#'   transition and the chain telescopes exactly:
#'   `gamma = alpha_finest * prod(beta)`. Works for all three diversity types
#'   (neutral, phylogenetic, functional); see the partitioning vignette for the
#'   shared construction. A single tree `depth` or functional `tau` is used
#'   across every scale in a partition, so `alpha * prod(beta) = gamma` exactly
#'   and every beta is at least one. Grouping variables are resolved against
#'   `metadata` when supplied, otherwise against the calling environment.
#'   Linear unit weighting changes results relative to the older `p_u^q`
#'   formula in unbalanced designs for `q != 1` (a few per cent in neutral
#'   examples; balanced equal-weight designs are identical). With functional
#'   diversity it can also shift `q = 1` slightly (about 0.003 in the gut
#'   data), because a unit's measure-weighted mass need not equal its sample
#'   weight. For
#'   functional diversity, a vector of `tau` values runs the full hierarchy
#'   once per threshold; per-scale thresholds are not supported.
#' @param metadata Optional per-sample `data.frame` supplying the variables
#'   named in `hierarchy`; rows are matched to the count-table columns by name
#'   when possible, otherwise by position.
#' @param weights Sample weights for hierarchical partitioning: `"equal"`
#'   (default), `"abundance"` (proportional to each sample's total count), or
#'   a positive numeric vector with one value per sample. A named vector is
#'   matched to count-table columns; an unnamed vector uses column order.
#'   Weights are normalised to sum to one. Requires `hierarchy`.
#' @param depth Optional positive phylogenetic reference depth for hierarchical
#'   partitioning. On an ultrametric tree, use only branch segments within the
#'   last `depth` units before the tips. Must not exceed the tree height.
#'   `NULL` (default) uses the current shared mean-depth normaliser. Requires
#'   `hierarchy` and a phylogenetic `tree`.
#' @param out Output shape: `"tibble"` (default) returns a long-format
#'   `data.frame` with columns `q`, `component`, `value`; `"matrix"` returns the
#'   legacy matrix (orders in rows, `alpha`/`gamma`/`beta` in columns). With
#'   `hierarchy`, `"tibble"` returns one row per `(q, scale)` and `"matrix"`
#'   returns `alpha`, one `beta_<level>` per nesting level, and `gamma`. A
#'   multi-`tau` hierarchy returns an extra `tau` column for `"tibble"`, or a
#'   named list of matrices for `"matrix"`.
#'
#' @return A long-format `data.frame` of class `hill_partition` (default) with a
#'   `plot()` method, or a matrix with columns `alpha`, `gamma`, `beta` and
#'   diversity orders in rows when `out = "matrix"`. With `hierarchy`, a
#'   `hill_hierarchy` long-format `data.frame` (with its own `plot()` method) or
#'   the corresponding wide matrix.
#'
#' @seealso [hilldiv()], [hilldiss()], [hillsim()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hillpart(counts)
#' plot(hillpart(counts))
#'
#' # Multi-scale partitioning across a nested design.
#' set.seed(1)
#' tab <- matrix(rpois(12 * 8, 5), nrow = 12,
#'               dimnames = list(paste0("t", 1:12), paste0("s", 1:8)))
#' md <- data.frame(region = rep(c("N", "S"), each = 4),
#'                  site = rep(c("a", "b", "c", "d"), each = 2),
#'                  row.names = paste0("s", 1:8))
#' hillpart(tab, hierarchy = ~ region / site, metadata = md)
#' @export
hillpart <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL,
                     hierarchy = NULL, metadata = NULL,
                     type = c("auto", "neutral", "phylogenetic", "functional"),
                     out = c("tibble", "matrix"), weights = "equal",
                     depth = NULL) {
  out <- match.arg(out)
  type <- match.arg(type)
  if (is.null(hierarchy) && !missing(weights)) {
    cli::cli_abort("{.arg weights} requires {.arg hierarchy}.")
  }
  if (is.null(hierarchy) && !missing(depth)) {
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
    levs <- paste(c("sample", spec$level_names, "total"), collapse = " < ")
    cli::cli_inform("Partitioning {type} Hill numbers across scales
                     {.val {levs}}.")
    if (length(tau) > 1L) {
      runs <- lapply(tau, function(tau_value) {
        prep <- hier_part_prep(x$counts, type = type, dist = x$dist,
                               tau = tau_value, sample_weights = sample_weights)
        df <- hier_partition(prep, spec$groupings, q, spec$level_names)
        df$tau <- tau_value
        df
      })
      if (out == "matrix") {
        matrices <- lapply(runs, .hier_to_matrix, q = q)
        names(matrices) <- make.unique(paste0("tau_", tau))
        return(matrices)
      }
      df <- do.call(rbind, runs)
      rownames(df) <- NULL
      return(new_hill_result(df, "hill_hierarchy", type))
    }
    prep <- hier_part_prep(x$counts, type = type, tree = x$tree,
                           dist = x$dist, tau = tau,
                           sample_weights = sample_weights, depth = depth)
    df <- hier_partition(prep, spec$groupings, q, spec$level_names)
    if (out == "matrix") return(.hier_to_matrix(df, q))
    return(new_hill_result(df, "hill_hierarchy", type))
  }

  cli::cli_inform("Partitioning {type} Hill numbers of
                   {.val {paste0('q', q)}}.")
  mat <- hill_partition(x$counts, q = q, type = type,
                        tree = x$tree, dist = x$dist, tau = tau)
  if (out == "matrix") return(mat)
  new_hill_result(.hill_longify(mat, q, "component"), "hill_partition", type)
}
