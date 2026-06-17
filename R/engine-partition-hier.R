#' Multi-scale (hierarchical) neutral Hill-number partitioning engine
#'
#' Internal compute engine for nested, multi-level diversity partitioning. Given
#' a sampling hierarchy (e.g. samples nested in sites nested in regions), it
#' returns one beta per hierarchical transition such that the chain telescopes
#' exactly:
#'
#' \deqn{{}^qD_\gamma = {}^qD_{\alpha} \cdot \prod_k \beta_k}
#'
#' i.e. gamma equals the finest alpha times the product of the per-level betas.
#'
#' The construction uses **global equal sample weights** (every sample weighted
#' \eqn{1/n}). At an aggregation scale \eqn{k} whose units are \eqn{u} (each a
#' pool of \eqn{n_u} samples), the Jost-consistent diversity is
#'
#' \deqn{A_k = \frac{\left(\sum_{u,i} (m_{u,i}/n)^q\right)^{1/(1-q)}}
#'                  {\left(\sum_{u} (n_u/n)^q\right)^{1/(1-q)}}}
#'
#' where \eqn{m_{u,i}} is the summed relative abundance of taxon \eqn{i} over
#' the samples in \eqn{u}. This is exactly Jost's (2007) weighted alpha (the
#' denominator is the numbers-equivalent of the unit weights), so each
#' \eqn{\beta_k = A_k / A_{k-1}} is a proper Hill beta lying in
#' \eqn{[1, n_{\text{units below}}]} and independent of \eqn{A}. At the finest
#' scale (each sample its own unit) \eqn{A_0} reduces to the per-sample alpha of
#' [hill_partition()]; at the coarsest scale (one pooled unit) \eqn{A_K} reduces
#' to gamma. The intermediate scales interpolate, giving the nested betas.
#'
#' @param pi Numeric matrix of relative abundances (taxa x samples), already
#'   tss-normalised so each column sums to one.
#' @param groupings Ordered list of grouping vectors, finest to coarsest, one
#'   entry per *named* hierarchical level (the implicit per-sample and pooled
#'   scales are added here). Each vector has length `ncol(pi)`.
#' @param q Numeric vector of diversity orders.
#' @param level_names Character names for the grouping levels, aligned with
#'   `groupings` (finest to coarsest).
#'
#' @return A long-format `data.frame` with columns `q`, `scale` (ordered factor,
#'   finest to coarsest), `n_units`, `diversity` (\eqn{A_k}) and `beta`
#'   (\eqn{A_k/A_{k-1}}; `NA` at the finest scale).
#' @keywords internal
#' @noRd
hier_partition_neutral <- function(pi, groupings, q, level_names) {
  pi <- as.matrix(pi)
  n <- ncol(pi)

  # Full ordered list of scales, finest -> coarsest: each sample on its own,
  # then the named grouping levels, then everything pooled into one unit.
  sample_g <- seq_len(n)
  total_g <- rep(1L, n)
  scales <- c(list(sample_g), groupings, list(total_g))
  scale_names <- c("sample", level_names, "total")

  # Per-scale, sample-independent set-up: the taxa x units pooled-abundance
  # matrix M (m_{u,i}) and the unit sizes n_u. Done once per scale, reused
  # across q.
  prep <- lapply(scales, function(g) {
    g <- as.factor(g)
    units <- levels(g)
    M <- vapply(units, function(u) rowSums(pi[, g == u, drop = FALSE]),
                numeric(nrow(pi)))
    list(M = as.vector(M / n), nu = tabulate(g) / n, n_units = length(units))
  })

  n_scales <- length(prep)
  out <- vector("list", length(q) * n_scales)
  row <- 1L
  for (qi in seq_along(q)) {
    qv <- q[qi]
    A <- vapply(prep, function(p) {
      .hill_from_props(p$M, qv) / .hill_from_props(p$nu, qv)
    }, numeric(1))
    beta <- c(NA_real_, A[-1] / A[-n_scales])
    for (s in seq_len(n_scales)) {
      out[[row]] <- data.frame(
        q = qv, scale = scale_names[s], n_units = prep[[s]]$n_units,
        diversity = A[s], beta = beta[s], stringsAsFactors = FALSE
      )
      row <- row + 1L
    }
  }
  df <- do.call(rbind, out)
  df$scale <- factor(df$scale, levels = scale_names)
  rownames(df) <- NULL
  df
}

#' Resolve a hierarchy formula into ordered nested grouping vectors
#'
#' Turns a nesting formula such as `~ region / site` into the ordered list of
#' grouping vectors (finest to coarsest) consumed by [hier_partition_neutral()].
#' Variables are resolved against `metadata` (a per-sample data.frame) when
#' supplied, otherwise against the formula's environment. To honour nesting,
#' each named level is the interaction of its variable with all coarser ones,
#' so reused labels (a `site` value appearing in two regions) stay distinct.
#'
#' @param hierarchy A one-sided nesting formula, coarsest to finest, e.g.
#'   `~ region / site`.
#' @param metadata Optional data.frame with one row per sample; rows are aligned
#'   to `sample_ids` by name when both are available.
#' @param sample_ids Column names of the count table (sample identifiers).
#'
#' @return A list with `groupings` (finest-to-coarsest list of factors) and
#'   `level_names` (their names).
#' @keywords internal
#' @noRd
.parse_hierarchy <- function(hierarchy, metadata, sample_ids) {
  if (!inherits(hierarchy, "formula")) {
    cli::cli_abort("{.arg hierarchy} must be a formula, e.g.
                    {.code ~ region / site}.")
  }
  vars <- all.vars(hierarchy)                       # coarsest to finest
  if (length(vars) < 1) {
    cli::cli_abort("{.arg hierarchy} must name at least one grouping variable.")
  }
  n <- length(sample_ids)

  resolve <- function(v) {
    if (!is.null(metadata)) {
      md <- as.data.frame(metadata)
      if (!is.null(rownames(md)) && !is.null(sample_ids) &&
          setequal(intersect(rownames(md), sample_ids), sample_ids)) {
        md <- md[sample_ids, , drop = FALSE]          # align to sample order
      }
      if (!v %in% names(md)) {
        cli::cli_abort("Variable {.field {v}} not found in {.arg metadata}.")
      }
      val <- md[[v]]
    } else {
      val <- eval(as.name(v), envir = environment(hierarchy))
    }
    if (length(val) != n) {
      cli::cli_abort("Hierarchy variable {.field {v}} has length
                      {length(val)} but there are {n} samples.")
    }
    val
  }

  vals <- lapply(vars, resolve)

  # Cumulative interactions, coarsest -> finest: the j-th named level is the
  # interaction of variables 1..j, so finer levels are nested within coarser.
  cum <- vector("list", length(vars))
  for (j in seq_along(vars)) {
    cum[[j]] <- interaction(vals[seq_len(j)], drop = TRUE, sep = ":")
  }

  list(groupings = rev(cum), level_names = rev(vars))
}

# Pivot the long hierarchy frame to a wide matrix (q in rows; columns alpha,
# one beta per transition -- every scale above sample, including the top
# among-units beta -- and gamma). Then alpha * prod(beta_*) == gamma.
.hier_to_matrix <- function(df, q) {
  scales <- levels(df$scale)
  betas <- setdiff(scales, "sample")            # every transition above sample
  cols <- c("alpha", paste0("beta_", betas), "gamma")
  out <- matrix(NA_real_, nrow = length(q), ncol = length(cols),
                dimnames = list(paste0("q", q), cols))
  for (qi in seq_along(q)) {
    rows <- df[df$q == q[qi], ]
    out[qi, "alpha"] <- rows$diversity[rows$scale == "sample"]
    out[qi, "gamma"] <- rows$diversity[rows$scale == "total"]
    for (lv in betas) {
      out[qi, paste0("beta_", lv)] <- rows$beta[rows$scale == lv]
    }
  }
  out
}
