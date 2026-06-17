#' Multi-scale (hierarchical) Hill-number partitioning engine
#'
#' Internal compute engine for nested, multi-level diversity partitioning of
#' neutral, phylogenetic or functional Hill numbers. Given a sampling hierarchy
#' (e.g. samples nested in sites nested in regions), it returns one beta per
#' hierarchical transition such that the chain telescopes exactly:
#'
#' \deqn{{}^qD_\gamma = {}^qD_{\alpha} \cdot \prod_k \beta_k}
#'
#' i.e. gamma equals the finest alpha times the product of the per-level betas.
#'
#' **Unified construction.** All three diversity types share one formula under
#' global equal sample weights (every sample weighted \eqn{1/n}). Each type
#' supplies (i) a per-feature, per-sample contribution \eqn{c_{ij}} (features are
#' taxa for neutral/functional, branches for phylogenetic), (ii) a per-feature
#' *measure* \eqn{w_i} (1 for neutral, branch length \eqn{L_i} for phylogenetic,
#' attribute contribution \eqn{v_i} for functional), and (iii) a grand-total
#' normaliser \eqn{C} (the number of samples \eqn{n} for neutral, the summed
#' tree depth \eqn{T_+ = \sum_j T_j} for phylogenetic, the total count
#' \eqn{n_+} for functional). At an aggregation scale \eqn{k} whose units are
#' \eqn{u} (each pooling \eqn{n_u} of the \eqn{n} samples, with pooled
#' contribution \eqn{m_{u,i} = \sum_{j \in u} c_{ij}}), the diversity is
#'
#' \deqn{A_k = \frac{\left(\sum_{u,i} w_i (m_{u,i}/C)^q\right)^{1/(1-q)}}
#'                  {\left(\sum_{u} (n_u/n)^q\right)^{1/(1-q)}}.}
#'
#' The numerator is the (measure-weighted) Hill number of the joint
#' feature-by-unit distribution; the denominator is the numbers-equivalent of
#' the unit weights, i.e. Jost's (2007) weighted alpha. Hence each
#' \eqn{\beta_k = A_k/A_{k-1}} is a proper Hill beta, \eqn{\geq 1} and
#' independent of \eqn{A}. The finest scale (each sample its own unit) gives
#' \eqn{A_0 = } the single-level alpha of [hill_partition()]; the coarsest scale
#' (one pooled unit) gives \eqn{A_K = } gamma; intermediate scales interpolate,
#' yielding the nested betas. Because \eqn{C} and \eqn{w_i} are fixed once over
#' the whole table (one shared \eqn{T_+}; one shared \eqn{\tau} and one shared
#' \eqn{v_i}), the chain telescopes for every type.
#'
#' @section Limitations:
#' * Sample weighting is fixed to equal per sample (`n_u / n`). This is what
#'   keeps every beta independent of alpha for all `q` (Jost 2007); abundance- or
#'   effort-weighting is not currently offered.
#' * **Phylogenetic.** Alpha and gamma are returned in effective-branch-length
#'   (PD) units and share a single mean tree depth `T_+` across *all* scales, so
#'   the betas telescope. For ultrametric trees the finest alpha equals the
#'   single-level [hillpart()] phylogenetic alpha exactly; for non-ultrametric
#'   trees it still telescopes but the shared `T_+` makes per-scale PD values
#'   depth-pooled rather than per-sample.
#' * **Functional.** `tau` is capped once over the whole table and the attribute
#'   contributions `v_i` are computed once from the fully pooled data, so they
#'   are constant across scales (a prerequisite for telescoping). Supplying a
#'   per-subset `tau` is therefore not meaningful here.
#'
#' @param prep List from [hier_part_prep()] with `contrib` (features x samples),
#'   `measure` (length-features vector \eqn{w_i}) and `C` (normaliser).
#' @param groupings Ordered list of grouping vectors, finest to coarsest, one
#'   entry per *named* hierarchical level (the implicit per-sample and pooled
#'   scales are added here). Each vector has length `ncol(contrib)`.
#' @param q Numeric vector of diversity orders.
#' @param level_names Character names for the grouping levels, aligned with
#'   `groupings` (finest to coarsest).
#'
#' @return A long-format `data.frame` with columns `q`, `scale` (ordered factor,
#'   finest to coarsest), `n_units`, `diversity` (\eqn{A_k}) and `beta`
#'   (\eqn{A_k/A_{k-1}}; `NA` at the finest scale).
#' @keywords internal
#' @noRd
hier_partition <- function(prep, groupings, q, level_names) {
  contrib <- as.matrix(prep$contrib)
  w <- prep$measure
  C <- prep$C
  n <- ncol(contrib)

  # Full ordered list of scales, finest -> coarsest: each sample on its own,
  # then the named grouping levels, then everything pooled into one unit.
  scales <- c(list(seq_len(n)), groupings, list(rep(1L, n)))
  scale_names <- c("sample", level_names, "total")

  # Per-scale, q-independent set-up: the joint feature-by-unit proportions z and
  # the matching per-feature measure, plus the unit weights n_u / n.
  prep_scales <- lapply(scales, function(g) {
    g <- as.factor(g)
    units <- levels(g)
    m <- vapply(units, function(u) rowSums(contrib[, g == u, drop = FALSE]),
                numeric(nrow(contrib)))            # features x units
    list(z = as.vector(m / C),
         w = rep(w, times = length(units)),
         nu = tabulate(g) / n,
         n_units = length(units))
  })

  n_scales <- length(prep_scales)
  out <- vector("list", length(q) * n_scales)
  row <- 1L
  for (qi in seq_along(q)) {
    qv <- q[qi]
    A <- vapply(prep_scales, function(p) {
      .hill_from_measure(p$z, p$w, qv) / .hill_from_props(p$nu, qv)
    }, numeric(1))
    beta <- c(NA_real_, A[-1] / A[-n_scales])
    for (s in seq_len(n_scales)) {
      out[[row]] <- data.frame(
        q = qv, scale = scale_names[s], n_units = prep_scales[[s]]$n_units,
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

#' Type-specific set-up for hierarchical partitioning
#'
#' Reduces a count table to the `(contrib, measure, C)` triple consumed by
#' [hier_partition()], mirroring the per-type maths of [part_prep()] but exposing
#' the shared, table-wide quantities (one `T_+`; one `tau` and one `v_i`) so the
#' nested chain telescopes.
#'
#' @inheritParams hier_partition
#' @param counts Numeric count table (features-as-taxa x samples), aligned.
#' @param type One of `"neutral"`, `"phylogenetic"`, `"functional"`.
#' @param tree A `phylo` tree (phylogenetic only).
#' @param dist A distance matrix (functional only).
#' @param tau Functional threshold; defaults to `max(dist)`.
#'
#' @return List with `contrib`, `measure`, `C`.
#' @keywords internal
#' @noRd
hier_part_prep <- function(counts, type = "neutral", tree = NULL, dist = NULL,
                           tau = NULL) {
  counts <- as.matrix(counts)
  switch(type,
    neutral = {
      pi <- tss(counts)
      list(contrib = pi, measure = rep(1, nrow(pi)), C = ncol(pi))
    },
    phylogenetic = {
      if (is.null(tree)) {
        cli::cli_abort("A {.cls phylo} {.arg tree} is required for phylogenetic
                        partitioning.")
      }
      ba <- branch_abundance(tree, tss(counts))    # Li (edges), ai (edges x samp)
      list(contrib = ba$ai, measure = ba$Li,
           C = sum(ba$Li * rowSums(ba$ai)))         # T_+ = sum_j T_j
    },
    functional = {
      if (is.null(dist)) {
        cli::cli_abort("A {.arg dist} matrix is required for functional
                        partitioning.")
      }
      dij <- as.matrix(dist)
      if (is.null(tau)) tau <- max(dij)
      dij[dij > tau] <- tau
      aik <- (1 - dij / tau) %*% counts             # taxa x samples
      vi <- rowSums(counts) / rowSums(aik)          # attribute contribution v_i
      vi[!is.finite(vi)] <- 0                        # absent taxa: 0/0 -> 0
      list(contrib = aik, measure = vi, C = sum(counts))
    },
    cli::cli_abort("Unknown diversity type {.val {type}}.")
  )
}

# Measure-weighted effective-number transform: (sum_i w_i z_i^q)^(1/(1-q)), with
# the q = 1 limit exp(-sum w_i z_i log z_i). `z` are the joint proportions (they
# sum to 1 over nonzero entries), `w` the per-feature measure (1 for neutral,
# branch length for phylogenetic, attribute contribution for functional). With
# w == 1 this is exactly `.hill_from_props()`. Zero-abundance cells are dropped
# so q = 0 and the q = 1 log limit stay well defined.
.hill_from_measure <- function(z, w, qvalue) {
  keep <- z != 0
  z <- z[keep]
  w <- w[keep]
  if (length(z) == 0) return(0)
  if (qvalue == 1) {
    exp(-sum(w * z * log(z)))
  } else {
    sum(w * z^qvalue)^(1 / (1 - qvalue))
  }
}

#' Resolve a hierarchy formula into ordered nested grouping vectors
#'
#' Turns a nesting formula such as `~ region / site` into the ordered list of
#' grouping vectors (finest to coarsest) consumed by [hier_partition()].
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
