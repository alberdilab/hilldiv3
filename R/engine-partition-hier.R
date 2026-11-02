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
#' global sample weights \eqn{p_j} (equal by default, \eqn{1/n}). Each type
#' supplies (i) a per-feature, per-sample contribution \eqn{c_{ij}} (features
#' are taxa for neutral/functional, branches for phylogenetic), (ii) a
#' per-feature *measure* \eqn{w_i} (1 for neutral, branch length \eqn{L_i} for
#' phylogenetic, attribute contribution \eqn{v_i} for functional), and (iii) a
#' grand-total normaliser \eqn{C} (the number of samples \eqn{n} for neutral,
#' the summed tree depth \eqn{T_+ = \sum_j T_j} for phylogenetic, the total
#' count \eqn{n_+} for functional under equal weights). At an aggregation scale
#' \eqn{k} whose units are \eqn{u} (each pooling samples with pooled
#' contribution \eqn{m_{u,i} = \sum_{j \in u} c_{ij}}), let
#' \eqn{p_u = \sum_{j \in u} p_j} and \eqn{x_{u,i} = m_{u,i}/(p_u C)}. The
#' linear-weighted power sum and diversity are, for \eqn{q \ne 1},
#'
#' \deqn{A_k = \left(\sum_u p_u \sum_i w_i x_{u,i}^q\right)^{1/(1-q)}.}
#'
#' For \eqn{q = 1}, \eqn{A_k = \exp(\sum_u p_u H_u)} with
#' \eqn{H_u = -\sum_i w_i x_{u,i}\log x_{u,i}}. This linear unit weighting
#' gives a Routledge/Tsallis-type decomposition: coarsening averages the
#' \eqn{x_{u,i}} with weights \eqn{p_u}, so convexity (or concavity when
#' \eqn{q < 1}) gives each \eqn{\beta_k = A_k/A_{k-1} \geq 1} for
#' \eqn{q \geq 0}. The finest scale (each sample its own unit) gives
#' \eqn{A_0 = } the single-level alpha of [hill_partition()] under equal sample
#' weights; the coarsest scale (one pooled unit) gives \eqn{A_K = } gamma;
#' intermediate scales interpolate,
#' yielding the nested betas. Because \eqn{C} and \eqn{w_i} are fixed once over
#' the whole table (one shared tree-depth normaliser; one shared \eqn{\tau}
#' and one shared \eqn{v_i}), the chain telescopes for every type.
#'
#' For a transition from child units \eqn{v} to parent units \eqn{u}, the
#' optional overlap ceiling is the beta of a disjoint-feature construction:
#' \eqn{N_k = [\sum_u p_u^{1-q}\sum_{v\in u}p_v^q S_v /
#' \sum_v p_v S_v]^{1/(1-q)}} for \eqn{q\ne 1}, where
#' \eqn{S_v=\sum_i w_i x_{v,i}^q}. At \eqn{q=1}, its continuous limit is
#' \eqn{N_k=\exp[\sum_v p_v M_v\log(p_{u(v)}/p_v)]}, where
#' \eqn{M_v=\sum_i w_i x_{v,i}}. The measure mass \eqn{M_v} matters for
#' functional diversity and non-ultrametric trees.
#'
#' @section Limitations:
#' * Unequal sample weighting uses \eqn{p_u = \sum_{j \in u} p_j}; contributions
#'   are scaled by \eqn{p_j}. The normaliser is \eqn{1} for neutral,
#'   \eqn{\sum_j p_j T_j} for phylogenetic, and \eqn{1} for functional after
#'   computing \eqn{v_i} from the weighted pooled proportions.
#' * **Phylogenetic.** Alpha and gamma are returned in effective-branch-length
#'   (PD) units and share a single tree-depth normaliser across *all* scales, so
#'   the betas telescope. An explicit `depth` truncates an ultrametric tree to
#'   the last `depth` units before its tips. With equal sample weights and an
#'   ultrametric tree,
#'   the finest alpha equals the single-level [hillpart()] phylogenetic alpha
#'   exactly; for non-ultrametric trees the shared normaliser makes per-scale
#'   PD values depth-pooled rather than per-sample.
#' * **Functional.** Each `tau` is capped once over the whole table and the
#'   attribute contributions `v_i` are computed once from the fully pooled
#'   data, so they are constant across scales (a prerequisite for telescoping).
#'   A per-scale `tau` is not supported.
#'
#' @param prep List from [hier_part_prep()] with `contrib` (features x samples),
#'   `measure` (length-features vector \eqn{w_i}), `C` (normaliser) and optional
#'   `sample_weights` (normalised sample weights; `NULL` for equal weights).
#' @param groupings Ordered list of grouping vectors, finest to coarsest, one
#'   entry per *named* hierarchical level (the implicit per-sample and pooled
#'   scales are added here). Each vector has length `ncol(contrib)`.
#' @param q Numeric vector of diversity orders.
#' @param level_names Character names for the grouping levels, aligned with
#'   `groupings` (finest to coarsest).
#' @param with_ceiling Whether to include the data-dependent transition ceiling
#'   `n_eff` used by the overlap metrics.
#'
#' @return A long-format `data.frame` with columns `q`, `scale` (ordered factor,
#'   finest to coarsest), `n_units`, `diversity` (\eqn{A_k}) and `beta`
#'   (\eqn{A_k/A_{k-1}}; `NA` at the finest scale), plus `n_eff` when requested.
#' @keywords internal
#' @noRd
hier_partition <- function(prep, groupings, q, level_names,
                           with_ceiling = FALSE) {
  contrib <- as.matrix(prep$contrib)
  w <- prep$measure
  C <- prep$C
  n <- ncol(contrib)
  sample_weights <- prep$sample_weights

  # Full ordered list of scales, finest -> coarsest: each sample on its own,
  # then the named grouping levels, then everything pooled into one unit.
  scales <- c(list(seq_len(n)), groupings, list(rep(1L, n)))
  scale_names <- c("sample", level_names, "total")

  # Per-scale, q-independent set-up: pooled unit weights and the contribution
  # per unit, divided by that unit's share of the grand normaliser.
  prep_scales <- lapply(scales, function(g) {
    g <- as.factor(g)
    units <- levels(g)
    m <- vapply(units, function(u) rowSums(contrib[, g == u, drop = FALSE]),
                numeric(nrow(contrib)))            # features x units
    pu <- if (is.null(sample_weights)) tabulate(g) / n else
      vapply(units, function(u) sum(sample_weights[g == u]), numeric(1))
    list(x = sweep(m / C, 2, pu, `/`), pu = pu, group = as.integer(g),
         n_units = length(units))
  })

  n_scales <- length(prep_scales)
  if (with_ceiling) {
    parent_of_child <- lapply(seq.int(2L, n_scales), function(s) {
      child <- prep_scales[[s - 1L]]$group
      parent <- prep_scales[[s]]$group
      vapply(seq_along(prep_scales[[s - 1L]]$pu), function(v) {
        unique(parent[child == v])
      }, integer(1))
    })
  }
  out <- vector("list", length(q) * n_scales)
  row <- 1L
  for (qi in seq_along(q)) {
    qv <- q[qi]
    scores <- lapply(prep_scales, function(p) {
      vapply(seq_along(p$pu), function(u) {
        .measure_power_sum(p$x[, u], w, qv)
      }, numeric(1))
    })
    A <- vapply(seq_len(n_scales), function(s) {
      weighted_score <- sum(prep_scales[[s]]$pu * scores[[s]])
      if (qv == 1) exp(weighted_score) else
        weighted_score^(1 / (1 - qv))
    }, numeric(1))
    beta <- c(NA_real_, A[-1] / A[-n_scales])
    if (with_ceiling) {
      n_eff <- c(NA_real_, vapply(seq.int(2L, n_scales), function(s) {
        child <- prep_scales[[s - 1L]]
        parent <- prep_scales[[s]]
        parent_p <- parent$pu[parent_of_child[[s - 1L]]]
        if (qv == 1) {
          mass <- colSums(sweep(child$x, 1, w, `*`))
          exp(sum(child$pu * mass * log(parent_p / child$pu)))
        } else {
          numerator <- sum(parent_p^(1 - qv) * child$pu^qv *
                             scores[[s - 1L]])
          denominator <- sum(child$pu * scores[[s - 1L]])
          (numerator / denominator)^(1 / (1 - qv))
        }
      }, numeric(1)))
    }
    for (s in seq_len(n_scales)) {
      out[[row]] <- data.frame(
        q = qv, scale = scale_names[s], n_units = prep_scales[[s]]$n_units,
        diversity = A[s], beta = beta[s], stringsAsFactors = FALSE
      )
      if (with_ceiling) out[[row]]$n_eff <- n_eff[s]
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
#' Reduces a count table to the `(contrib, measure, C)` quantities consumed by
#' [hier_partition()], mirroring the per-type maths of [part_prep()] but
#' exposing the shared, table-wide quantities (one tree-depth normaliser; one
#' `tau` and one `v_i`) so the nested chain telescopes.
#'
#' @inheritParams hier_partition
#' @param counts Numeric count table (features-as-taxa x samples), aligned.
#' @param type One of `"neutral"`, `"phylogenetic"`, `"functional"`.
#' @param tree A `phylo` tree (phylogenetic only).
#' @param dist A distance matrix (functional only).
#' @param tau Functional threshold; defaults to `max(dist)`.
#' @param sample_weights Normalised sample weights, or `NULL` for equal weights.
#' @param depth Optional positive reference depth for an ultrametric tree.
#'
#' @return List with `contrib`, `measure`, `C` and optional `sample_weights`.
#' @keywords internal
#' @noRd
hier_part_prep <- function(counts, type = "neutral", tree = NULL, dist = NULL,
                           tau = NULL, sample_weights = NULL, depth = NULL) {
  counts <- as.matrix(counts)
  switch(type,
    neutral = {
      pi <- tss(counts)
      if (is.null(sample_weights)) {
        list(contrib = pi, measure = rep(1, nrow(pi)), C = ncol(pi))
      } else {
        weighted_pi <- sweep(pi, 2, sample_weights, `*`)
        list(contrib = weighted_pi,
             measure = rep(1, nrow(pi)), C = sum(weighted_pi),
             sample_weights = sample_weights)
      }
    },
    phylogenetic = {
      if (is.null(tree)) {
        cli::cli_abort("A {.cls phylo} {.arg tree} is required for phylogenetic
                        partitioning.")
      }
      if (!is.null(depth)) {
        if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) ||
            !is.finite(depth) || depth <= 0) {
          cli::cli_abort("{.arg depth} must be one finite, positive number.")
        }
        if (!ape::is.ultrametric(tree)) {
          cli::cli_abort("{.arg depth} requires an ultrametric {.arg tree}.")
        }
        node_depth <- ape::node.depth.edgelength(tree)
        height <- max(node_depth[seq_along(tree$tip.label)])
        if (depth > height) {
          cli::cli_abort("{.arg depth} must not exceed the tree height ({height}).")
        }
      }
      ba <- branch_abundance(tree, tss(counts))   # Li, ai (edges x samples)
      if (!is.null(depth) && depth < height) {
        # Keep only each edge's overlap with the final depth units before tips.
        parent <- tree$edge[, 1]
        child <- tree$edge[, 2]
        ba$Li <- pmax(0, node_depth[child] - pmax(node_depth[parent],
                                                  height - depth))
      }
      if (is.null(sample_weights)) {
        list(contrib = ba$ai, measure = ba$Li,
             C = sum(ba$Li * rowSums(ba$ai)))       # T_+ = sum_j T_j
      } else {
        ai <- sweep(ba$ai, 2, sample_weights, `*`)
        list(contrib = ai, measure = ba$Li,
             C = sum(ba$Li * rowSums(ai)),          # sum_j p_j T_j
             sample_weights = sample_weights)
      }
    },
    functional = {
      if (is.null(dist)) {
        cli::cli_abort("A {.arg dist} matrix is required for functional
                        partitioning.")
      }
      sim <- .functional_similarity(dist, tau)
      if (is.null(sample_weights)) {
        aik <- sim %*% counts                       # taxa x samples
        vi <- rowSums(counts) / rowSums(aik)        # attribute contribution v_i
        vi[!is.finite(vi)] <- 0                      # absent taxa: 0/0 -> 0
        list(contrib = aik, measure = vi, C = sum(counts))
      } else {
        weighted_counts <- sweep(tss(counts), 2, sample_weights, `*`)
        aik <- sim %*% weighted_counts
        vi <- rowSums(weighted_counts) / rowSums(aik)
        vi[!is.finite(vi)] <- 0
        list(contrib = aik, measure = vi, C = sum(weighted_counts),
             sample_weights = sample_weights)
      }
    },
    cli::cli_abort("Unknown diversity type {.val {type}}.")
  )
}

# Measure-weighted power sum, or Shannon entropy at q = 1. Zero cells are
# dropped so q = 0 and the log limit are well defined.
.measure_power_sum <- function(x, w, qvalue) {
  keep <- x != 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(0)
  if (qvalue == 1) {
    -sum(w * x * log(x))
  } else {
    sum(w * x^qvalue)
  }
}

#' Resolve sample weights for hierarchical partitioning
#'
#' @param weights `"equal"`, `"abundance"`, or a positive numeric vector.
#' @param counts Count table with samples in columns.
#' @return Normalised numeric vector, or `NULL` for equal weights.
#' @keywords internal
#' @noRd
.hier_sample_weights <- function(weights, counts) {
  if (identical(weights, "equal")) return(NULL)
  n <- ncol(counts)
  if (identical(weights, "abundance")) {
    weights <- colSums(counts)
  } else if (!is.numeric(weights)) {
    cli::cli_abort("{.arg weights} must be {.val equal}, {.val abundance}, or a
                    positive numeric vector.")
  }
  if (length(weights) != n || any(!is.finite(weights)) || any(weights <= 0)) {
    cli::cli_abort("{.arg weights} must contain one finite, positive value per
                    sample.")
  }
  if (!is.null(names(weights))) {
    sample_ids <- colnames(counts)
    if (is.null(sample_ids) || anyNA(names(weights)) ||
        anyDuplicated(names(weights)) || !setequal(names(weights), sample_ids)) {
      cli::cli_abort("Named {.arg weights} must match the sample column names
                      exactly.")
    }
    weights <- weights[sample_ids]
  }
  if (all(weights == weights[1])) return(NULL)
  weights <- weights / max(weights)
  unname(weights / sum(weights))
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
