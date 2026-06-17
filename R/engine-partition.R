#' Hill-number diversity partitioning engine
#'
#' Internal compute engine for alpha / gamma / beta partitioning. Single source
#' of truth used by [hillpart()], [hilldiss()] and [hillsim()]. Inputs are
#' assumed validated and aligned.
#'
#' Phylogenetic partitioning follows Chiu, Jost & Chao (2014). A single mean
#' tree depth `T = sum_j T_j` (the sum of per-sample tree depths) is shared by
#' alpha and gamma so that `beta = gamma / alpha` lies in `[1, N]`. Alpha and
#' gamma are returned in PD units (effective branch length); dividing both by
#' `T / N` would give effective-species units but leaves beta unchanged. This
#' matches hilldiv2's `hillpart.phylogenetic`. Note the shared `T` differs from
#' the per-sample `T` used by [hill_alpha()], so partition alpha is not the mean
#' of [hilldiv()] phylogenetic values.
#'
#' @inheritParams hill_alpha
#'
#' @return Numeric matrix with columns `alpha`, `gamma`, `beta` (q in rows).
#' @keywords internal
#' @noRd
hill_partition <- function(p, q = c(0, 1, 2), type = "neutral",
                           tree = NULL, dist = NULL, tau = NULL) {
  p <- as.matrix(p)
  switch(type,
    # Drop all-zero taxa for the neutral path; for phylogenetic/functional the
    # full taxon set must be kept to stay aligned with the tree / distance
    # matrix (empty lineages are excluded internally instead).
    neutral = hill_part_neutral(p[rowSums(p) > 0, , drop = FALSE], q),
    phylogenetic = hill_part_phylo(p, q, tree),
    functional = hill_part_func(p, q, dist, tau),
    cli::cli_abort("Unknown diversity type {.val {type}}.")
  )
}

hill_part_neutral <- function(p, q) {
  N <- ncol(p)
  pi <- tss(p)
  gamma_props <- rowSums(pi / N)
  out <- .part_matrix(q)
  for (r in seq_along(q)) {
    qv <- q[r]
    alpha <- (1 / N) * .hill_from_props(as.vector(pi / N), qv)
    gamma <- .hill_from_props(gamma_props, qv)
    out[r, ] <- c(alpha, gamma, gamma / alpha)
  }
  out
}

hill_part_phylo <- function(p, q, tree) {
  if (is.null(tree)) {
    cli::cli_abort("A {.cls phylo} {.arg tree} is required for phylogenetic
                    partitioning.")
  }
  N <- ncol(p)
  pi <- tss(p)                          # normalise each sample to sum 1
  ba <- branch_abundance(tree, pi)
  Li <- ba$Li                           # branch lengths (edges)
  aij <- ba$ai                          # edges x samples (a_ik)
  ai <- rowSums(aij)                    # pooled per-edge abundance (a_i+)
  Tval <- sum(Li * ai)                  # shared T = sum_j T_j (Chiu et al. 2014)

  # Gamma uses pooled branches; alpha uses every (edge, sample) cell. Restrict
  # to nonzero entries so q = 0 and the q = 1 log limit stay well defined.
  keep_g <- ai != 0
  Lg <- Li[keep_g]
  ag <- ai[keep_g] / Tval
  Lmat <- matrix(Li, nrow = nrow(aij), ncol = N)
  keep_a <- aij != 0
  La <- Lmat[keep_a]
  aa <- aij[keep_a] / Tval

  out <- .part_matrix(q)
  for (r in seq_along(q)) {
    qv <- q[r]
    if (qv == 1) {
      alpha <- (1 / N) * exp(-sum(La * aa * log(aa)))
      gamma <- exp(-sum(Lg * ag * log(ag)))
    } else {
      alpha <- (1 / N) * sum(La * aa^qv)^(1 / (1 - qv))
      gamma <- sum(Lg * ag^qv)^(1 / (1 - qv))
    }
    out[r, ] <- c(alpha, gamma, gamma / alpha)
  }
  out
}

hill_part_func <- function(p, q, dist, tau) {
  # TODO(engine): port functional partitioning from hilldiv2 (hillpart.functional).
  .NotYetImplemented()
}

# Empty alpha/gamma/beta result matrix.
.part_matrix <- function(q) {
  out <- matrix(0, nrow = length(q), ncol = 3)
  rownames(out) <- paste0("q", q)
  colnames(out) <- c("alpha", "gamma", "beta")
  out
}
