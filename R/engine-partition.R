#' Hill-number diversity partitioning engine
#'
#' Internal compute engine for alpha / gamma / beta partitioning. Single source
#' of truth used by [hillpart()], [hilldiss()], [hillsim()] and [hillpair()].
#' Inputs are assumed validated and aligned.
#'
#' The engine is split into two steps so callers that partition many subsets of
#' the same samples (notably [hillpair()], which partitions every pair) pay the
#' per-type set-up cost once:
#' * `part_prep()` does the type-specific, sample-independent work -- per-sample
#'   normalisation, the `ape` tree traversal, or the functional similarity
#'   product -- over the whole table.
#' * `part_eval()` combines a chosen subset of columns into alpha/gamma/beta.
#'
#' Phylogenetic partitioning follows Chiu, Jost & Chao (2014). A single mean
#' tree depth `T = sum_j T_j` (the sum of per-sample tree depths over the chosen
#' subset) is shared by alpha and gamma so that `beta = gamma / alpha` lies in
#' `[1, N]`. Alpha and gamma are returned in PD units (effective branch length).
#' This matches hilldiv2's `hillpart.phylogenetic`. Note the shared `T` differs
#' from the per-sample `T` used by [hill_alpha()], so partition alpha is not the
#' mean of [hilldiv()] phylogenetic values.
#'
#' Functional partitioning follows Chiu & Chao (2014): the distance matrix is
#' capped at `tau`, turned into an attribute similarity `1 - d/tau`, and the
#' resulting effective abundances are pooled across the subset. Unlike the
#' neutral and phylogenetic paths it works on the raw counts normalised only by
#' the grand total (no per-sample `tss`), matching hilldiv2's
#' `hillpart.functional`.
#'
#' @inheritParams hill_alpha
#'
#' @return Numeric matrix with columns `alpha`, `gamma`, `beta` (q in rows).
#' @keywords internal
#' @noRd
hill_partition <- function(p, q = c(0, 1, 2), type = "neutral",
                           tree = NULL, dist = NULL, tau = NULL) {
  p <- as.matrix(p)
  prep <- part_prep(p, type, tree = tree, dist = dist, tau = tau)
  part_eval(prep, seq_len(ncol(p)), q)
}

# Sample-independent set-up, done once over the full count table.
part_prep <- function(p, type = "neutral", tree = NULL, dist = NULL,
                      tau = NULL) {
  p <- as.matrix(p)
  switch(type,
    neutral = list(type = "neutral", pi = tss(p)),
    phylogenetic = {
      if (is.null(tree)) {
        cli::cli_abort("A {.cls phylo} {.arg tree} is required for phylogenetic
                        partitioning.")
      }
      ba <- branch_abundance(tree, tss(p))
      list(type = "phylogenetic", Li = ba$Li, aij = ba$ai)
    },
    functional = {
      if (is.null(dist)) {
        cli::cli_abort("A {.arg dist} matrix is required for functional
                        partitioning.")
      }
      dij <- as.matrix(dist)
      if (is.null(tau)) tau <- max(dij)
      dij[dij > tau] <- tau
      list(type = "functional", p = p, aik = (1 - dij / tau) %*% p)
    },
    cli::cli_abort("Unknown diversity type {.val {type}}.")
  )
}

# Combine the chosen columns (sample indices) into an alpha/gamma/beta matrix.
part_eval <- function(prep, cols, q) {
  switch(prep$type,
    neutral = .part_eval_neutral(prep, cols, q),
    phylogenetic = .part_eval_phylo(prep, cols, q),
    functional = .part_eval_func(prep, cols, q)
  )
}

.part_eval_neutral <- function(prep, cols, q) {
  pi <- prep$pi[, cols, drop = FALSE]
  N <- length(cols)
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

.part_eval_phylo <- function(prep, cols, q) {
  Li <- prep$Li
  aij <- prep$aij[, cols, drop = FALSE]      # edges x subset
  N <- length(cols)
  ai <- rowSums(aij)                          # pooled per-edge abundance (a_i+)
  Tval <- sum(Li * ai)                        # shared T = sum_j T_j

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

.part_eval_func <- function(prep, cols, q) {
  p <- prep$p[, cols, drop = FALSE]
  aik <- prep$aik[, cols, drop = FALSE]       # taxa x subset (a_ik)
  N <- length(cols)
  aiplus <- rowSums(aik)                       # pooled per-taxon (a_i+)
  vi <- rowSums(p) / aiplus                    # attribute contributions v_i
  nplus <- sum(p)

  # Gamma uses pooled taxa; alpha uses every (taxon, sample) cell. Restrict to
  # nonzero entries so q = 0 and the q = 1 log limit stay well defined.
  keep_g <- aiplus != 0
  vg <- vi[keep_g]
  ag <- aiplus[keep_g] / nplus
  vmat <- matrix(vi, nrow = nrow(aik), ncol = N)
  keep_a <- aik != 0
  va <- vmat[keep_a]
  aa <- aik[keep_a] / nplus

  out <- .part_matrix(q)
  for (r in seq_along(q)) {
    qv <- q[r]
    if (qv == 1) {
      alpha <- (1 / N) * exp(-sum(va * aa * log(aa)))
      gamma <- exp(-sum(vg * ag * log(ag)))
    } else {
      alpha <- (1 / N) * sum(va * aa^qv)^(1 / (1 - qv))
      gamma <- sum(vg * ag^qv)^(1 / (1 - qv))
    }
    out[r, ] <- c(alpha, gamma, gamma / alpha)
  }
  out
}

# Empty alpha/gamma/beta result matrix.
.part_matrix <- function(q) {
  out <- matrix(0, nrow = length(q), ncol = 3)
  rownames(out) <- paste0("q", q)
  colnames(out) <- c("alpha", "gamma", "beta")
  out
}
