#' Hill-number alpha diversity engine
#'
#' Internal compute engine for per-sample (alpha) Hill numbers. This is the
#' single source of truth for the diversity maths; the user-facing [hilldiv()]
#' wrapper validates/aligns inputs and then calls this. Inputs here are assumed
#' to be already validated and aligned.
#'
#' The `q = 1` limit (`exp(-sum(p log p))`) is defined here once and reused for
#' all diversity types.
#'
#' @param p Numeric matrix of relative abundances (taxa x samples). Columns need
#'   not be normalised; they are normalised internally via [tss()].
#' @param q Numeric vector of diversity orders (>= 0).
#' @param type One of "neutral", "phylogenetic", "functional".
#' @param tree A `phylo` tree (required when `type = "phylogenetic"`).
#' @param dist A distance matrix (required when `type = "functional"`).
#' @param tau Functional threshold; defaults to `max(dist)`.
#'
#' @return Numeric matrix of diversity values (q in rows, samples in columns).
#' @keywords internal
#' @noRd
hill_alpha <- function(p, q = c(0, 1, 2), type = "neutral",
                       tree = NULL, dist = NULL, tau = NULL) {
  p <- tss(as.matrix(p))
  switch(type,
    neutral = hill_alpha_neutral(p, q),
    phylogenetic = hill_alpha_phylo(p, q, tree),
    functional = hill_alpha_func(p, q, dist, tau),
    cli::cli_abort("Unknown diversity type {.val {type}}.")
  )
}

# Generic effective-number transform of a set of proportions for a given q.
# `props` must sum to ~1 over its nonzero entries.
.hill_from_props <- function(props, qvalue) {
  props <- props[props != 0]
  if (length(props) == 0) return(0)
  if (qvalue == 1) {
    exp(-sum(props * log(props)))
  } else {
    sum(props^qvalue)^(1 / (1 - qvalue))
  }
}

hill_alpha_neutral <- function(p, q) {
  res <- vapply(q, function(qv) apply(p, 2, .hill_from_props, qvalue = qv),
                numeric(ncol(p)))
  # vapply drops to a vector with a single sample; force samples x q.
  res <- matrix(res, nrow = ncol(p), ncol = length(q))
  .shape_alpha(t(res), q, colnames(p))
}

hill_alpha_phylo <- function(p, q, tree) {
  ba <- branch_abundance(tree, p)
  Li <- ba$Li
  ai <- ba$ai                       # edges x samples
  Tbar <- colSums(Li * ai)          # per-sample tree depth T
  present <- colSums(p != 0)

  res <- matrix(0, nrow = length(q), ncol = ncol(p))
  for (r in seq_along(q)) {
    qv <- q[r]
    res[r, ] <- vapply(seq_len(ncol(p)), function(j) {
      if (present[j] == 0) return(0)
      if (present[j] == 1) return(1)
      a <- ai[, j]
      Tj <- Tbar[j]
      props <- (Li * a / Tj)
      .hill_from_props(props, qv) / Tj
    }, numeric(1))
  }
  .shape_alpha(res, q, colnames(p))
}

hill_alpha_func <- function(p, q, dist, tau) {
  dij <- as.matrix(dist)
  if (is.null(tau)) tau <- max(dij)
  dij[dij > tau] <- tau
  sim <- 1 - dij / tau

  res <- matrix(0, nrow = length(q), ncol = ncol(p))
  for (j in seq_len(ncol(p))) {
    vec <- p[, j]
    a <- as.vector(sim %*% vec)
    keep <- a != 0
    a <- a[keep]
    v <- vec[keep] / a
    for (r in seq_along(q)) {
      qv <- q[r]
      if (qv == 1) {
        res[r, j] <- exp(sum(-v * a * log(a)))
      } else {
        res[r, j] <- sum(v * a^qv)^(1 / (1 - qv))
      }
    }
  }
  .shape_alpha(res, q, colnames(p))
}

# Attach dimnames to an alpha result matrix.
.shape_alpha <- function(res, q, samples) {
  res <- as.matrix(res)
  rownames(res) <- paste0("q", q)
  colnames(res) <- samples
  res
}
