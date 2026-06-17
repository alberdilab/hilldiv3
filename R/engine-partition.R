#' Hill-number diversity partitioning engine
#'
#' Internal compute engine for alpha / gamma / beta partitioning. Single source
#' of truth used by [hillpart()], [hilldiss()] and [hillsim()]. Inputs are
#' assumed validated and aligned.
#'
#' For phylogenetic diversity a single reference `T` (computed from an even pool
#' of all present tips) is used across samples, matching the
#' diversity-partitioning approach in hilldiv2's `hilldiv.phylogenetic.multi`.
#'
#' @inheritParams hill_alpha
#'
#' @return Numeric matrix with columns `alpha`, `gamma`, `beta` (q in rows).
#' @keywords internal
#' @noRd
hill_partition <- function(p, q = c(0, 1, 2), type = "neutral",
                           tree = NULL, dist = NULL, tau = NULL) {
  p <- as.matrix(p)
  # Drop all-zero taxa (rows) so empty lineages don't distort partitioning.
  p <- p[rowSums(p) > 0, , drop = FALSE]
  switch(type,
    neutral = hill_part_neutral(p, q),
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
  # TODO(engine): port and unify the phylogenetic partitioning from hilldiv2
  # (hillpart.phylogenetic) using branch_abundance(); use a single reference T.
  .NotYetImplemented()
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
