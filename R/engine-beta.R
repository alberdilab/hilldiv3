#' Convert beta diversity to (dis)similarity metrics
#'
#' Internal single source of truth for the four Chiu et al. (2014) overlap /
#' turnover (dis)similarity metrics derived from a beta value. hilldiv2 defined
#' these twice (once in `hilldiss`, once as complements in `hillsim`); here they
#' are defined once and the similarity is simply `1 - dissimilarity`.
#'
#' The `q -> 1` cases use the analytic limit rather than the `0.999999`
#' numerical approximation used in hilldiv2.
#'
#' @param beta Numeric beta value.
#' @param N Number of samples (assemblages).
#' @param q Diversity order corresponding to `beta`.
#'
#' @return Named numeric vector with elements `S`, `C`, `U`, `V`
#'   (dissimilarities).
#' @keywords internal
#' @noRd
beta_to_dissim <- function(beta, N, q) {
  invb <- 1 / beta

  # S: Jaccard-type turnover.
  S <- 1 - ((invb - 1 / N) / (1 - 1 / N))
  # V: Sorensen-type turnover.
  V <- 1 - ((N - beta) / (N - 1))

  # C and U are the Sorensen- and Jaccard-type overlap complements. Both
  # converge to the same log-based limit as q -> 1 (verified by
  # L'Hopital on the hilldiv2 expressions): C = U = log(beta) / log(N).
  if (q == 1) {
    C <- log(beta) / log(N)
    U <- C
  } else {
    C <- 1 - ((invb^(q - 1) - (1 / N)^(q - 1)) / (1 - (1 / N)^(q - 1)))
    U <- 1 - ((invb^(1 - q) - (1 / N)^(1 - q)) / (1 - (1 / N)^(1 - q)))
  }

  c(S = S, C = C, U = U, V = V)
}

#' @keywords internal
#' @noRd
beta_to_sim <- function(beta, N, q) {
  1 - beta_to_dissim(beta, N, q)
}
