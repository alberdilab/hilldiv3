#' Hill numbers computation
#'
#' Compute neutral, phylogenetic or functional Hill numbers (alpha diversity)
#' from a single sample or a count table. The diversity type is inferred from
#' the inputs: counts only -> neutral; counts + `tree` -> phylogenetic;
#' counts + `dist` -> functional.
#'
#' @param data Counts: a numeric vector (one sample), a matrix/data.frame
#'   (taxa x samples), a `phyloseq` object or a `TreeSummarizedExperiment`.
#' @param q Numeric vector of diversity orders (>= 0). Defaults to
#'   `c(0, 1, 2)` (richness, Shannon, Simpson).
#' @param tree A phylogenetic tree of class `phylo` whose tip labels match the
#'   taxa in `data`.
#' @param dist A functional distance matrix (or `dist`) over the taxa.
#' @param tau Optional functional distance threshold. Defaults to `max(dist)`.
#'
#' @return A matrix of Hill numbers with diversity orders in rows (`q0`, `q1`,
#'   ...) and samples in columns.
#'
#' @references
#' Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
#' based on Hill numbers. Phil. Trans. R. Soc. B, 365, 3599-3609.\cr\cr
#' Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
#' numbers to DNA-based diversity analyses. Mol. Ecol. Resour., 19, 804-817.
#'
#' @seealso [hillpart()], [hilldiss()], [hillprof()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hilldiv(counts)
#' hilldiv(counts, q = c(0, 1, 2))
#' @export
hilldiv <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL) {
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q)
  type <- attr(x, "type")
  cli::cli_inform("Computing {type} Hill numbers of {.val {paste0('q', q)}}.")
  hill_alpha(x$counts, q = q, type = type,
             tree = x$tree, dist = x$dist, tau = tau)
}
