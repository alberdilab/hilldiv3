#' Hill numbers diversity partitioning
#'
#' Partition neutral, phylogenetic or functional Hill-number diversity into
#' alpha, gamma and beta components across a set of samples.
#'
#' @inheritParams hilldiv
#' @param data A count table (taxa x samples) or a supported object; a single
#'   sample is not meaningful for partitioning.
#'
#' @return A matrix with columns `alpha`, `gamma`, `beta` and diversity orders
#'   in rows.
#'
#' @seealso [hilldiv()], [hilldiss()], [hillsim()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hillpart(counts)
#' @export
hillpart <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL) {
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q)
  type <- attr(x, "type")
  cli::cli_inform("Partitioning {type} Hill numbers of {.val {paste0('q', q)}}.")
  hill_partition(x$counts, q = q, type = type,
                 tree = x$tree, dist = x$dist, tau = tau)
}
