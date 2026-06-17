#' Convert a trait table into a distance matrix
#'
#' Build a pairwise functional distance matrix from a table of taxon traits,
#' suitable as the `dist` argument of [hilldiv()] and friends.
#'
#' @param traits A table with taxa (OTUs/ASVs/MAGs) in rows and traits in
#'   columns. Traits may be continuous, binary or proportional.
#' @param method Distance metric passed to [cluster::daisy()]: `"gower"`
#'   (default), `"euclidean"` or `"manhattan"`.
#'
#' @return A numeric distance matrix.
#' @examples
#' traits <- data.frame(body = c(1, 0.2, 0.9), diet = c(0L, 1L, 1L),
#'                      row.names = c("t1", "t2", "t3"))
#' traits2dist(traits)
#' @export
traits2dist <- function(traits, method = c("gower", "euclidean", "manhattan")) {
  rlang::check_installed("cluster")
  method <- match.arg(method)
  traits <- as.data.frame(traits)
  # Drop constant columns, which carry no distance information.
  non_constant <- vapply(traits, function(col) length(unique(col)) > 1, logical(1))
  traits <- traits[, non_constant, drop = FALSE]
  d <- cluster::daisy(traits, metric = method, warnType = FALSE)
  as.matrix(d)
}
